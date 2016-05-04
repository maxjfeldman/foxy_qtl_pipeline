# Script to make timeseries QTL plots
library(qtl)
library(funqtl)
library(ggplot2)

# Read in arguments
args <- commandArgs(TRUE)

# Directory program launched from
directory<-getwd()
# Base directory of trait
job <- args[1]
# Get treatments
t1 <- args[2]
t2 <- args[3]

# Get directory path and names of tables
dir.job<-paste(directory, job, sep="/")
dir.job.t1<-paste(dir.job, "/", "timeseries.", t1, sep="")
dir.job.t2<-paste(dir.job, "/", "timeseries.", t2, sep="")
# Lets make a set of directories to store the data
dir.create(dir.job.t1)
dir.create(dir.job.t2)

# Make table name
st.name<-paste(job, 'concatenated_summary_table.csv', sep="_")
st.path<-paste(dir.job, st.name, sep="/")
st<-read.csv(st.path, header=F)
colnames(st)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
st<-st[st$type == 'raw',]
st.t1<-st[st$treatment == t1, ]
st.t2<-st[st$treatment == t2, ]


###### This is a function to get unique qtl from a qtl summary table
remove_dup_qtl<-function(temp){
  all_qtl<-sort(table(temp$marker), decreasing=T)
  if (length(all_qtl) == 1) {
    m.names<<-c(m.names, names(all_qtl)[1])
    # <<- means change the global variable (chr<<-max) changes the global variable chr to local variable max
    chr<<-c(chr,unique(temp[temp$marker == names(all_qtl)[1],'chr']))
    pos<<-c(pos,unique(temp[temp$marker == names(all_qtl)[1],'pos']))
    print(chr) 
    print(pos)
  }
  if (length(all_qtl) > 1) {
    name<-names(all_qtl)[1]
    ave.pos<-mean(temp[temp$marker == name, 'pos'])
    m.names<<-c(m.names, names(all_qtl)[1])
    chr<<-c(chr,unique(temp[temp$marker == names(all_qtl)[1],'chr']))
    pos<<-c(pos,unique(temp[temp$marker == names(all_qtl)[1],'pos']))
    max.pos<-ave.pos+10
    min.pos<-ave.pos-10
    temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
    print(ave.pos) 
    print(chr) 
    print(pos)
    remove_dup_qtl(temp)
  }
}

###### Lets use that function to get a list of unique qtl found in t1
chrs<-sort(unique(st.t1$chr))
m.names<-c()
chr<-c()
pos<-c()
for(ch in 1:length(chrs)) {
  temp<-st.t1[st.t1$chr == chrs[ch],]
  temp$marker<-as.character(temp$marker)
  remove_dup_qtl(temp)
}

# Combine the marker names, chromosome and positional info into a data.frame
unique_qtl_t1<-as.data.frame(cbind(m.names, chr, pos))

###### Lets use that function to get a list of unique qtl found in t2
chrs<-sort(unique(st.t2$chr))
m.names<-c()
chr<-c()
pos<-c()
for(ch in 1:length(chrs)) {
  temp<-st.t2[st.t2$chr == chrs[ch],]
  temp$marker<-as.character(temp$marker)
  remove_dup_qtl(temp)
}

# Combine the marker names, chromosome and positional info into a data.frame
unique_qtl_t2<-as.data.frame(cbind(m.names, chr, pos))

# Lets load in the cross object for raw phenotypes 
raw_cr.obj_path<-paste(dir.job, 'cross.object.raw.Rdata', sep="/")
load(raw_cr.obj_path)

# Make a regrex search term for each treatment
t1_meta<-paste(t1, '*', sep="")
t2_meta<-paste(t2, '*', sep="")

cols<-grep(t1_meta, phenames(fg.cr.obj))
assign(paste(t1, 'cols', sep="."), cols)
cols<-grep(t2_meta, phenames(fg.cr.obj))
assign(paste(t2, 'cols', sep="."), cols)

fg.cr.obj <- calc.genoprob(fg.cr.obj, step=0)

#################################################################################
# Treatment 1 (t1)
#################################################################################
pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][2])
}


out<- scanone(fg.cr.obj, pheno.col = get(paste(t1, 'cols', sep=".")), method="hk")
eff <- geteffects(fg.cr.obj, pheno.cols=get(paste(t1, 'cols', sep=".")))
setwd(dir.job.t1)

pdf(paste('so_timeseries_qtl_',job, ".", t1, '.pdf', sep=""))
plotlod(out, eff, get(paste(t1, 'cols', sep=".")), gap=15, ylab="Time")
dev.off()

out.F <- scanoneF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), method="hk")
o.perm.F <- scanoneF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), method = "hk", n.perm=1000)
max.perm.F<-summary(o.perm.F, alpha=0.05)
setwd(dir.job.t1)
qtl.table.so<-summary(out.F, perms=o.perm.F, alpha=0.05, pvalues=TRUE)
write.csv(qtl.table.so, file=paste('qtl.table.so.F_', job, ".", t1,'.csv', sep=""), quote=F, row.names=F)

par(mfrow=c(2,1))
pdf(paste('so_timeseries_qtl_slod_and_mlod_', job, ".", t1, '.pdf' ,sep=""))
plot(out.F, main="The SLOD curve", bandcol="gray90")
abline(h=max.perm.F[1], col="red", lty=3)
plot(out.F, lodcolumn=2, main="The MLOD curve", bandcol="gray90")
abline(h=max.perm.F[2], col="red", lty=3)
dev.off()
par(mfrow=c(1,1))

# MQM
# First SLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), max.qtl=9, usec= "slod", method="hk", penalties=c(max.perm.F[1],0,0))

if(out.mqm.F$n.qtl > 1 ) {

chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

# This is on the chopping block

#if (length(chr) == 1) {
#  temp<-summary(out.F)
#  temp<-temp[temp$chr != chr,]
#  temp<-temp[order(temp$slod, decreasing=T),]
#  chr<-c(chr, temp[1,'chr'])
#  pos<-c(pos, temp[1,'pos'])
#}

#if (length(chr) == 0) {
#  temp<-summary(out.F)
#  temp<-temp[order(temp$slod, decreasing=T),]
#  chr<-temp[1:2, 'chr']
#  pos<-temp[1:2, 'pos']
#}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t1, 'slod', sep="_"), out.mqm.F)
assign(paste('chr', t1, 'slod', sep="_"), chr)
assign(paste('pos', t1, 'slod', sep="_"), pos)
assign(paste('qtl', t1, 'slod', sep="_"), qtl)
assign(paste('my.formula', t1, 'slod', sep="_"), my.formula)
#assign(paste('lodmat.F', t1, 'slod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(2,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_slod.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="SLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
#dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod', t1, 'slod', sep="_"), refqtlslod)

#slodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

#for(i in 1:length(slodeff)) {
#  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(slodeff[[1]])
#slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
#colnames(slodeff) <- nam

#n.col<-ncol(slodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t1, ".pdf", sep=""))
# Draw a plot of the intercept
#plot(days, slodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, slodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(slodeff)[i], side=3, line=0.5)
#}

#dev.off()
#par(mfrow=c(1,1))

#assign(paste('slodeff', t1, 'slod', sep="_"), slodeff)

#write.csv(slodeff, file=paste('slodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)

# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)

# Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
intercept.slod<-c()
summary.table.slod<-c()
for(i in get(paste(t1, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlslod, formula=paste('my.formula', t1, 'slod', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'slod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.slod$est$est[[1]]
  intercept.slod<-c(intercept.slod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.slod<-as.data.frame(out.fitqtl.slod$result.full[,2])
  full.mdl.var.slod$prop.variance<-c(100*(full.mdl.var.slod[1,1]/full.mdl.var.slod[3,1]), 100*(full.mdl.var.slod[2,1]/full.mdl.var.slod[3,1]), 100)
  colnames(full.mdl.var.slod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.slod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.slod<-rbind(dropone.mdl.var.slod, full.mdl.var.slod[3,1])
  dropone.mdl.var.slod$prop.variance<-c(100*(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-1),1]/full.mdl.var.slod[3,1]), 100)
  colnames(dropone.mdl.var.slod)<-c("variance", "prop.variance")
  
  rownames(dropone.mdl.var.slod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),1]-sum(dropone.mdl.var.slod[(1:nrow(dropone.mdl.var.slod)-1),1]), dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),2]-sum(dropone.mdl.var.slod[(1:nrow(dropone.mdl.var.slod)-1),2]))
  dropone.mdl.var.slod<-rbind(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-1),], dropone.residual, dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),])
  
  
  rownames(dropone.mdl.var.slod)[(nrow(dropone.mdl.var.slod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.slod<-refqtlslod$n.qtl
  CI_L<-c()
  CI_R<-c()
  for(q in 1:n.QTL.slod){
    lod_int.slod<-lodint(refqtlslod, qtl.index=q)
    lod_int.slod$marker<-rownames(lod_int.slod)
    L<-as.data.frame(lod_int.slod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.slod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  }
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.slod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.slod<-as.data.frame(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.slod$result.drop[,4]
  slod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(slod.markers, marker.prop.var.slod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.slod<-rbind(summary.table.slod, summary.table)
}

assign(paste('summary.table.slod', t1, sep="."), summary.table.slod)
write.csv(summary.table.slod, file=paste('summary.table.slod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)

names(intercept.slod)<-days
assign(paste('intercept.slod', t1, sep="."), intercept.slod)
write.csv(intercept.slod, file=paste('intercept.slod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)


}


if(out.mqm.F$n.qtl == 1){
  
  chr<-out.mqm.F$chr
  pos<-out.mqm.F$pos
  
  #if (length(chr) == 1) {
  #  temp<-summary(out.F)
  #  temp<-temp[temp$chr != chr,]
  #  temp<-temp[order(temp$slod, decreasing=T),]
  #  chr<-c(chr, temp[1,'chr'])
  #  pos<-c(pos, temp[1,'pos'])
  #}
  
  #if (length(chr) == 0) {
  #  temp<-summary(out.F)
  #  temp<-temp[order(temp$slod, decreasing=T),]
  #  chr<-temp[1:2, 'chr']
  #  pos<-temp[1:2, 'pos']
  #}
  
  qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
  Qs<-paste('Q', 1:length(pos), sep="")
  my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
  #lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")
  
  assign(paste('out.mqm.F', t1, 'slod', sep="_"), out.mqm.F)
  assign(paste('chr', t1, 'slod', sep="_"), chr)
  assign(paste('pos', t1, 'slod', sep="_"), pos)
  assign(paste('qtl', t1, 'slod', sep="_"), qtl)
  assign(paste('my.formula', t1, 'slod', sep="_"), my.formula)
  #assign(paste('lodmat.F', t1, 'slod', sep="_"), lodmat.F)
  
  # Make plots of the QTL LOD profile and the LOD profile over time
  #par(mfrow=c(2,1))
  #pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_slod.pdf', sep=""))
  #plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="slod")
  refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
  #plotLodProfile(refqtlslod)
  #dev.off()
  par(mfrow=c(1,1))
  
  assign(paste('refqtlslod', t1, 'slod', sep="_"), refqtlslod)
  #slodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))
  
  #for(i in 1:length(slodeff)) {
  #  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
  #                                 method="hk", get.ests=TRUE,
  #                                 dropone=FALSE))$ests[,1]*c(1,2,2)
  #}
  
  pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')
  
  days<-c()
  for(p in 1:length(pname_list)) {
    days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
  }
  
  #nam <- names(slodeff[[1]])
  #slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
  #colnames(slodeff) <- nam
  
  #n.col<-ncol(slodeff)
  
  # Lets plot the effect size over time
  #par(mfrow=c(1,n.col))
  #pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t1, ".pdf", sep=""))
  # Draw a plot of the intercept
  #plot(days, slodeff[,1], lwd=2, type="l",
  #     xlab="Days after planting",
  #     ylab="Height (cm)", col="red")
  #mtext("baseline curve", side=3, line=0.5)
  # Now add plot of effect size for each QTL
  #for (i in 2:n.col) {
  #  plot(days, slodeff[,i], lwd=2, type="l",
  #       xlab="Days after planting",
  #       ylab="QTL effect (cm)", col="red")
  #  mtext(colnames(slodeff)[i], side=3, line=0.5)
  #}
  
  #dev.off()
  #par(mfrow=c(1,1))
  
  #assign(paste('slodeff', t1, 'slod', sep="_"), slodeff)
  #write.csv(slodeff, file=paste('slodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)
  # get marker names
  m.names<-find.marker(fg.cr.obj, chr, pos)
  
  # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
  summary.table.slod<-c()
  intercept.slod<-c()
  for(i in get(paste(t1, 'cols', sep="."))) {
    out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlslod, formula=paste('my.formula', t1, 'slod', sep="_"), method='hk', get.ests=T)
    assign(paste('out.fitqtl', 'slod', sep="."), out.fitqtl)
    
    # Get value of height for the intercept
    int<-out.fitqtl.slod$est$est[[1]]
    intercept.slod<-c(intercept.slod, int)
    
    # Proportioning of variance for scanone results (full model)
    full.mdl.var.slod<-as.data.frame(out.fitqtl.slod$result.full[,2])
    full.mdl.var.slod$prop.variance<-c(100*(full.mdl.var.slod[1,1]/full.mdl.var.slod[3,1]), 100*(full.mdl.var.slod[2,1]/full.mdl.var.slod[3,1]), 100)
    colnames(full.mdl.var.slod)<-c("variance", "prop.variance")
    
    
    # Get confidence interval (95%) for each QTL:
    n.QTL.slod<-refqtlslod$n.qtl
    CI_L<-c()
    CI_R<-c()
    for(q in 1:n.QTL.slod){
      lod_int.slod<-lodint(refqtlslod, qtl.index=q)
      lod_int.slod$marker<-rownames(lod_int.slod)
      L<-as.data.frame(lod_int.slod[1,c(4,1,2,3)])
      colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
      R<-as.data.frame(lod_int.slod[3,c(4,1,2,3)])
      colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
      CI_L<-rbind(CI_L, L)
      CI_R<-rbind(CI_R, R)
    }
    CI<-cbind(CI_L, CI_R)
    
    # Get proportion of variance for only marker
    marker.prop.var.slod<-full.mdl.var.slod$prop.variance[1]
    # Get the fx size for only marker
    fx_size.a<-out.fitqtl.slod$est$est[[2]][1]
    fx_size.se<-c(0)
    # Build a table that summarizes the results
    lods<-out.fitqtl.slod$lod
    slod.markers<-cbind(chr, pos, lods)
    phe.name<-colnames(fg.cr.obj$pheno)[i]
    summary.table<-c(slod.markers, marker.prop.var.slod, as.numeric(fx_size.a), fx_size.se, CI, phe.name)
    summary.table<-data.frame(matrix(unlist(summary.table), nrow=1, byrow=T),stringsAsFactors=FALSE)
    rownames(summary.table)<-m.names
    colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
    summary.table.slod<-rbind(summary.table.slod, summary.table)
  }
  
  assign(paste('summary.table.slod', t1, sep="."), summary.table.slod)
  write.csv(summary.table.slod, file=paste('summary.table.slod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)
 
   names(intercept.slod)<-days
  assign(paste('intercept.slod', t1, sep="."), intercept.slod)
  write.csv(intercept.slod, file=paste('intercept.slod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)
  
}


if(out.mqm.F$n.qtl < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t1, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t1, '_cross.object.raw.Rdata', sep=""))

# MQM
# Now MLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), max.qtl=9, usec= "mlod", method="hk", penalties=c(max.perm.F[2],0,0))

if(out.mqm.F$n.qtl > 1) {
 
chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

#if (length(chr) == 1) {
#  temp<-summary(out.F)
#  temp<-temp[temp$chr != chr,]
#  temp<-temp[order(temp$mlod, decreasing=T),]
#  chr<-c(chr, temp[1,'chr'])
#  pos<-c(pos, temp[1,'pos'])
#}

#if (length(chr) == 0) {
#  temp<-summary(out.F)
#  temp<-temp[order(temp$mlod, decreasing=T),]
#  chr<-temp[1:2, 'chr']
#  pos<-temp[1:2, 'pos']
#}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t1, 'mlod', sep="_"), out.mqm.F)
assign(paste('chr', t1, 'mlod', sep="_"), chr)
assign(paste('pos', t1, 'mlod', sep="_"), pos)
assign(paste('qtl', t1, 'mlod', sep="_"), qtl)
assign(paste('my.formula', t1, 'mlod', sep="_"), my.formula)
#assign(paste('lodmat.F', t1, 'mlod', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(2,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_mlod.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
refqtlmlod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlmlod)
#dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlmlod', t1, 'mlod', sep="_"), refqtlmlod)


#mlodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

### Problem crops up here

#for(i in 1:length(mlodeff)) {
#  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(mlodeff[[1]])
#mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
#colnames(mlodeff) <- nam

#n.col<-ncol(mlodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_mlod_',job, ".", t1, '.pdf', sep=""))
# Draw a plot of the intercept
#plot(days, mlodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, mlodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(mlodeff)[i], side=3, line=0.5)
#}


#dev.off()
par(mfrow=c(1,1))

#assign(paste('mlodeff', t1, 'mlod', sep="_"), mlodeff)

#write.csv(mlodeff, file=paste('mlodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)

# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)

# Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
intercept.mlod<-c()
summary.table.mlod<-c()
for(i in get(paste(t1, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlmlod, formula=paste('my.formula', t1, 'mlod', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'mlod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.mlod$est$est[[1]]
  intercept.mlod<-c(intercept.mlod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.mlod<-as.data.frame(out.fitqtl.mlod$result.full[,2])
  full.mdl.var.mlod$prop.variance<-c(100*(full.mdl.var.mlod[1,1]/full.mdl.var.mlod[3,1]), 100*(full.mdl.var.mlod[2,1]/full.mdl.var.mlod[3,1]), 100)
  colnames(full.mdl.var.mlod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.mlod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.mlod<-rbind(dropone.mdl.var.mlod, full.mdl.var.mlod[3,1])
  dropone.mdl.var.mlod$prop.variance<-c(100*(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-1),1]/full.mdl.var.mlod[3,1]), 100)
  colnames(dropone.mdl.var.mlod)<-c("variance", "prop.variance")
  
  rownames(dropone.mdl.var.mlod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),1]-sum(dropone.mdl.var.mlod[(1:nrow(dropone.mdl.var.mlod)-1),1]), dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),2]-sum(dropone.mdl.var.mlod[(1:nrow(dropone.mdl.var.mlod)-1),2]))
  dropone.mdl.var.mlod<-rbind(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-1),], dropone.residual, dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),])
  
  
  rownames(dropone.mdl.var.mlod)[(nrow(dropone.mdl.var.mlod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.mlod<-refqtlmlod$n.qtl
  CI_L<-c()
  CI_R<-c()
  for(q in 1:n.QTL.mlod){
    lod_int.mlod<-lodint(refqtlmlod, qtl.index=q)
    lod_int.mlod$marker<-rownames(lod_int.mlod)
    L<-as.data.frame(lod_int.mlod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.mlod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  }
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mlod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.mlod<-as.data.frame(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.mlod$result.drop[,4]
  mlod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(mlod.markers, marker.prop.var.mlod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.mlod<-rbind(summary.table.mlod, summary.table)
}

assign(paste('summary.table.mlod', t1, sep="."), summary.table.mlod)
write.csv(summary.table.mlod, file=paste('summary.table.mlod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)

names(intercept.mlod)<-days
assign(paste('intercept.mlod', t1, sep="."), intercept.mlod)
write.csv(intercept.mlod, file=paste('intercept.mlod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)

}

if(out.mqm.F$n.qtl == 1){
  
  chr<-out.mqm.F$chr
  pos<-out.mqm.F$pos
  
  #if (length(chr) == 1) {
  #  temp<-summary(out.F)
  #  temp<-temp[temp$chr != chr,]
  #  temp<-temp[order(temp$mlod, decreasing=T),]
  #  chr<-c(chr, temp[1,'chr'])
  #  pos<-c(pos, temp[1,'pos'])
  #}
  
  #if (length(chr) == 0) {
  #  temp<-summary(out.F)
  #  temp<-temp[order(temp$mlod, decreasing=T),]
  #  chr<-temp[1:2, 'chr']
  #  pos<-temp[1:2, 'pos']
  #}
  
  qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
  Qs<-paste('Q', 1:length(pos), sep="")
  my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
  #lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")
  
  assign(paste('out.mqm.F', t1, 'mlod', sep="_"), out.mqm.F)
  assign(paste('chr', t1, 'mlod', sep="_"), chr)
  assign(paste('pos', t1, 'mlod', sep="_"), pos)
  assign(paste('qtl', t1, 'mlod', sep="_"), qtl)
  assign(paste('my.formula', t1, 'mlod', sep="_"), my.formula)
  #assign(paste('lodmat.F', t1, 'mlod', sep="_"), lodmat.F)
  
  # Make plots of the QTL LOD profile and the LOD profile over time
  #par(mfrow=c(2,1))
  #pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_mlod.pdf', sep=""))
  #plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
  refqtlmlod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
  #plotLodProfile(refqtlmlod)
  #dev.off()
  par(mfrow=c(1,1))
  
  assign(paste('refqtlmlod', t1, 'mlod', sep="_"), refqtlmlod)
  #mlodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))
  
  #for(i in 1:length(mlodeff)) {
  #  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
  #                                 method="hk", get.ests=TRUE,
  #                                 dropone=FALSE))$ests[,1]*c(1,2,2)
  #}
  
  pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')
  
  days<-c()
  for(p in 1:length(pname_list)) {
    days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
  }
  
  #nam <- names(mlodeff[[1]])
  #mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
  #colnames(mlodeff) <- nam
  
  #n.col<-ncol(mlodeff)
  
  # Lets plot the effect size over time
  #par(mfrow=c(1,n.col))
  #pdf(paste('mqm_timeseries_fx_size_mlod_', job, ".", t1, ".pdf", sep=""))
  # Draw a plot of the intercept
  #plot(days, mlodeff[,1], lwd=2, type="l",
  #     xlab="Days after planting",
  #     ylab="Height (cm)", col="red")
  #mtext("baseline curve", side=3, line=0.5)
  # Now add plot of effect size for each QTL
  #for (i in 2:n.col) {
  #  plot(days, mlodeff[,i], lwd=2, type="l",
  #       xlab="Days after planting",
  #       ylab="QTL effect (cm)", col="red")
  #  mtext(colnames(mlodeff)[i], side=3, line=0.5)
  #}
  
  #dev.off()
  #par(mfrow=c(1,1))
  
  #assign(paste('mlodeff', t1, 'mlod', sep="_"), mlodeff)
  #write.csv(mlodeff, file=paste('mlodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)
  # get marker names
  m.names<-find.marker(fg.cr.obj, chr, pos)
  
  # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
  summary.table.mlod<-c()
  intercept.mlod<-c()
  for(i in get(paste(t1, 'cols', sep="."))) {
    out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlmlod, formula=paste('my.formula', t1, 'mlod', sep="_"), method='hk', get.ests=T)
    assign(paste('out.fitqtl', 'mlod', sep="."), out.fitqtl)
    
    # Get value of height for the intercept
    int<-out.fitqtl.mlod$est$est[[1]]
    intercept.mlod<-c(intercept.mlod, int)
    
    # Proportioning of variance for scanone results (full model)
    full.mdl.var.mlod<-as.data.frame(out.fitqtl.mlod$result.full[,2])
    full.mdl.var.mlod$prop.variance<-c(100*(full.mdl.var.mlod[1,1]/full.mdl.var.mlod[3,1]), 100*(full.mdl.var.mlod[2,1]/full.mdl.var.mlod[3,1]), 100)
    colnames(full.mdl.var.mlod)<-c("variance", "prop.variance")
    
    
    # Get confidence interval (95%) for each QTL:
    n.QTL.mlod<-refqtlmlod$n.qtl
    CI_L<-c()
    CI_R<-c()
    for(q in 1:n.QTL.mlod){
      lod_int.mlod<-lodint(refqtlmlod, qtl.index=q)
      lod_int.mlod$marker<-rownames(lod_int.mlod)
      L<-as.data.frame(lod_int.mlod[1,c(4,1,2,3)])
      colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
      R<-as.data.frame(lod_int.mlod[3,c(4,1,2,3)])
      colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
      CI_L<-rbind(CI_L, L)
      CI_R<-rbind(CI_R, R)
    }
    CI<-cbind(CI_L, CI_R)
    
    # Get proportion of variance for only marker
    marker.prop.var.mlod<-full.mdl.var.mlod$prop.variance[1]
    # Get the fx size for only marker
    fx_size.a<-out.fitqtl.mlod$est$est[[2]][1]
    fx_size.se<-c(0)
    # Build a table that summarizes the results
    lods<-out.fitqtl.mlod$lod
    mlod.markers<-cbind(chr, pos, lods)
    phe.name<-colnames(fg.cr.obj$pheno)[i]
    summary.table<-c(mlod.markers, marker.prop.var.mlod, as.numeric(fx_size.a), fx_size.se, CI, phe.name)
    summary.table<-data.frame(matrix(unlist(summary.table), nrow=1, byrow=T),stringsAsFactors=FALSE)
    rownames(summary.table)<-m.names
    colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
    summary.table.mlod<-rbind(summary.table.mlod, summary.table)
  }
  
  assign(paste('summary.table.mlod', t1, sep="."), summary.table.mlod)
  write.csv(summary.table.mlod, file=paste('summary.table.mlod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)
  
  names(intercept.mlod)<-days
  assign(paste('intercept.mlod', t1, sep="."), intercept.mlod)
  write.csv(intercept.mlod, file=paste('intercept.mlod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)
  
}

if(out.mqm.F$n.qtl < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t1, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}
save.image(file=paste('timeseries_', t1, '_cross.object.raw.Rdata', sep=""))


######## Now lets do the same thing all unique qtl at all time points for t1
unique_qtl_t1$chr<-as.numeric(as.character(unique_qtl_t1$chr))
unique_qtl_t1$pos<-as.numeric(as.character(unique_qtl_t1$pos))
chr<-unique_qtl_t1$chr
pos<-unique_qtl_t1$pos
qtl<-makeqtl(fg.cr.obj, unique_qtl_t1$chr, unique_qtl_t1$pos, what=c("prob"))
Qs<-paste('Q', 1:length(unique_qtl_t1$pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('chr', t1, 'all_qtl', sep="_"), chr)
assign(paste('pos', t1, 'all_qtl', sep="_"), pos)
assign(paste('qtl', t1, 'all_qtl', sep="_"), qtl)
assign(paste('my.formula', t1, 'all_qtl', sep="_"), my.formula)
#assign(paste('lodmat.F', t1, 'all_qtl', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(1,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_all_qtl.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="All QTL")
#refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
#dev.off()
par(mfrow=c(1,1))

#assign(paste('refqtlslod', t1, 'all_qtl', sep="_"), refqtlslod)

#klodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

#for(i in 1:length(klodeff)) {
#  klodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(klodeff[[1]])
#klodeff <- matrix(unlist(klodeff), byrow=TRUE, ncol=length(nam))
#colnames(klodeff) <- nam

#n.col<-ncol(klodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_all_qtl_', job, ".", t1, ".pdf", sep=""))
# Draw a plot of the intercept
#plot(days, klodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, klodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(klodeff)[i], side=3, line=0.5)
#}


#dev.off()
#par(mfrow=c(1,1))

#assign(paste('klodeff', t1, 'all_qtl', sep="_"), klodeff)

#write.csv(klodeff, file=paste('klodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)
# Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location

# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)

intercept.klod<-c()
summary.table.klod<-c()
for(i in get(paste(t1, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=qtl, formula=paste('my.formula', t1, 'all_qtl', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'klod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.klod$est$est[[1]]
  intercept.klod<-c(intercept.klod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.klod<-as.data.frame(out.fitqtl.klod$result.full[,2])
  full.mdl.var.klod$prop.variance<-c(100*(full.mdl.var.klod[1,1]/full.mdl.var.klod[3,1]), 100*(full.mdl.var.klod[2,1]/full.mdl.var.klod[3,1]), 100)
  colnames(full.mdl.var.klod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.klod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.klod<-rbind(dropone.mdl.var.klod, full.mdl.var.klod[3,1])
  dropone.mdl.var.klod$prop.variance<-c(100*(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-1),1]/full.mdl.var.klod[3,1]), 100)
  colnames(dropone.mdl.var.klod)<-c("variance", "prop.variance")
  
  
  rownames(dropone.mdl.var.klod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),1]-sum(dropone.mdl.var.klod[(1:nrow(dropone.mdl.var.klod)-1),1]), dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),2]-sum(dropone.mdl.var.klod[(1:nrow(dropone.mdl.var.klod)-1),2]))
  dropone.mdl.var.klod<-rbind(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-1),], dropone.residual, dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),])
  
  
  rownames(dropone.mdl.var.klod)[(nrow(dropone.mdl.var.klod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.klod<-chr
  CI_L<-c()
  CI_R<-c()
  for(q in n.QTL.klod){
    print(q)
    lod_int.klod<-lodint(out.F, chr=q)
    lod_int.klod$marker<-rownames(lod_int.klod)
    L<-as.data.frame(lod_int.klod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.klod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  }
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.klod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.klod<-as.data.frame(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.klod$result.drop[,4]
  klod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(klod.markers, marker.prop.var.klod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.klod<-rbind(summary.table.klod, summary.table)
}

assign(paste('summary.table.klod', t1, sep="."), summary.table.klod)
write.csv(summary.table.klod, file=paste('summary.table.klod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)

names(intercept.klod)<-days
assign(paste('intercept.klod', t1, sep="."), intercept.klod)
write.csv(intercept.klod, file=paste('intercept.klod.', t1, "_", job, '.csv', sep=""), quote=F, row.names=T)

save.image(file=paste('timeseries_', t1, '_cross.object.raw.Rdata', sep=""))

#################################################################################
# Treatment 2 (t2)
#################################################################################
setwd(dir.job.t2)

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][2])
}

out<- scanone(fg.cr.obj, pheno.col = get(paste(t2, 'cols', sep=".")), method="hk")
eff <- geteffects(fg.cr.obj, pheno.cols=get(paste(t2, 'cols', sep=".")))
pdf(paste('so_timeseries_qtl_', job, ".", t2, '.pdf', sep=""))
plotlod(out, eff, get(paste(t2, 'cols', sep=".")), gap=15, ylab="Time")
dev.off()

out.F <- scanoneF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), method="hk")
o.perm.F <- scanoneF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), method = "hk", n.perm=1000)
max.perm.F<-summary(o.perm.F, alpha=0.05)

qtl.table.so<-summary(out.F, perms=o.perm.F, alpha=0.05, pvalues=TRUE)
write.csv(qtl.table.so, file=paste('qtl.table.so.F_', job, ".", t2,'.csv', sep=""), quote=F, row.names=F)

par(mfrow=c(2,1))
pdf(paste('so_timeseries_qtl_slod_and_mlod_', job, '.', t2, '.pdf', sep=""))
plot(out.F, main="The SLOD curve", bandcol="gray90")
abline(h=max.perm.F[1], col="red", lty=3)
plot(out.F, lodcolumn=2, main="The MLOD curve", bandcol="gray90")
abline(h=max.perm.F[2], col="red", lty=3)
dev.off()
par(mfrow=c(1,1))

# MQM
# First SLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), max.qtl=9, usec= "slod", method="hk", penalties=c(max.perm.F[1],0,0))

if(out.mqm.F$n.qtl > 1) {
chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

#if (length(chr) == 1) {
#  temp<-summary(out.F)
#  temp<-temp[temp$chr != chr,]
#  temp<-temp[order(temp$slod, decreasing=T),]
#  chr<-c(chr, temp[1,'chr'])
#  pos<-c(pos, temp[1,'pos'])
#}

#if (length(chr) == 0) {
#  temp<-summary(out.F)
#  temp<-temp[order(temp$slod, decreasing=T),]
#  chr<-temp[1:2, 'chr']
#  pos<-temp[1:2, 'pos']
#}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl = qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t2, 'slod', sep="_"), out.mqm.F)
assign(paste('chr', t2, 'slod', sep="_"), chr)
assign(paste('pos', t2, 'slod', sep="_"), pos)
assign(paste('qtl', t2, 'slod', sep="_"), qtl)
assign(paste('my.formula', t2, 'slod', sep="_"), my.formula)
#assign(paste('lodmat.F', t2, 'slod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(2,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_slod.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="SLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
#dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod', t2, 'slod', sep="_"), refqtlslod)

#slodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

#for(i in 1:length(slodeff)) {
#  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(slodeff[[1]])
#slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
#colnames(slodeff) <- nam

#n.col<-ncol(slodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
#plot(days, slodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, slodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(slodeff)[i], side=3, line=0.5)
#}

dev.off()
par(mfrow=c(1,1))

#assign(paste('slodeff', t2, 'slod', sep="_"), slodeff)

#write.csv(slodeff, file=paste('slodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)

# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)

# Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
intercept.slod<-c()
summary.table.slod<-c()
for(i in get(paste(t2, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlslod, formula=paste('my.formula', t2, 'slod', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'slod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.slod$est$est[[1]]
  intercept.slod<-c(intercept.slod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.slod<-as.data.frame(out.fitqtl.slod$result.full[,2])
  full.mdl.var.slod$prop.variance<-c(100*(full.mdl.var.slod[1,1]/full.mdl.var.slod[3,1]), 100*(full.mdl.var.slod[2,1]/full.mdl.var.slod[3,1]), 100)
  colnames(full.mdl.var.slod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.slod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.slod<-rbind(dropone.mdl.var.slod, full.mdl.var.slod[3,1])
  dropone.mdl.var.slod$prop.variance<-c(100*(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-1),1]/full.mdl.var.slod[3,1]), 100)
  colnames(dropone.mdl.var.slod)<-c("variance", "prop.variance")
  
  rownames(dropone.mdl.var.slod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),1]-sum(dropone.mdl.var.slod[(1:nrow(dropone.mdl.var.slod)-1),1]), dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),2]-sum(dropone.mdl.var.slod[(1:nrow(dropone.mdl.var.slod)-1),2]))
  dropone.mdl.var.slod<-rbind(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-1),], dropone.residual, dropone.mdl.var.slod[nrow(dropone.mdl.var.slod),])
  
  
  rownames(dropone.mdl.var.slod)[(nrow(dropone.mdl.var.slod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.slod<-refqtlslod$n.qtl
  CI_L<-c()
  CI_R<-c()
  for(q in 1:n.QTL.slod){
    lod_int.slod<-lodint(refqtlslod, qtl.index=q)
    lod_int.slod$marker<-rownames(lod_int.slod)
    L<-as.data.frame(lod_int.slod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.slod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  }
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.slod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.slod<-as.data.frame(dropone.mdl.var.slod[1:(nrow(dropone.mdl.var.slod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.slod$result.drop[,4]
  slod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(slod.markers, marker.prop.var.slod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.slod<-rbind(summary.table.slod, summary.table)
}

assign(paste('summary.table.slod', t2, sep="."), summary.table.slod)
write.csv(summary.table.slod, file=paste('summary.table.slod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)

assign(paste('intercept.slod', t2, sep="."), intercept.slod)
write.csv(intercept.slod, file=paste('intercept.slod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)

}

if(out.mqm.F$n.qtl == 1){
  
  chr<-out.mqm.F$chr
  pos<-out.mqm.F$pos
  
  #if (length(chr) == 1) {
  #  temp<-summary(out.F)
  #  temp<-temp[temp$chr != chr,]
  #  temp<-temp[order(temp$slod, decreasing=T),]
  #  chr<-c(chr, temp[1,'chr'])
  #  pos<-c(pos, temp[1,'pos'])
  #}
  
  #if (length(chr) == 0) {
  #  temp<-summary(out.F)
  #  temp<-temp[order(temp$slod, decreasing=T),]
  #  chr<-temp[1:2, 'chr']
  #  pos<-temp[1:2, 'pos']
  #}
  
  qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
  Qs<-paste('Q', 1:length(pos), sep="")
  my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
  #lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")
  
  assign(paste('out.mqm.F', t2, 'slod', sep="_"), out.mqm.F)
  assign(paste('chr', t2, 'slod', sep="_"), chr)
  assign(paste('pos', t2, 'slod', sep="_"), pos)
  assign(paste('qtl', t2, 'slod', sep="_"), qtl)
  assign(paste('my.formula', t2, 'slod', sep="_"), my.formula)
  #assign(paste('lodmat.F', t2, 'slod', sep="_"), lodmat.F)
  
  # Make plots of the QTL LOD profile and the LOD profile over time
  #par(mfrow=c(2,1))
  #pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_slod.pdf', sep=""))
  #plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="slod")
  refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
  #plotLodProfile(refqtlslod)
  #dev.off()
  par(mfrow=c(1,1))
  
  assign(paste('refqtlslod', t2, 'slod', sep="_"), refqtlslod)
  #slodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))
  
  #for(i in 1:length(slodeff)) {
  #  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
  #                                 method="hk", get.ests=TRUE,
  #                                 dropone=FALSE))$ests[,1]*c(1,2,2)
  #}
  
  pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')
  
  days<-c()
  for(p in 1:length(pname_list)) {
    days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
  }
  
  #nam <- names(slodeff[[1]])
  #slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
  #colnames(slodeff) <- nam
  
  #n.col<-ncol(slodeff)
  
  # Lets plot the effect size over time
  #par(mfrow=c(1,n.col))
  #pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t2, ".pdf", sep=""))
  # Draw a plot of the intercept
  #plot(days, slodeff[,1], lwd=2, type="l",
  #     xlab="Days after planting",
  #     ylab="Height (cm)", col="red")
  #mtext("baseline curve", side=3, line=0.5)
  # Now add plot of effect size for each QTL
  #for (i in 2:n.col) {
  #  plot(days, slodeff[,i], lwd=2, type="l",
  #       xlab="Days after planting",
  #       ylab="QTL effect (cm)", col="red")
  #  mtext(colnames(slodeff)[i], side=3, line=0.5)
  #}
  
  #dev.off()
  #par(mfrow=c(1,1))
  
  #assign(paste('slodeff', t2, 'slod', sep="_"), slodeff)
  #write.csv(slodeff, file=paste('slodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)
  # get marker names
  m.names<-find.marker(fg.cr.obj, chr, pos)
  
  # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
  summary.table.slod<-c()
  intercept.slod<-c()
  for(i in get(paste(t2, 'cols', sep="."))) {
    out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlslod, formula=paste('my.formula', t2, 'slod', sep="_"), method='hk', get.ests=T)
    assign(paste('out.fitqtl', 'slod', sep="."), out.fitqtl)
    
    # Get value of height for the intercept
    int<-out.fitqtl.slod$est$est[[1]]
    intercept.slod<-c(intercept.slod, int)
    
    # Proportioning of variance for scanone results (full model)
    full.mdl.var.slod<-as.data.frame(out.fitqtl.slod$result.full[,2])
    full.mdl.var.slod$prop.variance<-c(100*(full.mdl.var.slod[1,1]/full.mdl.var.slod[3,1]), 100*(full.mdl.var.slod[2,1]/full.mdl.var.slod[3,1]), 100)
    colnames(full.mdl.var.slod)<-c("variance", "prop.variance")
    
    
    # Get confidence interval (95%) for each QTL:
    n.QTL.slod<-refqtlslod$n.qtl
    CI_L<-c()
    CI_R<-c()
    for(q in 1:n.QTL.slod){
      lod_int.slod<-lodint(refqtlslod, qtl.index=q)
      lod_int.slod$marker<-rownames(lod_int.slod)
      L<-as.data.frame(lod_int.slod[1,c(4,1,2,3)])
      colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
      R<-as.data.frame(lod_int.slod[3,c(4,1,2,3)])
      colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
      CI_L<-rbind(CI_L, L)
      CI_R<-rbind(CI_R, R)
    }
    CI<-cbind(CI_L, CI_R)
    
    # Get proportion of variance for only marker
    marker.prop.var.slod<-full.mdl.var.slod$prop.variance[1]
    # Get the fx size for only marker
    fx_size.a<-out.fitqtl.slod$est$est[[2]][1]
    fx_size.se<-c(0)
    # Build a table that summarizes the results
    lods<-out.fitqtl.slod$lod
    slod.markers<-cbind(chr, pos, lods)
    phe.name<-colnames(fg.cr.obj$pheno)[i]
    summary.table<-c(slod.markers, marker.prop.var.slod, as.numeric(fx_size.a), fx_size.se, CI, phe.name)
    summary.table<-data.frame(matrix(unlist(summary.table), nrow=1, byrow=T),stringsAsFactors=FALSE)
    rownames(summary.table)<-m.names
    colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
    summary.table.slod<-rbind(summary.table.slod, summary.table)
  }
  
  assign(paste('summary.table.slod', t2, sep="."), summary.table.slod)
  write.csv(summary.table.slod, file=paste('summary.table.slod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)
  names(intercept.slod)<-days
  assign(paste('intercept.slod', t2, sep="."), intercept.slod)
  write.csv(intercept.slod, file=paste('intercept.slod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)
  
}

if(out.mqm.F$n.qtl < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t2, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))


# MQM
# Now MLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), max.qtl=9, usec= "mlod", method="hk", penalties=c(max.perm.F[2],0,0))

if(out.mqm.F$n.qtl > 1) {

chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

#if (length(chr) == 1) {
#  temp<-summary(out.F)
#  temp<-temp[temp$chr != chr,]
#  temp<-temp[order(temp$mlod, decreasing=T),]
#  chr<-c(chr, temp[1,'chr'])
#  pos<-c(pos, temp[1,'pos'])
#}

#if (length(chr) == 0) {
#  temp<-summary(out.F)
#  temp<-temp[order(temp$mlod, decreasing=T),]
#  chr<-temp[1:2, 'chr']
#  pos<-temp[1:2, 'pos']
#}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t2, 'mlod', sep="_"), out.mqm.F)
assign(paste('chr', t2, 'mlod', sep="_"), chr)
assign(paste('pos', t2, 'mlod', sep="_"), pos)
assign(paste('qtl', t2, 'mlod', sep="_"), qtl)
assign(paste('my.formula', t2, 'mlod', sep="_"), my.formula)
#assign(paste('lodmat.F', t2, 'mlod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(2,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_mlod.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
refqtlmlod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlmlod)
#dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlmlod', t2, 'mlod', sep="_"), refqtlmlod)
#mlodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

#for(i in 1:length(mlodeff)) {
#  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(mlodeff[[1]])
#mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
#colnames(mlodeff) <- nam

#n.col<-ncol(mlodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_mlod_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
#plot(days, mlodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, mlodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(mlodeff)[i], side=3, line=0.5)
#}

#dev.off()
#par(mfrow=c(1,1))

#assign(paste('mlodeff', t2, 'mlod', sep="_"), mlodeff)
#write.csv(mlodeff, file=paste('mlodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)
# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)

# Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
intercept.mlod<-c()
summary.table.mlod<-c()
for(i in get(paste(t2, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlmlod, formula=paste('my.formula', t2, 'mlod', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'mlod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.mlod$est$est[[1]]
  intercept.mlod<-c(intercept.mlod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.mlod<-as.data.frame(out.fitqtl.mlod$result.full[,2])
  full.mdl.var.mlod$prop.variance<-c(100*(full.mdl.var.mlod[1,1]/full.mdl.var.mlod[3,1]), 100*(full.mdl.var.mlod[2,1]/full.mdl.var.mlod[3,1]), 100)
  colnames(full.mdl.var.mlod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.mlod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.mlod<-rbind(dropone.mdl.var.mlod, full.mdl.var.mlod[3,1])
  dropone.mdl.var.mlod$prop.variance<-c(100*(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-1),1]/full.mdl.var.mlod[3,1]), 100)
  colnames(dropone.mdl.var.mlod)<-c("variance", "prop.variance")
  
  rownames(dropone.mdl.var.mlod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),1]-sum(dropone.mdl.var.mlod[(1:nrow(dropone.mdl.var.mlod)-1),1]), dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),2]-sum(dropone.mdl.var.mlod[(1:nrow(dropone.mdl.var.mlod)-1),2]))
  dropone.mdl.var.mlod<-rbind(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-1),], dropone.residual, dropone.mdl.var.mlod[nrow(dropone.mdl.var.mlod),])
  
  
  rownames(dropone.mdl.var.mlod)[(nrow(dropone.mdl.var.mlod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.mlod<-refqtlmlod$n.qtl
  CI_L<-c()
  CI_R<-c()
  for(q in 1:n.QTL.mlod){
    lod_int.mlod<-lodint(refqtlmlod, qtl.index=q)
    lod_int.mlod$marker<-rownames(lod_int.mlod)
    L<-as.data.frame(lod_int.mlod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.mlod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  }
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mlod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.mlod<-as.data.frame(dropone.mdl.var.mlod[1:(nrow(dropone.mdl.var.mlod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.mlod$result.drop[,4]
  mlod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(mlod.markers, marker.prop.var.mlod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.mlod<-rbind(summary.table.mlod, summary.table)
}

assign(paste('summary.table.mlod', t2, sep="."), summary.table.mlod)
write.csv(summary.table.mlod, file=paste('summary.table.mlod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)

names(intercept.mlod)<-days
assign(paste('intercept.mlod', t2, sep="."), intercept.mlod)
write.csv(intercept.mlod, file=paste('intercept.mlod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)


}

if(out.mqm.F$n.qtl == 1){
  
  chr<-out.mqm.F$chr
  pos<-out.mqm.F$pos
  
  #if (length(chr) == 1) {
  #  temp<-summary(out.F)
  #  temp<-temp[temp$chr != chr,]
  #  temp<-temp[order(temp$mlod, decreasing=T),]
  #  chr<-c(chr, temp[1,'chr'])
  #  pos<-c(pos, temp[1,'pos'])
  #}
  
  #if (length(chr) == 0) {
  #  temp<-summary(out.F)
  #  temp<-temp[order(temp$mlod, decreasing=T),]
  #  chr<-temp[1:2, 'chr']
  #  pos<-temp[1:2, 'pos']
  #}
  
  qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
  Qs<-paste('Q', 1:length(pos), sep="")
  my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
  #lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")
  
  assign(paste('out.mqm.F', t2, 'mlod', sep="_"), out.mqm.F)
  assign(paste('chr', t2, 'mlod', sep="_"), chr)
  assign(paste('pos', t2, 'mlod', sep="_"), pos)
  assign(paste('qtl', t2, 'mlod', sep="_"), qtl)
  assign(paste('my.formula', t2, 'mlod', sep="_"), my.formula)
  #assign(paste('lodmat.F', t2, 'mlod', sep="_"), lodmat.F)
  
  # Make plots of the QTL LOD profile and the LOD profile over time
  #par(mfrow=c(2,1))
  #pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_mlod.pdf', sep=""))
  #plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
  refqtlmlod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
  #plotLodProfile(refqtlmlod)
  #dev.off()
  par(mfrow=c(1,1))
  
  assign(paste('refqtlmlod', t2, 'mlod', sep="_"), refqtlmlod)
  #mlodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))
  
  #for(i in 1:length(mlodeff)) {
  #  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
  #                                 method="hk", get.ests=TRUE,
  #                                 dropone=FALSE))$ests[,1]*c(1,2,2)
  #}
  
  pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')
  
  days<-c()
  for(p in 1:length(pname_list)) {
    days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
  }
  
  #nam <- names(mlodeff[[1]])
  #mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
  #colnames(mlodeff) <- nam
  
  #n.col<-ncol(mlodeff)
  
  # Lets plot the effect size over time
  #par(mfrow=c(1,n.col))
  #pdf(paste('mqm_timeseries_fx_size_mlod_', job, ".", t2, ".pdf", sep=""))
  # Draw a plot of the intercept
  #plot(days, mlodeff[,1], lwd=2, type="l",
  #     xlab="Days after planting",
  #     ylab="Height (cm)", col="red")
  #mtext("baseline curve", side=3, line=0.5)
  # Now add plot of effect size for each QTL
  #for (i in 2:n.col) {
  #  plot(days, mlodeff[,i], lwd=2, type="l",
  #       xlab="Days after planting",
  #       ylab="QTL effect (cm)", col="red")
  #  mtext(colnames(mlodeff)[i], side=3, line=0.5)
  #}
  
  #dev.off()
  #par(mfrow=c(1,1))
  
  #assign(paste('mlodeff', t2, 'mlod', sep="_"), mlodeff)
  #write.csv(mlodeff, file=paste('mlodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)
  # get marker names
  m.names<-find.marker(fg.cr.obj, chr, pos)
  
  # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
  summary.table.mlod<-c()
  intercept.mlod<-c()
  for(i in get(paste(t2, 'cols', sep="."))) {
    out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=refqtlmlod, formula=paste('my.formula', t2, 'mlod', sep="_"), method='hk', get.ests=T)
    assign(paste('out.fitqtl', 'mlod', sep="."), out.fitqtl)
    
    # Get value of height for the intercept
    int<-out.fitqtl.mlod$est$est[[1]]
    intercept.mlod<-c(intercept.mlod, int)
    
    # Proportioning of variance for scanone results (full model)
    full.mdl.var.mlod<-as.data.frame(out.fitqtl.mlod$result.full[,2])
    full.mdl.var.mlod$prop.variance<-c(100*(full.mdl.var.mlod[1,1]/full.mdl.var.mlod[3,1]), 100*(full.mdl.var.mlod[2,1]/full.mdl.var.mlod[3,1]), 100)
    colnames(full.mdl.var.mlod)<-c("variance", "prop.variance")
    
    
    # Get confidence interval (95%) for each QTL:
    n.QTL.mlod<-refqtlmlod$n.qtl
    CI_L<-c()
    CI_R<-c()
    for(q in 1:n.QTL.mlod){
      lod_int.mlod<-lodint(refqtlmlod, qtl.index=q)
      lod_int.mlod$marker<-rownames(lod_int.mlod)
      L<-as.data.frame(lod_int.mlod[1,c(4,1,2,3)])
      colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
      R<-as.data.frame(lod_int.mlod[3,c(4,1,2,3)])
      colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
      CI_L<-rbind(CI_L, L)
      CI_R<-rbind(CI_R, R)
    }
    CI<-cbind(CI_L, CI_R)
    
    # Get proportion of variance for only marker
    marker.prop.var.mlod<-full.mdl.var.mlod$prop.variance[1]
    # Get the fx size for only marker
    fx_size.a<-out.fitqtl.mlod$est$est[[2]][1]
    fx_size.se<-c(0)
    # Build a table that summarizes the results
    lods<-out.fitqtl.mlod$lod
    mlod.markers<-cbind(chr, pos, lods)
    phe.name<-colnames(fg.cr.obj$pheno)[i]
    summary.table<-c(mlod.markers, marker.prop.var.mlod, as.numeric(fx_size.a), fx_size.se, CI, phe.name)
    summary.table<-data.frame(matrix(unlist(summary.table), nrow=1, byrow=T),stringsAsFactors=FALSE)
    rownames(summary.table)<-m.names
    colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
    summary.table.mlod<-rbind(summary.table.mlod, summary.table)
  }
  
  assign(paste('summary.table.mlod', t2, sep="."), summary.table.mlod)
  write.csv(summary.table.mlod, file=paste('summary.table.mlod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)
  names(intercept.mlod)<-days
  assign(paste('intercept.mlod', t2, sep="."), intercept.mlod)
  write.csv(intercept.mlod, file=paste('intercept.mlod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)
  
}


if(out.mqm.F$n.qtl < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t2, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))


######## Now lets do the same thing all unique qtl at all time points for t2
unique_qtl_t2$chr<-as.numeric(as.character(unique_qtl_t2$chr))
unique_qtl_t2$pos<-as.numeric(as.character(unique_qtl_t2$pos))
chr<-as.numeric(as.character(unique_qtl_t2$chr))
pos<-as.numeric(as.character(unique_qtl_t2$pos))

qtl<-makeqtl(fg.cr.obj, unique_qtl_t2$chr, unique_qtl_t2$pos, what=c("prob"))
Qs<-paste('Q', 1:length(unique_qtl_t2$pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
#lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('chr', t2, 'all_qtl', sep="_"), chr)
assign(paste('pos', t2, 'all_qtl', sep="_"), pos)
assign(paste('qtl', t2, 'all_qtl', sep="_"), qtl)
assign(paste('my.formula', t2, 'all_qtl', sep="_"), my.formula)
#assign(paste('lodmat.F', t2, 'all_qtl', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
#par(mfrow=c(1,1))
#pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_all_qtl.pdf', sep=""))
#plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="All QTL")
#refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
#dev.off()
par(mfrow=c(1,1))

#assign(paste('refqtlslod', t2, 'all_qtl', sep="_"), refqtlslod)

#klodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

#for(i in 1:length(klodeff)) {
#  klodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
#                                 method="hk", get.ests=TRUE,
#                                 dropone=FALSE))$ests[,1]*c(1,2,2)
#}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

#nam <- names(klodeff[[1]])
#klodeff <- matrix(unlist(klodeff), byrow=TRUE, ncol=length(nam))
#colnames(klodeff) <- nam

#n.col<-ncol(klodeff)

# Lets plot the effect size over time
#par(mfrow=c(1,n.col))
#pdf(paste('mqm_timeseries_fx_size_all_qtl_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
#plot(days, klodeff[,1], lwd=2, type="l",
#     xlab="Days after planting",
#     ylab="Height (cm)", col="red")
#mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
#for (i in 2:n.col) {
#  plot(days, klodeff[,i], lwd=2, type="l",
#       xlab="Days after planting",
#       ylab="QTL effect (cm)", col="red")
#  mtext(colnames(klodeff)[i], side=3, line=0.5)
#}


#dev.off()
#par(mfrow=c(1,1))

#assign(paste('klodeff', t2, 'all_qtl', sep="_"), klodeff)

#write.csv(klodeff, file=paste('klodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)


# get marker names
m.names<-find.marker(fg.cr.obj, chr, pos)
intercept.klod<-c()
summary.table.klod<-c()
for(i in get(paste(t2, 'cols', sep="."))) {
  out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=qtl, formula=paste('my.formula', t2, 'all_qtl', sep="_"), method='hk', get.ests=T)
  assign(paste('out.fitqtl', 'klod', sep="."), out.fitqtl)
  
  # Get value of height for the intercept
  int<-out.fitqtl.klod$est$est[[1]]
  intercept.klod<-c(intercept.klod, int)
  
  # Proportioning of variance for scanone results (full model)
  full.mdl.var.klod<-as.data.frame(out.fitqtl.klod$result.full[,2])
  full.mdl.var.klod$prop.variance<-c(100*(full.mdl.var.klod[1,1]/full.mdl.var.klod[3,1]), 100*(full.mdl.var.klod[2,1]/full.mdl.var.klod[3,1]), 100)
  colnames(full.mdl.var.klod)<-c("variance", "prop.variance")
  
  #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
  
  # Proportioning of variance for scanone results (individual QTL as proportion)
  dropone.mdl.var.klod<-as.data.frame(out.fitqtl$result.drop[,2])
  dropone.mdl.var.klod<-rbind(dropone.mdl.var.klod, full.mdl.var.klod[3,1])
  dropone.mdl.var.klod$prop.variance<-c(100*(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-1),1]/full.mdl.var.klod[3,1]), 100)
  colnames(dropone.mdl.var.klod)<-c("variance", "prop.variance")
  
  
  rownames(dropone.mdl.var.klod)<-c(m.names, "total")
  
  dropone.residual<-c(dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),1]-sum(dropone.mdl.var.klod[(1:nrow(dropone.mdl.var.klod)-1),1]), dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),2]-sum(dropone.mdl.var.klod[(1:nrow(dropone.mdl.var.klod)-1),2]))
  dropone.mdl.var.klod<-rbind(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-1),], dropone.residual, dropone.mdl.var.klod[nrow(dropone.mdl.var.klod),])
  
  
  rownames(dropone.mdl.var.klod)[(nrow(dropone.mdl.var.klod)-1)]<-c('residual')
  
  # Get confidence interval (95%) for each QTL:
  n.QTL.klod<-length(chr)
  CI_L<-c()
  CI_R<-c()
  for(q in chr){
    print(q)
    lod_int.klod<-lodint(out.F, chr=q)
    lod_int.klod$marker<-rownames(lod_int.klod)
    L<-as.data.frame(lod_int.klod[1,c(4,1,2,3)])
    colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
    R<-as.data.frame(lod_int.klod[3,c(4,1,2,3)])
    colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
    CI_L<-rbind(CI_L, L)
    CI_R<-rbind(CI_R, R)
  } 
  CI<-cbind(CI_L, CI_R)
  
  # Lets get addative effect sizes from the fitqtl model of scanone qtl
  fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.klod)[3])
  fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
  fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
  # Get proportion of variance per marker
  marker.prop.var.klod<-as.data.frame(dropone.mdl.var.klod[1:(nrow(dropone.mdl.var.klod)-2),2])
  
  # Build a table that summarizes the results
  lods<-out.fitqtl.klod$result.drop[,4]
  klod.markers<-cbind(chr, pos, lods)
  phe.name<-colnames(fg.cr.obj$pheno)[i]
  summary.table<-cbind(klod.markers, marker.prop.var.klod, fx_size.a, fx_size.se, CI, phe.name)
  rownames(summary.table)<-m.names
  colnames(summary.table)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI), 'trait')
  summary.table.klod<-rbind(summary.table.klod, summary.table)
}


assign(paste('summary.table.klod', t2, sep="."), summary.table.klod)
write.csv(summary.table.klod, file=paste('summary.table.klod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)

names(intercept.klod)<-days
assign(paste('intercept.klod', t2, sep="."), intercept.klod)

write.csv(intercept.klod, file=paste('intercept.klod.', t2, "_", job, '.csv', sep=""), quote=F, row.names=T)


save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))
