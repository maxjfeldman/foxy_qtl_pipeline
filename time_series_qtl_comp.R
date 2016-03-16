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
if(length(out.mqm.F) > 0 ) {

chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

if (length(chr) == 1) {
  temp<-summary(out.F)
  temp<-temp[temp$chr != chr,]
  temp<-temp[order(temp$slod, decreasing=T),]
  chr<-c(chr, temp[1,'chr'])
  pos<-c(pos, temp[1,'pos'])
}

if (length(chr) == 0) {
  temp<-summary(out.F)
  temp<-temp[order(temp$slod, decreasing=T),]
  chr<-temp[1:2, 'chr']
  pos<-temp[1:2, 'pos']
}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t1, 'slod', sep="_"), out.mqm.F)
assign(paste('chr', t1, 'slod', sep="_"), chr)
assign(paste('pos', t1, 'slod', sep="_"), pos)
assign(paste('qtl', t1, 'slod', sep="_"), qtl)
assign(paste('my.formula', t1, 'slod', sep="_"), my.formula)
assign(paste('lodmat.F', t1, 'slod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_slod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="SLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod ', t1, 'slod', sep="_"), refqtlslod)

slodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

for(i in 1:length(slodeff)) {
  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(slodeff[[1]])
slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
colnames(slodeff) <- nam

n.col<-ncol(slodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t1, ".pdf", sep=""))
# Draw a plot of the intercept
plot(days, slodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, slodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(slodeff)[i], side=3, line=0.5)
}

dev.off()
par(mfrow=c(1,1))

assign(paste('slodeff ', t1, 'slod', sep="_"), slodeff)

write.csv(slodeff, file=paste('slodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)
}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t1, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t1, '_cross.object.raw.Rdata', sep=""))

# MQM
# Now MLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), max.qtl=9, usec= "mlod", method="hk", penalties=c(max.perm.F[2],0,0))

if(length(out.mqm.F) > 0) {
 
chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

if (length(chr) == 1) {
  temp<-summary(out.F)
  temp<-temp[temp$chr != chr,]
  temp<-temp[order(temp$mlod, decreasing=T),]
  chr<-c(chr, temp[1,'chr'])
  pos<-c(pos, temp[1,'pos'])
}

if (length(chr) == 0) {
  temp<-summary(out.F)
  temp<-temp[order(temp$mlod, decreasing=T),]
  chr<-temp[1:2, 'chr']
  pos<-temp[1:2, 'pos']
}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t1, 'mlod', sep="_"), out.mqm.F)
assign(paste('chr', t1, 'mlod', sep="_"), chr)
assign(paste('pos', t1, 'mlod', sep="_"), pos)
assign(paste('qtl', t1, 'mlod', sep="_"), qtl)
assign(paste('my.formula', t1, 'mlod', sep="_"), my.formula)
assign(paste('lodmat.F', t1, 'mlod', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_mlod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlmlod ', t1, 'mlod', sep="_"), refqtlslod)


mlodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

for(i in 1:length(mlodeff)) {
  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')
days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(mlodeff[[1]])
mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
colnames(mlodeff) <- nam

n.col<-ncol(mlodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_mlod_',job, ".", t1, '.pdf', sep=""))
# Draw a plot of the intercept
plot(days, mlodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, mlodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(mlodeff)[i], side=3, line=0.5)
}


dev.off()
par(mfrow=c(1,1))
assign(paste('mlodeff ', t1, 'mlod', sep="_"), mlodeff)
write.csv(mlodeff, file=paste('mlodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)

}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t1, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}
save.image(file=paste('timeseries_', t1, '_cross.object.raw.Rdata', sep=""))


######## Now lets do the same thing all unique qtl at all time points for t1
unique_qtl_t1$chr<-as.numeric(as.character(unique_qtl_t1$chr))
unique_qtl_t1$pos<-as.numeric(as.character(unique_qtl_t1$pos))
qtl<-makeqtl(fg.cr.obj, unique_qtl_t1$chr, unique_qtl_t1$pos, what=c("prob"))
Qs<-paste('Q', 1:length(unique_qtl_t1$pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t1, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('chr', t1, 'all_qtl', sep="_"), chr)
assign(paste('pos', t1, 'all_qtl', sep="_"), pos)
assign(paste('qtl', t1, 'all_qtl', sep="_"), qtl)
assign(paste('my.formula', t1, 'all_qtl', sep="_"), my.formula)
assign(paste('lodmat.F', t1, 'all_qtl', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(1,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t1, '_all_qtl.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="All QTL")
#refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t1, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

#assign(paste('refqtlslod ', t1, 'all_qtl', sep="_"), refqtlslod)

klodeff <- vector("list", length(get(paste(t1, 'cols', sep="."))))

for(i in 1:length(klodeff)) {
  klodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t1, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t1, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(klodeff[[1]])
klodeff <- matrix(unlist(klodeff), byrow=TRUE, ncol=length(nam))
colnames(klodeff) <- nam

n.col<-ncol(klodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_all_qtl_', job, ".", t1, ".pdf", sep=""))
# Draw a plot of the intercept
plot(days, klodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, klodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(klodeff)[i], side=3, line=0.5)
}


dev.off()
par(mfrow=c(1,1))

assign(paste('klodeff', t1, 'all_qtl', sep="_"), klodeff)

write.csv(klodeff, file=paste('klodeff_', t1, "_", job, '.csv', sep=""), quote=F, row.names=F)
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

if(length(out.mqm.F) > 0) {
chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

if (length(chr) == 1) {
  temp<-summary(out.F)
  temp<-temp[temp$chr != chr,]
  temp<-temp[order(temp$slod, decreasing=T),]
  chr<-c(chr, temp[1,'chr'])
  pos<-c(pos, temp[1,'pos'])
}

if (length(chr) == 0) {
  temp<-summary(out.F)
  temp<-temp[order(temp$slod, decreasing=T),]
  chr<-temp[1:2, 'chr']
  pos<-temp[1:2, 'pos']
}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t2, 'slod', sep="_"), out.mqm.F)
assign(paste('chr', t2, 'slod', sep="_"), chr)
assign(paste('pos', t2, 'slod', sep="_"), pos)
assign(paste('qtl', t2, 'slod', sep="_"), qtl)
assign(paste('my.formula', t2, 'slod', sep="_"), my.formula)
assign(paste('lodmat.F', t2, 'slod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_slod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="SLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod ', t2, 'slod', sep="_"), refqtlslod)

slodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

for(i in 1:length(slodeff)) {
  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(slodeff[[1]])
slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
colnames(slodeff) <- nam

n.col<-ncol(slodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_slod_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
plot(days, slodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, slodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(slodeff)[i], side=3, line=0.5)
}

dev.off()
par(mfrow=c(1,1))

assign(paste('slodeff ', t2, 'slod', sep="_"), slodeff)

write.csv(slodeff, file=paste('slodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)

}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t2, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))


# MQM
# Now MLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), max.qtl=9, usec= "mlod", method="hk", penalties=c(max.perm.F[2],0,0))

if(length(out.mqm.F) > 0) {

chr<-out.mqm.F$chr
pos<-out.mqm.F$pos

if (length(chr) == 1) {
  temp<-summary(out.F)
  temp<-temp[temp$chr != chr,]
  temp<-temp[order(temp$mlod, decreasing=T),]
  chr<-c(chr, temp[1,'chr'])
  pos<-c(pos, temp[1,'pos'])
}

if (length(chr) == 0) {
  temp<-summary(out.F)
  temp<-temp[order(temp$mlod, decreasing=T),]
  chr<-temp[1:2, 'chr']
  pos<-temp[1:2, 'pos']
}

qtl<-makeqtl(fg.cr.obj, chr, pos, what=c("prob"))
Qs<-paste('Q', 1:length(pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', t2, 'mlod', sep="_"), out.mqm.F)
assign(paste('chr', t2, 'mlod', sep="_"), chr)
assign(paste('pos', t2, 'mlod', sep="_"), pos)
assign(paste('qtl', t2, 'mlod', sep="_"), qtl)
assign(paste('my.formula', t2, 'mlod', sep="_"), my.formula)
assign(paste('lodmat.F', t2, 'mlod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_mlod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod ', t2, 'mlod', sep="_"), refqtlslod)
mlodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

for(i in 1:length(mlodeff)) {
  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(mlodeff[[1]])
mlodeff <- matrix(unlist(mlodeff), byrow=TRUE, ncol=length(nam))
colnames(mlodeff) <- nam

n.col<-ncol(mlodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_mlod_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
plot(days, mlodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, mlodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(mlodeff)[i], side=3, line=0.5)
}

dev.off()
par(mfrow=c(1,1))

assign(paste('mlodeff ', t2, 'mlod', sep="_"), mlodeff)
write.csv(mlodeff, file=paste('mlodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)
}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job.t2, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))


######## Now lets do the same thing all unique qtl at all time points for t2
unique_qtl_t2$chr<-as.numeric(as.character(unique_qtl_t2$chr))
unique_qtl_t2$pos<-as.numeric(as.character(unique_qtl_t2$pos))
qtl<-makeqtl(fg.cr.obj, unique_qtl_t2$chr, unique_qtl_t2$pos, what=c("prob"))
Qs<-paste('Q', 1:length(unique_qtl_t2$pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = get(paste(t2, 'cols', sep=".")), formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('chr', t2, 'all_qtl', sep="_"), chr)
assign(paste('pos', t2, 'all_qtl', sep="_"), pos)
assign(paste('qtl', t2, 'all_qtl', sep="_"), qtl)
assign(paste('my.formula', t2, 'all_qtl', sep="_"), my.formula)
assign(paste('lodmat.F', t2, 'all_qtl', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(1,1))
pdf(paste('mqm_timeseries_qtl_', job, ".", t2, '_all_qtl.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="All QTL")
#refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = get(paste(t2, 'cols', sep=".")), usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

#assign(paste('refqtlslod ', t2, 'all_qtl', sep="_"), refqtlslod)

klodeff <- vector("list", length(get(paste(t2, 'cols', sep="."))))

for(i in 1:length(klodeff)) {
  klodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(get(paste(t2, 'cols', sep=".")))-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[get(paste(t2, 'cols', sep="."))], '_')

days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

nam <- names(klodeff[[1]])
klodeff <- matrix(unlist(klodeff), byrow=TRUE, ncol=length(nam))
colnames(klodeff) <- nam

n.col<-ncol(klodeff)

# Lets plot the effect size over time
par(mfrow=c(1,n.col))
pdf(paste('mqm_timeseries_fx_size_all_qtl_', job, ".", t2, ".pdf", sep=""))
# Draw a plot of the intercept
plot(days, klodeff[,1], lwd=2, type="l",
     xlab="Days after planting",
     ylab="Height (cm)", col="red")
mtext("baseline curve", side=3, line=0.5)
# Now add plot of effect size for each QTL
for (i in 2:n.col) {
  plot(days, klodeff[,i], lwd=2, type="l",
       xlab="Days after planting",
       ylab="QTL effect (cm)", col="red")
  mtext(colnames(klodeff)[i], side=3, line=0.5)
}


dev.off()
par(mfrow=c(1,1))

assign(paste('klodeff', t2, 'all_qtl', sep="_"), klodeff)

write.csv(klodeff, file=paste('klodeff_', t2, "_", job, '.csv', sep=""), quote=F, row.names=F)
save.image(file=paste('timeseries_', t2, '_cross.object.raw.Rdata', sep=""))
