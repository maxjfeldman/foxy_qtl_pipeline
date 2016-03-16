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

# Get directory path and names of tables
dir.job<-paste(directory, job, sep="/")
dir.job.diff<-paste(dir.job, '/', 'timeseries.diff', sep="")
dir.create(dir.job.diff)

# Make table name
st.name<-paste(job, 'concatenated_summary_table.csv', sep="_")
st.path<-paste(dir.job, st.name, sep="/")
st<-read.csv(st.path, header=F)
colnames(st)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')
st<-st[st$type == 'comp',]


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


###### Lets use that function to get a list of unique qtl found in diff
chrs<-sort(unique(st$chr))
m.names<-c()
chr<-c()
pos<-c()
for(ch in 1:length(chrs)) {
  temp<-st[st$chr == chrs[ch],]
  temp$marker<-as.character(temp$marker)
  remove_dup_qtl(temp)
}

# Combine the marker names, chromosome and positional info into a data.frame
unique_qtl<-as.data.frame(cbind(m.names, chr, pos))

# Lets load in the cross object for diff phenotypes 
diff_cr.obj_path<-paste(dir.job, 'cross.object.diff.Rdata', sep="/")
load(diff_cr.obj_path)

cols<-(2:length(phenames(fg.cr.obj)))

fg.cr.obj <- calc.genoprob(fg.cr.obj, step=0)

pname_list<-strsplit(phenames(fg.cr.obj)[cols], '_')
days<-c()
for(p in 1:length(pname_list)) {
  days<-c(days,pname_list[[p]][length(pname_list[[1]])-1])
}

out<- scanone(fg.cr.obj, pheno.col = cols, method="hk")
eff <- geteffects(fg.cr.obj, pheno.cols=cols)
setwd(dir.job.diff)
pdf(paste('so_timeseries_qtl_',job, ".", 'diff.pdf', sep=""))
plotlod(out, eff, cols, gap=15, ylab="Time")
dev.off()

out.F <- scanoneF(fg.cr.obj, pheno.cols = cols, method="hk")
o.perm.F <- scanoneF(fg.cr.obj, pheno.cols = cols, method = "hk", n.perm=1000)
max.perm.F<-summary(o.perm.F, alpha=0.05)

qtl.table.so<-summary(out.F, perms=o.perm.F, alpha=0.05, pvalues=TRUE)
write.csv(qtl.table.so, file=paste('qtl.table.so.F_', job, ".",'diff.csv', sep=""), quote=F, row.names=F)

par(mfrow=c(2,1))
pdf(paste('so_timeseries_qtl_slod_and_mlod_', job, ".", 'diff.pdf' ,sep=""))
plot(out.F, main="The SLOD curve", bandcol="gray90")
abline(h=max.perm.F[1], col="red", lty=3)
plot(out.F, lodcolumn=2, main="The MLOD curve", bandcol="gray90")
abline(h=max.perm.F[2], col="red", lty=3)
dev.off()
par(mfrow=c(1,1))


# MQM
# First SLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = cols, max.qtl=9, usec= "slod", method="hk", penalties=c(max.perm.F[1],0,0))

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
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = cols, formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', 'diff', 'slod', sep="_"), out.mqm.F)
assign(paste('chr', 'diff', 'slod', sep="_"), chr)
assign(paste('pos', 'diff', 'slod', sep="_"), pos)
assign(paste('qtl', 'diff', 'slod', sep="_"), qtl)
assign(paste('my.formula', 'diff', 'slod', sep="_"), my.formula)
assign(paste('lodmat.F', 'diff', 'slod', sep="_"), lodmat.F)

# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, '.diff', '_slod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="SLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = cols, usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlslod ', 'diff', 'slod', sep="_"), refqtlslod)

slodeff <- vector("list", length(cols))

for(i in 1:length(slodeff)) {
  slodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(cols)-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[cols], '_')

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
pdf(paste('mqm_timeseries_fx_size_slod_', job, ".diff.pdf", sep=""))
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

assign(paste('slodeff ', 'diff', 'slod', sep="_"), slodeff)
write.csv(slodeff, file=paste('slodeff_', 'diff', "_", job, '.csv', sep=""), quote=F, row.names=F)

}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job, '/mqm_timeseries_qtl_', job, '.diff', '_slod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', 'diff', '_cross.object.diff.Rdata', sep=""))

  
# MQM
# Now MLOD
out.mqm.F<-stepwiseqtlF(fg.cr.obj, pheno.cols = cols, max.qtl=9, usec= "mlod", method="hk", penalties=c(max.perm.F[2],0,0))


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
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = cols, formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('out.mqm.F', 'diff', 'mlod', sep="_"), out.mqm.F)
assign(paste('chr', 'diff', 'mlod', sep="_"), chr)
assign(paste('pos', 'diff', 'mlod', sep="_"), pos)
assign(paste('qtl', 'diff', 'mlod', sep="_"), qtl)
assign(paste('my.formula', 'diff', 'mlod', sep="_"), my.formula)
assign(paste('lodmat.F', 'diff', 'mlod', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time
par(mfrow=c(2,1))
pdf(paste('mqm_timeseries_qtl_', job, '.diff_mlod.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="MLOD")
refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = cols, usec = "mlod", qtl= qtl, method = "hk", keeplodprofile = T)
plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

assign(paste('refqtlmlod ', 'diff', 'mlod', sep="_"), refqtlslod)


mlodeff <- vector("list", length(cols))

for(i in 1:length(mlodeff)) {
  mlodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(cols)-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[cols], '_')

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
pdf(paste('mqm_timeseries_fx_size_mlod_',job,'.diff.pdf', sep=""))
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

assign(paste('mlodeff ', 'diff', 'mlod', sep="_"), mlodeff)
write.csv(mlodeff, file=paste('mlodeff_', 'diff', "_", job, '.csv', sep=""), quote=F, row.names=F)
}

if(length(out.mqm.F) < 1) {
  write("No QTL detected using stepwiseqtl funqtl fxn.", file=paste(dir.job, '/mqm_timeseries_qtl_', job, '.diff', '_mlod.txt' ,sep=""))
}

save.image(file=paste('timeseries_', 'diff', '_cross.object.diff.Rdata', sep=""))

######## Now lets do the same thing all unique qtl at all time points for t1
unique_qtl$chr<-as.numeric(as.character(unique_qtl$chr))
unique_qtl$pos<-as.numeric(as.character(unique_qtl$pos))
qtl<-makeqtl(fg.cr.obj, unique_qtl$chr, unique_qtl$pos, what=c("prob"))
Qs<-paste('Q', 1:length(unique_qtl$pos), sep="")
my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
lodmat.F<-getprofile(fg.cr.obj, qtl =  qtl, pheno.cols = cols, formula = my.formula, method = "hk", verbose = F, tpy="comb")

assign(paste('chr', 'diff', 'all_qtl', sep="_"), chr)
assign(paste('pos', 'diff', 'all_qtl', sep="_"), pos)
assign(paste('qtl', 'diff', 'all_qtl', sep="_"), qtl)
assign(paste('my.formula', 'diff', 'all_qtl', sep="_"), my.formula)
assign(paste('lodmat.F', 'diff', 'all_qtl', sep="_"), lodmat.F)


# Make plots of the QTL LOD profile and the LOD profile over time

par(mfrow=c(1,1))
pdf(paste('mqm_timeseries_qtl_', job, 'diff_all_qtl.pdf', sep=""))
plotprofile(lodmat.F, mval = 8, col=heat.colors(100)[100:1], main="All QTL")
#refqtlslod <- refineqtlF(fg.cr.obj, pheno.cols = cols, usec = "slod", qtl= qtl, method = "hk", keeplodprofile = T)
#plotLodProfile(refqtlslod)
dev.off()
par(mfrow=c(1,1))

#assign(paste('refqtlslod ', 'diff', 'all_qtl', sep="_"), refqtlslod)

klodeff <- vector("list", length(cols))

for(i in 1:length(klodeff)) {
  klodeff[[i]] <- summary(fitqtl(fg.cr.obj, phe=i+(min(cols)-1), qtl=qtl,
                                 method="hk", get.ests=TRUE,
                                 dropone=FALSE))$ests[,1]*c(1,2,2)
}

pname_list<-strsplit(phenames(fg.cr.obj)[cols], '_')

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
pdf(paste('mqm_timeseries_fx_size_all_qtl_', job, ".diff.pdf", sep=""))
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

assign(paste('klodeff', 'diff', 'all_qtl', sep="_"), klodeff)

write.csv(mlodeff, file=paste('klodeff_', 'diff', "_", job, '.csv', sep=""), quote=F, row.names=F)
save.image(file=paste('timeseries_', 'diff', '_cross.object.diff.Rdata', sep=""))