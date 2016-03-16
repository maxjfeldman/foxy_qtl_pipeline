# This is a function for plotting comparisons between the results of treatment specific QTL analysis
# Basically, plot the treatment who's name is first in alphabetical order in orange (strong penalty) or gold (permissive penalty)
# Plot the second treatment (second in alphabetical order) in dark blue (strong penalty) or light blue (permissive penalty)

library(qtl)
library(ggplot2)

args <- commandArgs(TRUE)
# Get phenotype as an argument
trait <- args[1] # Trait
cond1 <- args[2] # Treatment level 1
cond2 <- args[3] # Treatment level 2
outstem <- args[4] # Base directory
dir<-getwd()
dir.stem <- paste(dir, outstem, sep="/")
dir.cond1 <- paste(dir, outstem, trait, cond1, sep="/")
dir.cond2 <- paste(dir, outstem, trait, cond2, sep="/")

pheno1<-paste(cond1, trait, sep=".")
pheno2<-paste(cond2, trait, sep=".")

# Lets load in the data
cr.obj1<-paste('cross.obj_', pheno1, '_raw.Rdata', sep="")
cr.obj2<-paste('cross.obj_', pheno2, '_raw.Rdata', sep="")
load(paste(dir.cond1, cr.obj1, sep="/"))
load(paste(dir.cond2, cr.obj2, sep="/"))

# Prepare a text string
fname<-c('comparison')

# Get a directory to write to
cmpdir<-paste(dir.stem,trait,fname, sep="/")

# What is the y-lim value?
limit.comp.so<-max(get(paste('limit.so', pheno1, sep=".")), get(paste('limit.so', pheno2, sep="."))) + 0.5

pdf(file=paste(cmpdir, '/scanone.qtl.', trait, '.pdf', sep=""))
plot(get(paste('out.so',pheno1, sep=".")), get(paste('out.so', pheno2, sep=".")), col=c("orange", "blue"), main=trait, ylab=c("LOD score"), ylim=c(0,limit.comp.so), bandcol="gray70")
abline(h=get(paste('max.perm.so', pheno1, sep=".")), col=c("orange"))
abline(h=get(paste('max.perm.so', pheno2, sep=".")), col=c("red"))
dev.off()

# If both treatments show QTLs in MQM lite
if (length(get(paste('stepout.a.lite', pheno1, sep="."))) > 0 & length(get(paste('stepout.a.lite', pheno2, sep="."))) > 0) { 
  p1<-get(paste('limit.so', pheno1, sep="."))
  p2<-get(paste('limit.so', pheno2, sep="."))

  if (p1 >= p2) {

    pdf(file=paste(cmpdir, '/mqm.qtl.', trait, '.pdf', sep=""))

    plotLodProfile(get(paste('stepout.a.lite', pheno1, sep=".")), col=c("yellow"), main=trait, showallchr=T)
    if (length(get(paste('stepout.a', pheno1, sep="."))) > 0) {
      plotLodProfile(get(paste('stepout.a', pheno1, sep=".")), col=c("orange"), add=T, showallchr=T)
    }
  
    plotLodProfile(get(paste('stepout.a.lite', pheno2, sep=".")), col=c("light blue"), add=T, showallchr=T)
    if (length(get(paste('stepout.a', pheno2, sep="."))) > 0) {
      plotLodProfile(get(paste('stepout.a', pheno2, sep=".")), col=c("navy"), add=T, showallchr=T)  
    }

    abline(h=get(paste('max.perm.so', pheno1, sep=".")), col=c("orange"))
    abline(h=get(paste('max.perm.so', pheno2, sep=".")), col=c("red"))
    dev.off()
  } 
  
  if (p1 < p2) {
      pdf(file=paste(cmpdir, '/mqm.qtl.', trait, '.pdf', sep=""))  

      plotLodProfile(get(paste('stepout.a.lite', pheno2, sep=".")), col=c("light blue"), main=trait, showallchr=T)
      if (length(get(paste('stepout.a', pheno2, sep="."))) > 0) {
        plotLodProfile(get(paste('stepout.a', pheno2, sep=".")), col=c("navy"), add=T, showallchr=T)  
      }

      plotLodProfile(get(paste('stepout.a.lite', pheno1, sep=".")), col=c("yellow"), add=T, showallchr=T)
      if (length(get(paste('stepout.a', pheno1, sep="."))) > 0) {
        plotLodProfile(get(paste('stepout.a', pheno1, sep=".")), col=c("orange"), add=T, showallchr=T)
      }
      
      abline(h=get(paste('max.perm.so', pheno1, sep=".")), col=c("orange"))
      abline(h=get(paste('max.perm.so', pheno2, sep=".")), col=c("red"))
      dev.off()
  }
}

# If only the first treatment show QTLs in MQM lite
if (length(get(paste('stepout.a.lite', pheno1, sep="."))) > 0 & length(get(paste('stepout.a.lite', pheno2, sep="."))) == 0) { 
  
  pdf(file=paste(cmpdir, '/mqm.qtl.', trait, '.pdf', sep=""))
  plotLodProfile(get(paste('stepout.a.lite', pheno1, sep=".")), col=c("yellow"), main=trait, showallchr=T)
  if (length(get(paste('stepout.a', pheno1, sep="."))) > 0) {
    plotLodProfile(get(paste('stepout.a', pheno1, sep=".")), col=c("orange"), add=T, showallchr=T)
  }
  abline(h=get(paste('max.perm.so', pheno1, sep=".")), col=c("red"))
  dev.off()
}

if (length(get(paste('stepout.a.lite', pheno1, sep="."))) == 0 & length(get(paste('stepout.a.lite', pheno2, sep="."))) > 0) { 
  
  pdf(file=paste(cmpdir, '/mqm.qtl.', trait, '.pdf', sep=""))
  plotLodProfile(get(paste('stepout.a.lite', pheno2, sep=".")), col=c("light blue"), main=trait, showallchr=T)
  if (length(get(paste('stepout.a', pheno2, sep="."))) > 0) {
    plotLodProfile(get(paste('stepout.a', pheno2, sep=".")), col=c("navy"), add=T, showallchr=T)
  }
  abline(h=get(paste('max.perm.so', pheno2, sep=".")), col=c("orange"))
  dev.off()
}


phe1<-get(paste('phevalues', pheno1, sep="."))
phe1<-as.data.frame(phe1)
colnames(phe1)[2]<-c('phenotype')
phe1$cond<-rep(cond1, nrow(phe1))
phe2<-get(paste('phevalues', pheno2, sep="."))
phe2<-as.data.frame(phe2)
colnames(phe2)[2]<-c('phenotype')
phe2$cond<-rep(cond2, nrow(phe2))
cmpphe<-rbind(phe1,phe2)

tt.result<-t.test(phe1[,2],phe2[,2])
pval<-tt.result[3]

pdf(file=paste(cmpdir, '/boxplot.', trait, '.pdf', sep=""))
boxplot(cmpphe$phenotype~cmpphe$cond, col=c("orange", "blue"), main = paste(trait, ' p-value = ', pval, sep=""))
dev.off()

pdf(file=paste(cmpdir, '/histogram.', trait, '.pdf', sep=""))
ggplot(cmpphe, aes(phenotype, fill=cond)) + geom_density(alpha = 0.2) + ggtitle(paste(trait, ' p-value = ', pval, sep=""))
dev.off()