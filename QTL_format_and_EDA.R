library(lme4)
library(ggplot2)
library(lattice)
args <- commandArgs(TRUE)

# First argument is phenotype .csv file ($input in perl/python program)
trait.file<-args[1]
dir<-getwd()
# Second argument is output directory ($outstem in perl/python program)
outstem<-args[2]
# Get path to write results to
dir.outstem<-paste(dir, outstem, sep="/")
# Read file and re-format keeping only interesting columns
pheno<-read.csv(trait.file, na.strings='NA')
year<-unique(pheno[,3])
exp<-unique(pheno[,2])
pheno<-pheno[,c(7,4,5,8:length(colnames(pheno)))]
colnames(pheno)[c(1,2)]<-c("id", "treatment")
pheno[pheno == "."] <- NA
colnames(pheno)[5:ncol(pheno)]<-paste(colnames(pheno)[5:ncol(pheno)] , exp, sep="_")

tester<-sort(table(pheno$id), decreasing=T)
len<-length(tester)


traits<-colnames(pheno)[5:ncol(pheno)]
# How many columns? (i.e. phenotypes)
n.cols<-length(colnames(pheno))

# Make sure all phenotypes are numeric (i.e. not factor)
for (q in 5:n.cols) {
   pheno[,q]<-as.numeric(as.character(pheno[,q]))
}

# Get average value of phenotype across plots
pheno.ag<-aggregate(pheno[,5:n.cols], by=list(pheno$id, pheno$treatment), mean, na.action = na.pass, na.rm=TRUE)
colnames(pheno.ag)<-c("id", "treatment", colnames(pheno)[5:length(colnames(pheno))])
treatments<-levels(pheno.ag$treatment)

if (length(treatments) < 2) {
  pheno.ag.final<-pheno.ag[,-2]
  out.name<-paste('single_treatment.phe.csv', sep=".")
  out.file<-paste(dir.outstem, out.name, sep="/")
  write.table(pheno.ag.final, file=out.file, append=F, quote=F, sep=",", row.names=F)
}

# If more than 1 treatment is present calculate difference between treatments. 
if (length(treatments) > 1) {
  diff.pheno<-unique(pheno.ag[,1])
  diff.pheno<-as.data.frame(diff.pheno)
  colnames(diff.pheno)[1]<-c('id')
  for(i in 1:length(treatments)) {
    t<-treatments[i]
    n.cols<-length(colnames(pheno.ag))
    c.names<-paste(t,colnames(pheno.ag)[3:n.cols], sep=".")
    pheno.ag.final<-pheno.ag[pheno.ag$treatment == t,c(1,3:n.cols)]
    n.cols<-length(colnames(pheno.ag.final))
    colnames(pheno.ag.final)[2:n.cols]<-c.names
    out.name<-paste(t, 'phe.csv', sep=".")
    out.file<-paste(dir.outstem, out.name, sep="/")
    write.table(pheno.ag.final, file=out.file, append=F, quote=F, sep=",", row.names=F)
    diff.pheno<-merge(diff.pheno, pheno.ag.final, by = c('id'))
  }
  # Calculate the difference between phenotypes and write to a dataframe 
  diff<-as.data.frame(diff.pheno[,1])
  colnames(diff)[1]<-c('id')
  # First loop through treatments
  for(j in 1:(length(treatments)-1)) {
    # Subset the dataframe to include only 2 treatments at once
    t1<-diff.pheno[,c(1, grep(treatments[j], colnames(diff.pheno)))]
    t2<-diff.pheno[,c(1, grep(treatments[j+1], colnames(diff.pheno)))]
    # Next loop through traits, You have two treatments (t1 & t2), now for each phenotype subtract values between treatments (tp1 - tp2)
    for(k in 1:length(traits)) {
      # Subset the treatment specific dataframe by phenotype (trait)
      tp1<-t1[,c(1, grep(traits[k], colnames(t1)))]
      tp2<-t2[,c(1, grep(traits[k], colnames(t2)))]
      tp<-merge(tp1, tp2, by = "id")
      # Calculate the absolute difference
      tp$diff<-tp[,2]-tp[,3]
      # Calculate the relative difference as defined by wikipedia
      tp$rel_diff<-abs(tp[,2] - tp[,3])/(max(abs(tp[,2]), abs(tp[,3]), na.rm=TRUE))
      # Calculate the ratio between the two, first make sure the second trait != 0
      col3<-colnames(tp)[3]
      # Cannot divide by zero so set any values equal to zero to some value close to zero
      tp[col3 == 0, 3] <-c(0.001)
      tp$ratio <- abs(tp[,2]/tp[,3])
      colnames(tp)[4]<-paste(traits[k], 'diff', sep="=")
      colnames(tp)[4]<-paste(colnames(tp)[4], paste(treatments[j], treatments[j+1], sep="-"), sep=".")
      colnames(tp)[5]<-paste(traits[k], 'rel_diff', sep="=")
      colnames(tp)[5]<-paste(colnames(tp)[5], paste(treatments[j], treatments[j+1], sep="-"), sep=".")
      colnames(tp)[6]<-paste(traits[k], 'ratio', sep="=")
      colnames(tp)[6]<-paste(colnames(tp)[6], paste(treatments[j], treatments[j+1], sep="-"), sep=".")
      diff<-cbind(diff, tp[,4:6])
    }
  }

  # Write .csv table
  diff.file.name<-paste(dir.outstem,"diff_phenotype_by_treatment.csv", sep="/")
  write.table(diff, file=diff.file.name, append=F, quote=F, sep=",", row.names=F)  
}

# Save workspace
image.name<-paste(dir.outstem, 'format_and_EDA.Rdata', sep="/")
save.image(file=image.name)

# If >1 phenotype evaluate phenotype correlation using PCC
if (length(traits) > 1) {
    
cc.pheno.ag<-pheno.ag[complete.cases(pheno.ag),]
cor.mat<-cor(as.matrix(cc.pheno.ag[,3:length(colnames(cc.pheno.ag))]))
rownames(cor.mat)<-colnames(cor.mat)
trait_cor.table<-paste(dir.outstem,'trait_correlation_table.csv', sep="/")
write.table(cor.mat, file=trait_cor.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)

trait_cor.plot<-paste(dir.outstem,'trait_correlation_plot.png', sep="/")
png(filename=trait_cor.plot)
levelplot(cor.mat, scales=list(x=list(rot=90)), ylab="", xlab="", at=seq(-1,1,0.1))
invisible(dev.off())

}
save.image(image.name)