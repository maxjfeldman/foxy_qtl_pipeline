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
checker<-as.integer(len/5)
if(tester[checker] > 1){
  # If more than >1 treatment (pheno[,2]) and >1 plot/environment (pheno[,3]) calculate heritability
  if (((length(unique(pheno[,2])) > 1 )) & (length(unique(pheno[,3])) > 1)) {
  # Create variables to store values
   H2<-c()
   p2<-c()
   t2<-c()
   e2<-c()
   gxp2<-c()
   gxt2<-c()

   # For each phenotype calculate variance
   for(i in 5:length(colnames(pheno))){
     # Use only RILs with all measurements for each phenotype
     cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
     # Build linear model each cofactor is a random effect
     model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|plot)+(1|id:treatment)+(1|id:plot), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
     # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
     re<-as.numeric(VarCorr(model))
     res<-attr(VarCorr(model), "sc")^2
     # Extract individual components (order will remain the same)
     gxp.var<-re[1]
     gxt.var<-re[2]
     geno.var<-re[3]
     plot.var<-re[4]
     treat.var<-re[5]
     # Total variance is sum of all variances
     tot.var<-sum(re, res)
     # Get proportion of variance
     h<-geno.var/tot.var
     p<-plot.var/tot.var
     t<-treat.var/tot.var
     e<-res/tot.var
     gxp<-gxp.var/tot.var
     gxt<-gxt.var/tot.var
     # Append variables to a vector of variables
     H2<-c(H2,h)
     p2<-c(p2,p)
     t2<-c(t2,t)
     e2<-c(e2,e)
     gxp2<-c(gxp2, gxp)
     gxt2<-c(gxt2, gxt)
   }

   variance<-rbind(H2, t2, p2, gxt2, gxp2, e2)
   colnames(variance)<-colnames(pheno)[5:length(pheno)]
   rownames(variance)<-c('Genotype', 'Treatment', 'Plot', 'G x Treatment', 'G x Plot', 'Error')
   her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
   write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)

   traits<-colnames(variance)
   types<-rownames(variance)

   variance.l<-c()
   t<-c()
   for(i in 1:length(traits)) {
     t<-rep(traits[i], length(types))
     rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
     variance.l<-rbind(variance.l, rf)
   }
   colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
   rownames(variance.l)<-c(1:nrow(variance.l))
   p = ggplot(data=variance.l, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
   p=p + geom_bar(width = 1, stat="identity")
   p=p+facet_grid(facets=. ~ Trait)
   p=p+xlab("Trait")
   # Print out table as.png
   her.plot<-paste(dir.outstem,'heritability_plot.png', sep="/")
   png(filename=her.plot)
   p
   invisible(dev.off())
 
   #### Now partition just genetic variance in each treatment
   i.treat<-unique(pheno$treatment)
 
   for (t in 1:length(i.treat)) {
     # Create variables to store values
     treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
     H2<-c()
     e2<-c()
 
     # For each treatment.phenotype calculate variance
     for(i in 5:length(colnames(treatment.pheno))){
       # Use only RILs with all measurements for each treatment.phenotype
       cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
       # Build linear model each cofactor is a random effect
       model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
       # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
       re<-as.numeric(VarCorr(model))
       res<-attr(VarCorr(model), "sc")^2
       # Extract individual components (order will remain the same)
       geno.var<-re[1]
       # Total variance is sum of all variances
       tot.var<-sum(re, res)
       # Get proportion of variance
       h<-geno.var/tot.var
       e<-res/tot.var
       # Append variables to a vector of variables
       H2<-c(H2,h)
       e2<-c(e2,e)
     }
 
     variance<-rbind(H2, e2)
     colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
     colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
     rownames(variance)<-c('Genotype', 'Error')
     her.table.name<-paste(i.treat[t], "_heritability_table.csv", sep="")
     her.table.treat<-paste(dir.outstem, her.table.name, sep="/")
     write.table(variance, file=her.table.treat, append=F, quote=F, sep=",", row.names=T, col.names=NA)
 
   }
 
  }

  traits<-colnames(pheno)[5:ncol(pheno)]


  # If more than >1 treatment (pheno[,2] and only 1 plot/environment (pheno[,3]) calculate heritability
  if (((length(unique(pheno[,2])) > 1 )) & (length(unique(pheno[,3])) == 1)) {
    # Create variables to store values
    H2<-c()
    t2<-c()
    e2<-c()
    gxt2<-c()
  
    # For each phenotype calculate variance
    for(i in 5:length(colnames(pheno))){
      # Use only RILs with all measurements for each phenotype
      cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.pheno[,4]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      gxt.var<-re[1]
      geno.var<-re[2]
      treat.var<-re[3]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      t<-treat.var/tot.var
      e<-res/tot.var
      gxt<-gxt.var/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      t2<-c(t2,t)
      e2<-c(e2,e)
      gxt2<-c(gxt2, gxt)
    }
  
    variance<-rbind(H2, t2, gxt2, e2)
    colnames(variance)<-colnames(pheno)[5:length(pheno)]
    rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
    her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
    write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  
    traits<-colnames(variance)
    types<-rownames(variance)
  
    variance.l<-c()
    t<-c()
    for(i in 1:length(traits)) {
      t<-rep(traits[i], length(types))
      rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
      variance.l<-rbind(variance.l, rf)
    }
    colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
    rownames(variance.l)<-c(1:nrow(variance.l))
    p = ggplot(data=variance.l, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
    p=p + geom_bar(width = 1, stat="identity")
    p=p+facet_grid(facets=. ~ Trait)
    p=p+xlab("Trait")
    # Print out table as.png
    her.plot<-paste(dir.outstem,'heritability_plot.png', sep="/")
    png(filename=her.plot)
    p
    invisible(dev.off())
  
  
    #### Now partition just genetic variance in each treatment
    i.treat<-unique(pheno$treatment)
  
    for (t in 1:length(i.treat)) {
      # Create variables to store values
      treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
      H2<-c()
      e2<-c()
    
      # For each treatment.phenotype calculate variance
      for(i in 5:length(colnames(treatment.pheno))){
        # Use only RILs with all measurements for each treatment.phenotype
        cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:3,i)]
        # Build linear model each cofactor is a random effect
        model<-lmer(cc.treatment.pheno[,4]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
        # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
        re<-as.numeric(VarCorr(model))
        res<-attr(VarCorr(model), "sc")^2
        # Extract individual components (order will remain the same)
        geno.var<-re[1]
        # Total variance is sum of all variances
        tot.var<-sum(re, res)
        # Get proportion of variance
        h<-geno.var/tot.var
        e<-res/tot.var
        # Append variables to a vector of variables
        H2<-c(H2,h)
        e2<-c(e2,e)
      }
    
      variance<-rbind(H2, e2)
      colnames(treatment.pheno)[5:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[5:length(treatment.pheno)], sep="_")
      colnames(variance)<-colnames(treatment.pheno)[5:length(treatment.pheno)]
      rownames(variance)<-c('Genotype', 'Error')
      her.table.name<-paste(i.treat[t], "_heritability_table.csv", sep="")
      her.table.treat<-paste(dir.outstem, her.table.name, sep="/")
      write.table(variance, file=her.table.treat, append=F, quote=F, sep=",", row.names=T, col.names=NA)
    
    }
  
  }


  # If 1 treatment (pheno[,2]) and >1 plot/environment (pheno[,3]) calculate heritability
  if (((length(unique(pheno[,2])) == 1 )) & (length(unique(pheno[,3])) > 1)) {
    # Create variables to store values
    H2<-c()
    p2<-c()
    e2<-c()
    gxp2<-c()
  
    # For each phenotype calculate variance
    for(i in 5:length(colnames(pheno))){
      # Use only RILs with all measurements for each phenotype
      cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.pheno[,4]~(1|id)+(1|plot)+(1|id:plot), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      gxp.var<-re[1]
      geno.var<-re[2]
      plot.var<-re[3]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      p<-plot.var/tot.var
      e<-res/tot.var
      gxp<-gxp.var/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      p2<-c(p2,p)
      e2<-c(e2,e)
      gxp2<-c(gxp2, gxp)
    }
  
    variance<-rbind(H2, p2, gxp2, e2)
    colnames(variance)<-colnames(pheno)[5:length(pheno)]
    rownames(variance)<-c('Genotype', 'Plot', 'G x Plot', 'Error')
    her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
    write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  
    traits<-colnames(variance)
    types<-rownames(variance)
  
    variance.l<-c()
    t<-c()
    for(i in 1:length(traits)) {
      t<-rep(traits[i], length(types))
      rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
      variance.l<-rbind(variance.l, rf)
    }
    colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
    rownames(variance.l)<-c(1:nrow(variance.l))
    p = ggplot(data=variance.l, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
    p=p + geom_bar(width = 1, stat="identity")
    p=p+facet_grid(facets=. ~ Trait)
    p=p+xlab("Trait")
    # Print out table as.png
    her.plot<-paste(dir.outstem,'heritability_plot.png', sep="/")
    png(filename=her.plot)
    p
    invisible(dev.off())
  
    #### Now partition just genetic
      # Create variables to store values
      her.pheno<-pheno
      H2<-c()
      e2<-c()
    
      # For each her.phenotype calculate variance
      for(i in 5:length(colnames(her.pheno))){
        # Use only RILs with all measurements for each her.phenotype
        cc.her.pheno<-her.pheno[complete.cases(her.pheno[,i]),c(1:3,i)]
        # Build linear model each cofactor is a random effect
        model<-lmer(cc.her.pheno[,4]~(1|id), data=cc.her.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
        # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
        re<-as.numeric(VarCorr(model))
        res<-attr(VarCorr(model), "sc")^2
        # Extract individual components (order will remain the same)
        geno.var<-re[1]
        # Total variance is sum of all variances
        tot.var<-sum(re, res)
        # Get proportion of variance
        h<-geno.var/tot.var
        e<-res/tot.var
        # Append variables to a vector of variables
        H2<-c(H2,h)
        e2<-c(e2,e)
      }
    
      variance<-rbind(H2, e2)
      colnames(variance)<-colnames(her.pheno)[5:length(her.pheno)]
      rownames(variance)<-c('Genotype', 'Error')
      her.table.name<-paste("genotype", "_heritability_table.csv", sep="")
      her.table.treat<-paste(dir.outstem, her.table.name, sep="/")
      write.table(variance, file=her.table.treat, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  
  }


  # If 1 treatment (pheno[,2]) and 1 plot/environment (pheno[,3]) calculate heritability
  if (((length(unique(pheno[,2])) == 1 )) & (length(unique(pheno[,3])) == 1)) {
    # Create variables to store values
    H2<-c()
    e2<-c()
  
    # For each phenotype calculate variance
    for(i in 5:length(colnames(pheno))){
      # Use only RILs with all measurements for each phenotype
      cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:3,i)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.pheno[,4]~(1|id), data=cc.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
  
    variance<-rbind(H2, e2)
    colnames(variance)<-colnames(pheno)[5:length(pheno)]
    rownames(variance)<-c('Genotype', 'Error')
    her.table<-paste(dir.outstem,"heritability_table.csv", sep="/")
    write.table(variance, file=her.table, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  
    traits<-colnames(variance)
    types<-rownames(variance)
  
    variance.l<-c()
    t<-c()
    for(i in 1:length(traits)) {
      t<-rep(traits[i], length(types))
      rf<-cbind(t, types, data.frame(variance[,i], stringsAsFactors=F))
      variance.l<-rbind(variance.l, rf)
    }
    colnames(variance.l)<-c('Trait', 'Type', 'Proportion_explained')
    rownames(variance.l)<-c(1:nrow(variance.l))
    p = ggplot(data=variance.l, aes(x=factor(1), y=Proportion_explained, fill = factor(Type)))
    p=p + geom_bar(width = 1, stat="identity")
    p=p+facet_grid(facets=. ~ Trait)
    p=p+xlab("Trait")
    # Print out table as.png
    her.plot<-paste(dir.outstem,'heritability_plot.png', sep="/")
    png(filename=her.plot)
    p
    invisible(dev.off())
  
    #### Now partition just genetic
    # Create variables to store values
    her.pheno<-pheno
    H2<-c()
    e2<-c()
  
    # For each her.phenotype calculate variance
    for(i in 5:length(colnames(her.pheno))){
      # Use only RILs with all measurements for each her.phenotype
      cc.her.pheno<-her.pheno[complete.cases(her.pheno[,i]),c(1:3,i)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.her.pheno[,4]~(1|id), data=cc.her.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
  
    variance<-rbind(H2, e2)
    colnames(variance)<-colnames(her.pheno)[5:length(her.pheno)]
    rownames(variance)<-c('Genotype', 'Error')
    her.table.name<-paste("genotype", "_heritability_table.csv", sep="")
    her.table.treat<-paste(dir.outstem, her.table.name, sep="/")
    write.table(variance, file=her.table.treat, append=F, quote=F, sep=",", row.names=T, col.names=NA)
  }

}
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
      tp$diff<-tp[,2]-tp[,3]
      colnames(tp)[4]<-traits[k]
      diff<-cbind(diff, tp[,traits[k]])
      colnames(diff)[k+1]<-paste(traits[k], paste(treatments[j], treatments[j+1], sep="-"), sep=".")
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