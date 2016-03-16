# This is a Rscript to make a file to input phenotypes into R/qtl using the "csvs" input option
# It also makes a RIL cross object, removes individuals with greater than 10% of their genotypes missing
# Calculates probability of correct genotype call

library(qtl)

args <- commandArgs(TRUE)
map_name<-args[1]
outstem<-args[2]
cross_type<-args[3]

dir<-getwd()
dir.outstem<-paste(dir, outstem, sep="/")
pheno.raw<-paste(dir.outstem, "qtl.phenotypes.raw.csv", sep="/")

if (file.exists(pheno.raw)) { warning("The qtl.phenotypes.raw.csv file already exists!") } else {
  # The genetic map is specified as a command line argument
  if (file.exists(map_name)) { warning("The map_name file exists!") }
  map <- read.csv(map_name)
  # Get the names of the RILs
  names <- map[2:nrow(map),1]
  # Put the names in a data.frame
  pheno.r <- as.data.frame(names)
  # Give the column which contains the RIL ID numbers a column header caled 'id'
  colnames(pheno.r)[1] <- c("id")
  # Find all phenotype files (all phenotype files must have the suffix ".phe.csv"
  setwd(dir.outstem)
  files.r <- list.files(pattern=c("*.phe.csv"))
  # For each file in the current directory with the "*.phe.csv" suffix, read the file in
  for(i in 1:length(files.r)) {
    label <- files.r[i]; 
    temp <- read.csv(file=files.r[i]); 
     # Write all phenotypes to a single data.frame
    pheno.r <- merge(pheno.r, temp, by.x="id", by.y="id", all.x=T)
    rm(label, temp, i);
  }
  # Replace all empty entries (in this case ".") with NA
  pheno.r[pheno.r == "."] <- NA
  # Write the phenotype file to current directory
  write.csv(pheno.r, row.names=F, file=paste(dir.outstem, "qtl.phenotypes.raw.csv", sep="/"), quote=F)
}

infile.r<-paste(dir.outstem, "qtl.phenotypes.raw.csv", sep="/")
outfile.r<-paste(dir.outstem, "cross.object.raw.Rdata", sep="/")
map_name<-paste(dir, map_name, sep="/")

if (file.exists(outfile.r)) { warning("A cross object for this phenotype set already exists!") } else {
  # Build cross object from mapping and phenotype file 
  
  if (length(grep('2013setariamapJGI.csv', map_name)) > 0) {
    cross.obj <- read.cross(format=c("csvs"), genfile=map_name, phefile=infile.r, crosstype=cross_type, estimate.map=T)
  }
  if (length(grep('2013setariamapJGI.csv', map_name)) == 0) {
    cross.obj <- read.cross(format=c("csvs"), genfile=map_name, phefile=infile.r, genotypes=c("AA","AB","BB"), crosstype=cross_type, estimate.map=F)
  }
  
  # Filter out individuals with greater than 10% missing genotype calls
  prop.called<-ntyped(cross.obj, what=c("ind"))
  less_than_10_per_missing<-names(prop.called[prop.called > (.9*sum(nmar(cross.obj)))])
  cross.obj<-subset(cross.obj, ind=less_than_10_per_missing)
  # Use jittermap to set very, very similar markers apart 
  cross.obj<-jittermap(cross.obj)
  # Calculate underlying genotype probabilities at markers
  fg.cr.obj<-calc.genoprob(cross.obj, error.prob=0.0005, map.function=c("kosambi"))
  #fg.cr.obj<-sim.geno(fg.cr.obj, error.prob=0.0005, map.function=c("kosambi"))

  save.image(file=outfile.r)
  
  # Get number of phenotypes for looping
  numphe.r<-nphe(fg.cr.obj)
  namephe.r<-colnames(fg.cr.obj$pheno)[2:length(colnames(fg.cr.obj$pheno))]
  phefile.r<-paste(dir.outstem, "phe.r.txt", sep="/")
  write(numphe.r, file=phefile.r)
  write(namephe.r, file=phefile.r, append=TRUE)
}

save.image(file=outfile.r)
