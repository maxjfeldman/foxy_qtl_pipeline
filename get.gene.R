# This is a program to return genes of interest from QTL intervals used in the case of comparisons
# This will need to be re-wired to enable use of multiple organisms
# Currently hardwired for Setaria viridis

library(qtl)

dir<-getwd()
warning(paste('PWD is: ', dir, sep=""))
dir.annotation<-paste(dir, 'annotation', sep = "/")
warning(paste('dir.annotation is: ', dir.annotation, sep=""))
# Read in genome annotation file
# Static path on server
annotation.path<-paste(dir.annotation, 'Sviridis_311_v1.1.gene.gff3', sep="/")
# Path to annotation file is:
warning(paste('Path to annotation is: ', annotation.path, sep=""))

if (file.exists(annotation.path)) {warning("The thing is here!")} else {warning("Can't find it!")}

genes<-read.table(annotation.path, sep="\t")
# Change columns to appropraite names
colnames(genes)<-c('scaffold', 'source', 'type', 'l_pos', 'r_pos', 'score', 'strand', 'phase', 'attr')
# Keep only protein coding genes
genes<-genes[genes$type == 'gene',]

# Read in peptide annotation file
# Static path on server
gene.fxn.path<-paste(dir.annotation, 'Sviridis_311_v1.1.annotation_info.txt', sep="/")
annotation<-read.table(gene.fxn.path, sep="\t", quote="", stringsAsFactors = FALSE, fill=TRUE)
colnames(annotation)<-c('pacid', 'locusName', 'transcriptName', 'peptideName', 'Pfam', 'Panther', 'KOG', 'KEGG', 'KO', 'GO', 'TAIR_ID', 'TAIR_GENE', 'TAIR_ANNOTATION', 'RICE_ID', 'RICE_GENE','RICE_ANNOTATION')

args <- commandArgs(TRUE)
dir.trait<-args[1]
trait<-args[2]

# specify directory path for MQM and Scanone results
path.mqm<-paste(dir.trait, '/mqm.out', sep="")
path.so<-paste(dir.trait, '/scanone.out', sep="")

# specify path of MQM file
mqm.qtl<-paste(dir.trait, '/mqm.out/', 'summary.table.mqm.', trait, ".csv", sep="")
so.qtl<-paste(dir.trait, '/scanone.out/', 'summary.table.so.', trait, ".csv", sep="")

if (file.exists(mqm.qtl)) {
  qtls.mqm<-read.csv(mqm.qtl)
  
  annotation.out<-c()
  for(rw in 1:nrow(qtls.mqm)) {    
    
    # Get genome coordinates from QTL analysis
    q<-qtls.mqm[rw,1]
    qtl.marker<-unlist(strsplit(as.character(q), split="_"))
    qtl.marker[1]<-sub("S", "", qtl.marker[1])
    qtl.marker[1]<-sprintf("%02d", as.numeric(qtl.marker[1]))
    qtl.marker[1]<-paste("Chr", qtl.marker[1], sep="_")
    qtl_chr<-qtl.marker[1]
    qtl_pos<-as.numeric(qtl.marker[2])
    lb<-qtls.mqm[rw,8]
    l.boundary <- unlist(strsplit(as.character(lb), split="_"))
    rb<-qtls.mqm[rw,12]
    r.boundary <- unlist(strsplit(as.character(rb), split="_"))
    # Get the genes within the confidence interval
    interval<-genes[((genes$scaffold == qtl_chr) & ((genes$l_pos >= as.numeric(l.boundary[2])) & (genes$r_pos <= as.numeric(r.boundary[2])))),]
    if (nrow(interval) == 0) {next}
    # Get the gene names within the interval
    sv.gene.id<-as.character(interval[,9])
    sv.gene.id<-sub('ID=', '', sv.gene.id)
    sv.gene.id<-sub(';Name=.*', '', sv.gene.id)
    sv.gene.id<-sub('.v1.1', '', sv.gene.id)
    conf_int_genes<-as.data.frame(sv.gene.id)
    conf_int_genes<-cbind(rep(as.character(q), nrow(conf_int_genes)), conf_int_genes,  interval$l_pos, interval$r_pos)
    colnames(conf_int_genes)<-c('qtl', 'gene', 'l_pos', 'r_pos')
    qtl.out<-c()
    if (length(sv.gene.id) < 1) {next;}
    for(n in 1:length(sv.gene.id)) {
      a<-annotation[annotation$name %in% sv.gene.id[n],]
      if(nrow(a) > 0) {
        #gn<-as.data.frame(rep(conf_int_genes[n,], nrow(a)), ncol=4)
        g<-conf_int_genes[n,]
        gn<-g[rep(seq_len(nrow(g)), nrow(a)),]
      } else {
        a[1,]<-c(rep('unknown', ncol(a)))
        gn<-as.data.frame(rep(conf_int_genes[n,], nrow(a)))
      }
      entry<-cbind(gn, a)
      qtl.out<-rbind(qtl.out, entry)
      
    }
    # Lets sort candidate genes by position to closest marker
    # First find average between L and R boundary of gene
    qtl.out$a_pos<-rowMeans(qtl.out[,c(3,4)])
    # Subtract the position of best marker from position of genes in CI
    qtl.out$diff_pos<-qtl_pos-qtl.out$a_pos
    qtl.out$diff_pos<-abs(qtl.out$diff_pos)
    qtl.out<-qtl.out[order(qtl.out$diff_pos),]
    qtl.out<-qtl.out[,c(1:4,21,22,5:20)]
    
    annotation.out<-rbind(annotation.out, qtl.out)
  }
  
  annotation.file<-paste(path.mqm, "/genes_in_qtl_CI_", trait, ".txt", sep="")
  write.table(annotation.out, file=annotation.file, quote=F, row.names=F,sep = "\t")
    
} else {
  annotation.file<-paste(path.mqm, "/genes_in_qtl_CI_",  trait, ".txt", sep="")
  write("No significant QTL in MQM", file=annotation.file)
}