# This Rscript performs scanone QTL analysis and writes plots and results to an already existing phenotype/scanone directory
# This script is run if the user is not doing a comparison. It also runs analysis on the numerical difference between traits.

args <- commandArgs(TRUE)
library(qtl)
# Get phenotype as an argument
i <- args[1]
outstem <- args[2]
i<-as.numeric(i)
trait<-args[3]
workspace<-args[4]
trait_model<-args[5]
qtl_method<-args[6]

dir<-getwd()
dir.stem <-paste(dir, outstem, sep="/")
load(paste(dir.stem, workspace, sep="/"))
map<-read.csv(map_name)
pname<-phenames(fg.cr.obj)[i]

# Are you analyzing a numerical difference trait? 
if (workspace == "cross.object.diff.Rdata") {
  dirpname<-paste(dir.stem, trait,'comparison', sep="/")
  dirscanone<-paste(dirpname, "scanone.out", sep="/")
  dirmqm<-paste(dirpname, "mqm.out", sep="/")
  session_image_name<-paste('cross.obj_', pname, "_diff.Rdata", sep="")

}

# Are you analyzing a single non-comparison trait?
if (workspace == "cross.object.raw.Rdata") {
  dirpname<-paste(dir.stem, pname, sep="/")
  dirscanone<-paste(dir.stem,pname,"scanone.out", sep="/")
  dirmqm<-paste(dir.stem,pname,"mqm.out", sep="/")
  session_image_name<-paste('cross.obj_', pname, "_raw.Rdata", sep="")

}

# Extract phenotype values, sort them then write them out to trait folder
phevalues<-fg.cr.obj$phe[,c(1,i)]
phevalues<-phevalues[complete.cases(phevalues),]
phevalues<-phevalues[order(phevalues[,2]),]
assign(paste("phevalues",pname,sep="."),phevalues)
write.csv(phevalues, file=paste(dirpname,"/ordered_phe_vals_", pname, ".csv", sep=""), quote=F, row.names=F)

# Perform scanone QTL analysis on a fill.geno cross object using method "hk" and default paramters
out.so <- scanone(fg.cr.obj, pheno.col=i, method=qtl_method, model=trait_model)
assign(paste("out.so",pname,sep="."), out.so)
# Perform significance test using permutation (1000)
operm.so<-scanone(fg.cr.obj, pheno.col=i, method=qtl_method, n.perm=1000, model=trait_model)
assign(paste("operm.so",pname,sep="."), operm.so)
# Get permutation significance threshold (p < 0.05)
max.perm.so<-max(summary(operm.so, alpha=0.05))
assign(paste("max.perm.so",pname,sep="."), max.perm.so)
# Find most signifiant QTL peak
max.out.so<-max(out.so)
max.lod.so<-max.out.so$lod
# Determine range of QTL plot
limit.so<-max(max.lod.so, max.perm.so)+0.5
assign(paste("limit.so",pname,sep="."), limit.so)

# Write significant QTLs to a csv file
qtl.table.so<-summary(out.so, perms=operm.so, alpha=0.05, pvalues=TRUE)
colnames(qtl.table.so)[1]<-c('chr')
table.name<-paste("qtl.table.so", pname, "csv", sep=".")
write.csv(as.data.frame(qtl.table.so), file=paste(dirscanone,table.name, sep="/"), quote=F, row.names=T)
assign(paste("qtl.table.so", pname, sep="."), qtl.table.so)

# Print results to file and include permutation result
pdf(file=paste(dirscanone,"/scanone.qtls.",pname, ".pdf", sep=""))
plot(out.so, ylim=c(0,limit.so), main=pname, bandcol="gray70")
abline(h=max.perm.so, col="red")
dev.off()

pdf(file=paste(dirpname,"/histogram_of_", pname, ".pdf", sep=""))
hist(phevalues[,2], col=c("green"), breaks=20, main=pname)
dev.off()

save.image(file=paste(dirpname,  session_image_name, sep="/" ))

# Do scantwo to look for QTL interactions/epistasis
out.s2<-scantwo(fg.cr.obj, pheno.col=i, method=qtl_method)
assign(paste('out.s2', pname, sep="."), out.s2)
pdf(file=paste(dirscanone,"/scantwo.qtls.",pname, ".pdf", sep=""))
plot(out.s2, main=pname)
dev.off()

# Do scantwo permutation analysis
operm.s2<-scantwo(fg.cr.obj, pheno.col=i, method=qtl_method, n.perm=100)
assign(paste('operm.s2', pname, sep="."), operm.s2)

# Obtain a permissive set of penalties during initial 'stepwise' QTL model search, we will use these as co-factors in MQM
# pen_lite are the permissive penalties
# pen_heavy are the ones used for significance testing

# If looking for epistatic intereactions you'll use this part (remove #s)
# pen_lite<-calc.penalties(operm.s2, alpha=0.25)
# pen_heavy<-calc.penalties(operm.s2, alpha=0.05)

# Right now we are just using an additive model
pen_lite<-c(max(summary(operm.so, alpha=0.25)), 0, 0)
pen_heavy<-c(max.perm.so, 0, 0)
assign(paste('pen_lite', pname, sep="."), pen_lite)
assign(paste('pen_heavy', pname, sep="."), pen_heavy)

# Get list of loci which interact significantly
sig.s2.hk<-summary(out.s2, perms=operm.s2, pvalues=T, alphas=c(0.05, 0.05, 0, 0.05, 0.05))
assign(paste('sig.s2.hk', pname, sep="."), sig.s2.hk)
write.csv(sig.s2.hk, file=paste(dirscanone,"/sig_interaction_qtl.", pname, ".csv", sep=""), quote=F)

best_markers<-summary(out.so)
best_markers<-best_markers[order(best_markers[,3], decreasing=T),]
assign(paste("best_markers", pname, sep="."), best_markers)
marker1<-find.marker(fg.cr.obj, chr=best_markers[1,1], pos=best_markers[1,2])
marker2<-find.marker(fg.cr.obj, chr=best_markers[2,1], pos=best_markers[2,2])
assign(paste("marker1", pname, sep="."), marker1)
assign(paste("marker2", pname, sep="."), marker2)

pdf(file=paste(dirscanone,"/QTL.interactionplot.",pname, ".pdf", sep=""))
effectplot(fg.cr.obj, pheno.col=i, mname1 = marker1, mname2 = marker2)
dev.off()

# Program throws error in effectscan() fxn if sim.geno() not called
fg.cr.obj_sim<-sim.geno(fg.cr.obj)

save.image(file=paste(dirpname, session_image_name, sep="/" ))

if (nrow(qtl.table.so) > 0 ) {
    markers<-row.names(qtl.table.so)
    small.map<-map[map$id %in% phevalues$id,c("id", markers)]
    small.map<-merge(small.map, phevalues, by=c('id'))
    small.map<-small.map[order(small.map[,pname]),]
    write.csv(small.map, file=paste(dirscanone,"/scanone_phe_vals_and_allele_counts_", pname, ".csv", sep=""), quote=F)
    assign(paste('scanone_phe_vals_and_allele_counts', pname, sep="."), small.map)
    
    # Sort list of qtl identified by scanone, and extract chromosome and position
    qtl.table.so<-qtl.table.so[order(qtl.table.so[,3], decreasing=T),]
    chrs<-qtl.table.so[qtl.table.so$lod > max.perm.so,1]
    pos<-qtl.table.so[qtl.table.so$lod > max.perm.so,2]
    qtl.markers<-rownames(qtl.table.so[qtl.table.so$lod > max.perm.so,])
    
    # Get effect size of scanone QTL
    pdf(file=paste(dirscanone,"/so.effect_size.",pname, ".pdf", sep=""))
    fx.so<-effectscan(fg.cr.obj_sim, pheno.col=i, get.se=T, draw=T)
    dev.off()
    
    # Write out the effect size for each marker as a .csv file
    write.csv(fx.so, file=paste(dirscanone,"/allele_effects_all_markers.", pname, ".csv", sep=""), quote=F)
    assign(paste('allele_effects_all_markers', pname, sep="."), fx.so)
    
    fx.qtl.marker<-fx.so[rownames(fx.so) %in% qtl.markers,]
    # Write out the effect size for qtl markers as a .csv file
    write.csv(fx.qtl.marker, file=paste(dirscanone,"/allele_effects_qtl_markers.", pname, ".csv", sep=""), quote=F)
    assign(paste('allele_effects_qtl_markers', pname, sep="."), fx.qtl.marker)
    
    # Make a plot of the PXG
    pdf(file=paste(dirscanone,"/so.phenotypeXgenotype.",pname, ".pdf", sep=""))
    plotPXG(fg.cr.obj, qtl.markers, pheno.col=i)
    dev.off()
    
    # Build a make QTL model based upon the information from scanone
    qtl<-makeqtl(fg.cr.obj, chrs, pos, what=c("prob"))
    
    # Make a formula that specifies this model
    Qs<-paste('Q', 1:length(qtl.markers), sep="")
    my.formula<-as.formula(paste("y~", paste(Qs, collapse="+")))
    assign(paste('so.formula', pname, sep="."), my.formula)
    
    # Fine-tune QTL location
    revqtl<-refineqtl(fg.cr.obj, pheno.col=i, qtl=qtl, formula = my.formula, method=qtl_method)
    assign(paste('so.revqtl', pname, sep="."), revqtl)
    
    # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
    out.fitqtl<-fitqtl(fg.cr.obj, pheno.col=i, qtl=revqtl, formula=my.formula, method=qtl_method, get.ests=T)
    assign(paste('so.out.fitqtl', pname, sep="."), out.fitqtl)
    
    # Proportioning of variance for scanone results (full model)
    full.mdl.var<-as.data.frame(out.fitqtl$result.full[,2])
    full.mdl.var$prop.variance<-c(100*(full.mdl.var[1,1]/full.mdl.var[3,1]), 100*(full.mdl.var[2,1]/full.mdl.var[3,1]), 100)
    colnames(full.mdl.var)<-c("variance", "prop.variance")
    assign(paste('so.full.mdl.var', pname, sep="."), full.mdl.var)
    write.csv(full.mdl.var, file=paste(dirscanone,"/prop.var.full.mdl.", pname, ".csv", sep=""), quote=F)

    if (length(qtl.markers) > 1) {
      
        # Proportioning of variance for scanone results (individual QTL as proportion)
        dropone.mdl.var<-as.data.frame(out.fitqtl$result.drop[,2])
        dropone.mdl.var<-rbind(dropone.mdl.var, full.mdl.var[3,1])
      
        # Partition variance within QTL
        dropone.mdl.var$prop.variance<-c(100*(dropone.mdl.var[1:(nrow(dropone.mdl.var)-1),1]/full.mdl.var[3,1]), 100)
        colnames(dropone.mdl.var)<-c("variance", "prop.variance")
        rownames(dropone.mdl.var)<-c(qtl.markers, "total")
        dropone.residual<-c(dropone.mdl.var[nrow(dropone.mdl.var),1]-sum(dropone.mdl.var[(1:nrow(dropone.mdl.var)-1),1]), dropone.mdl.var[nrow(dropone.mdl.var),2]-sum(dropone.mdl.var[(1:nrow(dropone.mdl.var)-1),2]))
        dropone.mdl.var<-rbind(dropone.mdl.var[1:(nrow(dropone.mdl.var)-1),], dropone.residual, dropone.mdl.var[nrow(dropone.mdl.var),])
        rownames(dropone.mdl.var)[(nrow(dropone.mdl.var)-1)]<-c('residual')
        assign(paste('so.dropone.mdl.var', pname, sep="."), dropone.mdl.var)
        write.csv(dropone.mdl.var, file=paste(dirscanone,"/prop.var.dropone.mdl.", pname, ".csv", sep=""), quote=F)
        
        # Get confidence interval (95%) for each QTL:
        n.QTL<-revqtl$n.qtl
        CI_L<-c()
        CI_R<-c()
        for(q in 1:n.QTL){
          lod_int<-lodint(revqtl, qtl.index=q)
          lod_int$marker<-rownames(lod_int)
          L<-as.data.frame(lod_int[1,c(4,1,2,3)])
          colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
          R<-as.data.frame(lod_int[3,c(4,1,2,3)])
          colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
          CI_L<-rbind(CI_L, L)
          CI_R<-rbind(CI_R, R)
        }
        CI<-cbind(CI_L, CI_R)
        assign(paste('so.confidence_interval', pname, sep="."), CI)
        write.csv(CI, file=paste(dirscanone,"/so.confidence_interval.", pname, ".csv", sep=""), quote=F)
        
        # Lets get addative effect sizes from the fitqtl model of scanone qtl
        fx_size_fitqtl<-as.data.frame(summary(out.fitqtl)[3])
        fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
        fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
    
        # Get proportion of variance per marker
        marker.prop.var<-as.data.frame(dropone.mdl.var[1:(nrow(dropone.mdl.var)-2),2])
    
        # Build a table that summarizes the results
        summary.table.so<-cbind(qtl.table.so[,1:ncol(qtl.table.so)], marker.prop.var, fx_size.a, fx_size.se, CI)
        colnames(summary.table.so)<-c('chr', 'pos', 'lod', 'pval', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
        assign(paste('summary.table.so', pname, sep="."), summary.table.so)
        write.csv(summary.table.so, file=paste(dirscanone,"/summary.table.so.", pname, ".csv", sep=""), quote=F)
    }

    if (length(qtl.markers) == 1) {
        # Full model incorporates the only QTL dropone analysis is empty
        write("No paritioning of variance among QTL if only 1 detected.", file=paste(dirscanone,"/prop.var.dropone.mdl.", pname, ".txt", sep=""))
        
        # Get confidence interval (95%) for the QTL:
        lod_int<-lodint(revqtl)
        lod_int$marker<-rownames(lod_int)
        L<-as.data.frame(lod_int[1,c(4,1,2,3)])
        colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
        R<-as.data.frame(lod_int[3,c(4,1,2,3)])
        colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
        CI<-cbind(L, R)
        assign(paste('so.confidence_interval', pname, sep="."), CI)
        write.csv(CI, file=paste(dirscanone,"/so.confidence_interval.", pname, ".csv", sep=""), quote=F)
        
        # Lets get addative effect sizes from the fitqtl model of scanone qtl
        fx_size_fitqtl<-as.data.frame(summary(out.fitqtl)[2])
        fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
        fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
        
        # Build a table that summarizes the results
        summary.table.so<-cbind(qtl.table.so[,1:ncol(qtl.table.so)], as.data.frame(full.mdl.var[1,2]), fx_size.a, fx_size.se, CI)
        colnames(summary.table.so)<-c('chr', 'pos', 'lod', 'pval', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
        assign(paste('summary.table.so', pname, sep="."), summary.table.so)
        write.csv(summary.table.so, file=paste(dirscanone,"/summary.table.so.", pname, ".csv", sep=""), quote=F)      
    }
    

####################           
# Stepwise MQM analysis
####################

    # If trait model is not normal or binary set to normal
    if (trait_model != 'normal' | trait_model != 'binary') {trait_model<-c("normal")}

    # Perform positive stepwise forward and reverse model selection given formula, fine tuned QTL location, and penalites
    # This is an additive model, run with two different penalties that are specified above
  
    stepout.a<-stepwiseqtl(fg.cr.obj, pheno.col=i, qtl=revqtl, formula=my.formula, method=qtl_method, penalties=pen_heavy, max.qtl=25, scan.pairs=T, additive.only=T, model=trait_model)
    assign(paste('stepout.a', pname, sep="."), stepout.a)
    stepout.a.lite<-stepwiseqtl(fg.cr.obj, pheno.col=i, qtl=revqtl, formula=my.formula, method=qtl_method, penalties=pen_lite, max.qtl=25, scan.pairs=T, additive.only=T, model=trait_model)
    assign(paste('stepout.a.lite', pname, sep="."), stepout.a.lite)
    save.image(paste(dirpname, session_image_name, sep="/"))

    if(length(stepout.a) > 0) {
      
      pdf(file=paste(dirmqm,"/mqm.LOD.profile.",pname, ".pdf", sep=""))
      plotLodProfile(stepout.a.lite, showallchr=T, col= "grey", main=pname)
      plotLodProfile(stepout.a, add=T, showallchr=T, col = "black")
      abline(h=max.perm.so, col="red")
      dev.off()  
      
      assign(paste('stepout.a', pname, sep="."), stepout.a)
      mqm.markers<-c()
      mqm.markers$chr<-stepout.a$chr
      mqm.markers$pos<-stepout.a$pos
      mqm.markers<-as.data.frame(mqm.markers)
    
      m.names<-c()
      for (m in 1:nrow(mqm.markers)) {
         m.name<-find.marker(fg.cr.obj, chr=mqm.markers[m,1], pos=mqm.markers[m,2])
         m.names<-c(m.names, m.name)
      }
    
    # Get lod scores for each significant marker
    mqm.lod.profile<-attr(stepout.a, "lodprofile")
    lods<-c()
    for(l in 1:length(mqm.lod.profile)) {
        d<-as.data.frame(mqm.lod.profile[l])
        colnames(d)<-c('chr', 'pos', 'lod')
        marker_name<-m.names[l]
        if(length(d[marker_name %in% rownames(d),3]) > 0)  {
            lod<-d[rownames(d) == marker_name,3]
            lods<-c(lods, lod)
        }
        if(length(d[marker_name %in% rownames(d),3]) == 0) {
            lod<-max(d[,3])
            lods<-c(lods, lod)
        }
    }
     # Add in the names of the markers and LOD scores
    rownames(mqm.markers)<-m.names
  
    mqm.markers$lod<-lods
    
    qtl.mqm<-makeqtl(fg.cr.obj, mqm.markers$chr, mqm.markers$pos, what=c("prob"))

    # Make a formula that specifies this model
    Qs.mqm<-paste('Q', 1:length(m.names), sep="")
    my.formula.mqm<-as.formula(paste("y~", paste(Qs.mqm, collapse="+")))
    assign(paste('mqm.formula', pname, sep="."), my.formula.mqm)

    # Fine-tune QTL location
    revqtl.mqm<-refineqtl(fg.cr.obj, pheno.col=i, qtl=qtl.mqm, formula = my.formula.mqm, method=qtl_method)
    assign(paste('mqm.revqtl', pname, sep="."), revqtl.mqm)

    # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
    out.fitqtl.mqm<-fitqtl(fg.cr.obj, pheno.col=i, qtl=revqtl.mqm, formula=my.formula.mqm, method=qtl_method, get.ests=T)
    assign(paste('mqm.out.fitqtl', pname, sep="."), out.fitqtl.mqm)

    # Proportioning of variance for scanone results (full model)
    full.mdl.var.mqm<-as.data.frame(out.fitqtl.mqm$result.full[,2])
    full.mdl.var.mqm$prop.variance<-c(100*(full.mdl.var.mqm[1,1]/full.mdl.var.mqm[3,1]), 100*(full.mdl.var.mqm[2,1]/full.mdl.var.mqm[3,1]), 100)
    colnames(full.mdl.var.mqm)<-c("variance", "prop.variance")
    assign(paste('mqm.full.mdl.var', pname, sep="."), full.mdl.var.mqm)
    write.csv(full.mdl.var.mqm, file=paste(dirmqm,"/prop.var.full.mdl.", pname, ".csv", sep=""), quote=F)

    if (length(m.names) > 1) {
  
        # Proportioning of variance for scanone results (individual QTL as proportion)
        dropone.mdl.var.mqm<-as.data.frame(out.fitqtl.mqm$result.drop[,2])
        dropone.mdl.var.mqm<-rbind(dropone.mdl.var.mqm, full.mdl.var.mqm[3,1])
        # dropone.mdl.var.mqm$prop.variance<-c(100*(dropone.mdl.var.mqm[1,1]/full.mdl.var.mqm[3,1]), 100*(dropone.mdl.var.mqm[2,1]/full.mdl.var.mqm[3,1]), 100)
  
  
        # Partition variance within QTL
        dropone.mdl.var.mqm$prop.variance<-c(100*(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-1),1]/full.mdl.var.mqm[3,1]), 100)
        colnames(dropone.mdl.var.mqm)<-c("variance", "prop.variance")
        rownames(dropone.mdl.var.mqm)<-c(m.names, "total")
        dropone.residual<-c(dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),1]-sum(dropone.mdl.var.mqm[(1:nrow(dropone.mdl.var.mqm)-1),1]), dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),2]-sum(dropone.mdl.var.mqm[(1:nrow(dropone.mdl.var.mqm)-1),2]))
        dropone.mdl.var.mqm<-rbind(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-1),], dropone.residual, dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),])
        rownames(dropone.mdl.var.mqm)[(nrow(dropone.mdl.var.mqm)-1)]<-c('residual')
        assign(paste('mqm.dropone.mdl.var', pname, sep="."), dropone.mdl.var.mqm)
        write.csv(dropone.mdl.var.mqm, file=paste(dirmqm,"/prop.var.dropone.mdl.", pname, ".csv", sep=""), quote=F)
  
        # Get confidence interval (95%) for each QTL:
        n.QTL.mqm<-revqtl.mqm$n.qtl
        CI_L<-c()
        CI_R<-c()
        for(q in 1:n.QTL.mqm){
             lod_int.mqm<-lodint(revqtl.mqm, qtl.index=q)
             lod_int.mqm$marker<-rownames(lod_int.mqm)
             L<-as.data.frame(lod_int.mqm[1,c(4,1,2,3)])
             colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
             R<-as.data.frame(lod_int.mqm[3,c(4,1,2,3)])
             colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
             CI_L<-rbind(CI_L, L)
             CI_R<-rbind(CI_R, R)
          }
          CI<-cbind(CI_L, CI_R)
          assign(paste('mqm.confidence_interval', pname, sep="."), CI)
          write.csv(CI, file=paste(dirmqm,"/mqm.confidence_interval.", pname, ".csv", sep=""), quote=F)
  
  
  
         # Lets get addative effect sizes from the fitqtl model of scanone qtl
         fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mqm)[3])
         fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
         fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
         # Get proportion of variance per marker
         marker.prop.var.mqm<-as.data.frame(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-2),2])
  
         # Build a table that summarizes the results
         summary.table.mqm<-cbind(mqm.markers[,1:ncol(mqm.markers)], marker.prop.var.mqm, fx_size.a, fx_size.se, CI)
         colnames(summary.table.mqm)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
         assign(paste('summary.table.mqm', pname, sep="."), summary.table.mqm)
         write.csv(summary.table.mqm, file=paste(dirmqm,"/summary.table.mqm.", pname, ".csv", sep=""), quote=F)
    }

    if (length(m.names) == 1) {
        # Full model incorporates the only QTL dropone analysis is empty
        write("No paritioning of variance among QTL if only 1 detected.", file=paste(dirmqm,"/prop.var.dropone.mdl.", pname, ".txt", sep=""))
  
        # Get confidence interval (95%) for the QTL:
        lod_int.mqm<-lodint(revqtl.mqm)
        lod_int.mqm$marker<-rownames(lod_int.mqm)
        L<-as.data.frame(lod_int.mqm[1,c(4,1,2,3)])
        colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
        R<-as.data.frame(lod_int.mqm[3,c(4,1,2,3)])
        colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
        CI<-cbind(L, R)
        assign(paste('mqm.confidence_interval', pname, sep="."), CI)
        write.csv(CI, file=paste(dirmqm,"/mqm.confidence_interval.", pname, ".csv", sep=""), quote=F)
  
        # Lets get addative effect sizes from the fitqtl model of scanone qtl
        fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mqm)[2])
        fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
        fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
  
        # Build a table that summarizes the results
        summary.table.mqm<-cbind(mqm.markers[,1:ncol(mqm.markers)], as.data.frame(full.mdl.var.mqm[1,2]), fx_size.a, fx_size.se, CI)
        colnames(summary.table.mqm)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
        assign(paste('summary.table.mqm', pname, sep="."), summary.table.mqm)
        write.csv(summary.table.mqm, file=paste(dirmqm,"/summary.table.mqm.", pname, ".csv", sep=""), quote=F)      
    }
  }

    if(length(stepout.a) == 0) {
      write("No significant QTL in MQM", file=paste(dirmqm,"/empty_summary.table.mqm.", pname, ".txt", sep=""))
      save.image(paste(dirpname, session_image_name, sep="/"))
      if(length(stepout.a.lite) > 0) {
        pdf(file=paste(dirmqm,"/mqm.LOD.profile.", pname, ".pdf", sep=""))
        plotLodProfile(stepout.a.lite, showallchr=T, col="grey", main=pname)
        abline(h=max.perm.so, col="red")
        dev.off()  
        save.image(paste(dirpname, session_image_name, sep="/"))
    }
save.image(file=paste(dirpname,  session_image_name, sep="/" ))
    
}

}

if (nrow(qtl.table.so) == 0 ) {
    write("No significant QTL", file=paste(dirscanone,"/scanone_phe_vals_and_allele_counts_", pname, ".txt", sep=""))
    # Perform co-factor search with a null QTL model
    
    # Perform positive stepwise forward and reverse model selection given formula, fine tuned QTL location, and penalites
    # This is an additive model, run with two different penalties that are specified above
    stepout.a<-stepwiseqtl(fg.cr.obj, pheno.col=i, penalties=pen_heavy, max.qtl=25, scan.pairs=T, additive.only=T, model=trait_model)
    assign(paste('stepout.a', pname, sep="."), stepout.a)
    stepout.a.lite<-stepwiseqtl(fg.cr.obj, pheno.col=i, penalties=pen_lite, max.qtl=25, scan.pairs=T, additive.only=T, model=trait_model)
    assign(paste('stepout.a.lite', pname, sep="."), stepout.a.lite)
    #save.image(file=paste(dirtrait, session_image_name, sep="/" ))
    save.image(file=paste(dirpname,  session_image_name, sep="/" ))
    
    ####################           
    # Stepwise MQM analysis
    ####################
    
    # Perform positive stepwise forward and reverse model selection given formula, fine tuned QTL location, and penalites
    if(length(stepout.a) > 0) {
      pdf(file=paste(dirmqm,"/mqm.LOD.profile.", pname, ".pdf", sep=""))
      plotLodProfile(stepout.a.lite, showallchr=T, col="grey", main=pname)
      plotLodProfile(stepout.a, add=T, showallchr=T, col = "black")
      abline(h=max.perm.so, col="red")
      dev.off()  
      
      assign(paste('stepout.a', pname, sep="."), stepout.a)
      mqm.markers<-c()
      mqm.markers$chr<-stepout.a$chr
      mqm.markers$pos<-stepout.a$pos
      mqm.markers<-as.data.frame(mqm.markers)
      
      mqm.lod.profile<-attr(stepout.a, "lodprofile")
      
      m.names<-c()
      for (m in 1:nrow(mqm.markers)) {
        m.name<-find.marker(fg.cr.obj, chr=mqm.markers[m,1], pos=mqm.markers[m,2])
        m.names<-c(m.names, m.name)
      }
      
      # Get lod scores for each significant marker
      mqm.lod.profile<-attr(stepout.a, "lodprofile")
      #mqm.lod.df<-c()
      lods<-c()
      for(l in 1:length(mqm.lod.profile)) {
          d<-as.data.frame(mqm.lod.profile[l])
          colnames(d)<-c('chr', 'pos', 'lod')
          marker_name<-m.names[l]
          if(length(d[marker_name %in% rownames(d),3]) > 0)  {
              lod<-d[rownames(d) == marker_name,3]
              lods<-c(lods, lod)
          }
          if(length(d[marker_name %in% rownames(d),3]) == 0) {
              lod<-max(d[,3])
              lods<-c(lods, lod)
          }
      }
      
      # Add in the names of the markers and LOD scores
      rownames(mqm.markers)<-m.names
      
      mqm.markers$lod<-lods
      
      qtl.mqm<-makeqtl(fg.cr.obj, mqm.markers$chr, mqm.markers$pos, what=c("prob"))
      
      # Make a formula that specifies this model
      Qs.mqm<-paste('Q', 1:length(m.names), sep="")
      my.formula.mqm<-as.formula(paste("y~", paste(Qs.mqm, collapse="+")))
      assign(paste('mqm.formula', pname, sep="."), my.formula.mqm)
      
      # Fine-tune QTL location
      revqtl.mqm<-refineqtl(fg.cr.obj, pheno.col=i, qtl=qtl.mqm, formula = my.formula.mqm, method=qtl_method)
      assign(paste('mqm.revqtl', pname, sep="."), revqtl.mqm)
      
      # Fit a mulitple QTL model based upon the formula derived above and the fine tuned QTL location
      out.fitqtl.mqm<-fitqtl(fg.cr.obj, pheno.col=i, qtl=revqtl.mqm, formula=my.formula.mqm, method=qtl_method, get.ests=T)
      assign(paste('mqm.out.fitqtl', pname, sep="."), out.fitqtl.mqm)
      
      # Proportioning of variance for scanone results (full model)
      full.mdl.var.mqm<-as.data.frame(out.fitqtl.mqm$result.full[,2])
      full.mdl.var.mqm$prop.variance<-c(100*(full.mdl.var.mqm[1,1]/full.mdl.var.mqm[3,1]), 100*(full.mdl.var.mqm[2,1]/full.mdl.var.mqm[3,1]), 100)
      colnames(full.mdl.var.mqm)<-c("variance", "prop.variance")
      assign(paste('mqm.full.mdl.var', pname, sep="."), full.mdl.var.mqm)
      write.csv(full.mdl.var.mqm, file=paste(dirmqm,"/prop.var.full.mdl.", pname, ".csv", sep=""), quote=F)
      

      if (length(m.names) > 1) {
        
        # Proportioning of variance for scanone results (individual QTL as proportion)
        dropone.mdl.var.mqm<-as.data.frame(out.fitqtl.mqm$result.drop[,2])
        dropone.mdl.var.mqm<-rbind(dropone.mdl.var.mqm, full.mdl.var.mqm[3,1])
        # dropone.mdl.var.mqm$prop.variance<-c(100*(dropone.mdl.var.mqm[1,1]/full.mdl.var.mqm[3,1]), 100*(dropone.mdl.var.mqm[2,1]/full.mdl.var.mqm[3,1]), 100)
        
        
        # Partition variance within QTL
        dropone.mdl.var.mqm$prop.variance<-c(100*(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-1),1]/full.mdl.var.mqm[3,1]), 100)
        colnames(dropone.mdl.var.mqm)<-c("variance", "prop.variance")
        rownames(dropone.mdl.var.mqm)<-c(m.names, "total")
        dropone.residual<-c(dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),1]-sum(dropone.mdl.var.mqm[(1:nrow(dropone.mdl.var.mqm)-1),1]), dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),2]-sum(dropone.mdl.var.mqm[(1:nrow(dropone.mdl.var.mqm)-1),2]))
        dropone.mdl.var.mqm<-rbind(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-1),], dropone.residual, dropone.mdl.var.mqm[nrow(dropone.mdl.var.mqm),])
        rownames(dropone.mdl.var.mqm)[(nrow(dropone.mdl.var.mqm)-1)]<-c('residual')
        assign(paste('mqm.dropone.mdl.var', pname, sep="."), dropone.mdl.var.mqm)
        write.csv(dropone.mdl.var.mqm, file=paste(dirmqm,"/prop.var.dropone.mdl.", pname, ".csv", sep=""), quote=F)
        
        # Get confidence interval (95%) for each QTL:
        n.QTL.mqm<-revqtl.mqm$n.qtl
        CI_L<-c()
        CI_R<-c()
        for(q in 1:n.QTL.mqm){
          lod_int.mqm<-lodint(revqtl.mqm, qtl.index=q)
          lod_int.mqm$marker<-rownames(lod_int.mqm)
          L<-as.data.frame(lod_int.mqm[1,c(4,1,2,3)])
          colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
          R<-as.data.frame(lod_int.mqm[3,c(4,1,2,3)])
          colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
          CI_L<-rbind(CI_L, L)
          CI_R<-rbind(CI_R, R)
        }
        CI<-cbind(CI_L, CI_R)
        assign(paste('mqm.confidence_interval', pname, sep="."), CI)
        write.csv(CI, file=paste(dirmqm,"/mqm.confidence_interval.", pname, ".csv", sep=""), quote=F)
        
        
        
        # Lets get addative effect sizes from the fitqtl model of scanone qtl
        fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mqm)[3])
        fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
        fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
        
        # Get proportion of variance per marker
        marker.prop.var.mqm<-as.data.frame(dropone.mdl.var.mqm[1:(nrow(dropone.mdl.var.mqm)-2),2])
        
        # Build a table that summarizes the results
        summary.table.mqm<-cbind(mqm.markers[,1:ncol(mqm.markers)], marker.prop.var.mqm, fx_size.a, fx_size.se, CI)
        colnames(summary.table.mqm)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
        assign(paste('summary.table.mqm', pname, sep="."), summary.table.mqm)
        write.csv(summary.table.mqm, file=paste(dirmqm,"/summary.table.mqm.", pname, ".csv", sep=""), quote=F)
      }
      
      if (length(m.names) == 1) {
        # Full model incorporates the only QTL dropone analysis is empty
        write("No paritioning of variance among QTL if only 1 detected.", file=paste(dirmqm,"/prop.var.dropone.mdl.", pname, ".txt", sep=""))
        
        # Get confidence interval (95%) for the QTL:
        lod_int.mqm<-lodint(revqtl.mqm)
        lod_int.mqm$marker<-rownames(lod_int.mqm)
        L<-as.data.frame(lod_int.mqm[1,c(4,1,2,3)])
        colnames(L)<-c("L.CI_marker", "L.CI_chr", "L.CI_pos", "L.CI_lod")
        R<-as.data.frame(lod_int.mqm[3,c(4,1,2,3)])
        colnames(R)<-c("R.CI_marker", "R.CI_chr", "R.CI_pos", "R.CI_lod")
        CI<-cbind(L, R)
        assign(paste('mqm.confidence_interval', pname, sep="."), CI)
        write.csv(CI, file=paste(dirmqm,"/mqm.confidence_interval.", pname, ".csv", sep=""), quote=F)
        
        # Lets get addative effect sizes from the fitqtl model of scanone qtl
        fx_size_fitqtl<-as.data.frame(summary(out.fitqtl.mqm)[2])
        fx_size.a<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),1])
        fx_size.se<-as.data.frame(fx_size_fitqtl[2:nrow(fx_size_fitqtl),2])
        
        # Build a table that summarizes the results
        summary.table.mqm<-cbind(mqm.markers[,1:ncol(mqm.markers)], as.data.frame(full.mdl.var.mqm[1,2]), fx_size.a, fx_size.se, CI)
        colnames(summary.table.mqm)<-c('chr', 'pos', 'lod', 'prop.var', 'additive.fx', 'additive.fx_se', colnames(CI))
        assign(paste('summary.table.mqm', pname, sep="."), summary.table.mqm)
        write.csv(summary.table.mqm, file=paste(dirmqm,"/summary.table.mqm.", pname, ".csv", sep=""), quote=F)      
      }
    }
    
    if(length(stepout.a) == 0) {
      write("No significant QTL in MQM", file=paste(dirmqm,"/empty_summary.table.mqm.", pname, ".txt", sep=""))
      
      if(length(stepout.a.lite) > 0) {
        pdf(file=paste(dirmqm,"/mqm.LOD.profile.", pname, ".pdf", sep=""))
        plotLodProfile(stepout.a.lite, showallchr=T, col="grey", main=pname)
        abline(h=max.perm.so, col="red")
        dev.off()  
     save.image(paste(dirpname, session_image_name, sep="/"))
    }
    
  save.image(file=paste(dirpname,  session_image_name, sep="/" ))

  }
}

rm(limit.so,max.out.so, max.lod.so, max.perm.so, operm.so, out.so, pname, i, dirscanone, dir, args)
save.image(paste(dirpname, session_image_name, sep="/"))
rm(list=ls())
