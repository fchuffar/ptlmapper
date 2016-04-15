## ------------------------------------------------------------------------
library(ptldata)
indiv = ptldata::indiv
genodata = ptldata::genodata
cells = ptldata::cells

## ------------------------------------------------------------------------
head(indiv)

## ------------------------------------------------------------------------
err = 0.2
nb_indiv = 80
idx = which(indiv$err==err)[1:nb_indiv]
cells = cells[idx]
indiv = indiv[idx,]
genodata = genodata[,idx]

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1), respect=TRUE)
for (m in c("m1","m2","m3")){  
  h = hist(indiv[[m]], nclass=20, 
    main=paste("moment", m, sep=""), xlab="", ylab="", 
    col=adjustcolor("grey", alpha=0.3))
  for (all in unique(indiv$all)) {
    h = hist(indiv[indiv$all==all,][[m]], plot=FALSE, breaks=h$breaks)    
    lines(h, col=adjustcolor(all+1, alpha=0.3))
  }
}

## ------------------------------------------------------------------------
genodata[1:10, 1:7]

## ----fig.height=6, fig.width=6-------------------------------------------
tetas = compute_teta(genodata)
matlab::imagesc(t(as.matrix(genodata)), 
  main="Parental origin of segregant genomes", 
  xlab="markers", ylab="segrgants")

## ----fig.height=6, fig.width=6-------------------------------------------
plot(density(tetas), 
  main=paste("Recombination Fraction Distribution"), 
  xlab=paste("mean(teta)=", signif(mean(tetas),3), sep=""))
abline(v=mean(tetas), lty=2)

## ------------------------------------------------------------------------
head(names(cells))

## ----fig.width=6, fig.height=6-------------------------------------------
# plot_dist(ptl_mapping_result, main=main, col=col)
library(ptlmapper)

## ------------------------------------------------------------------------
pheno_hists = build_pheno_hists(cells, bin_width=1)

## ------------------------------------------------------------------------
kd_matrix = build_kd_matrix(pheno_hists)
head(kd_matrix)

## ------------------------------------------------------------------------
mm_matrix = build_mmoments_matrix(pheno_hists, nb_moment=4)
head(mm_matrix)

## ------------------------------------------------------------------------
bckg = names(cells)
genodata$chromosome = 1
genodata$position = 1:nrow(genodata)
genodata$prob_name = rownames(genodata)
genodata$rec_fractions = 0.05
genodata_ptl = preprocess_genodata(genodata, bckg)

## ------------------------------------------------------------------------
kanto_analysis = ptl_scan(kd_matrix, genodata_ptl, nb_perm=3, method="kanto")
mmoments_analysis = ptl_scan(mm_matrix, genodata_ptl, nb_perm=3, method="mmoments")

## ------------------------------------------------------------------------
rqtl_analysis = rqtl_launch(genodata_ptl, pheno_hists, kanto_analysis, mmoments_analysis)

## ------------------------------------------------------------------------
ptl_mapping_result = ptl_mapping(genodata_ptl, cells, bckg, nb_perm=3, bin_width=1)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
plot_rqtl(ptl_mapping_result, which_pheno=1)
plot_rqtl(ptl_mapping_result, which_pheno=2)
plot_rqtl(ptl_mapping_result, which_pheno=3)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
plot_mds(ptl_mapping_result, main=main, col=col)
plot_wilks(ptl_mapping_result, main=main, method="kanto")
plot_can(ptl_mapping_result, marker_names[1], main=main, col=col)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
plot_wilks(ptl_mapping_result=ptl_mapping_result, pa=ptl_mapping_result$mmoments_analysis, main=paste(nb_indiv, "indiv (mm)"), method="mmoments")
plot_dist(ptl_mapping_result, main=main, col=col)

