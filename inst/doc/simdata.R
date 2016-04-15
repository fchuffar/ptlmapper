## ------------------------------------------------------------------------
library(ptldata)
indiv = ptldata::indiv
genodata = ptldata::genodata
cells = ptldata::cells

## ------------------------------------------------------------------------
err = 0.2
nb_indiv = 80
idx = which(indiv$err==err)[1:nb_indiv]
cells = cells[idx]
indiv = indiv[idx,]
genodata = genodata[,idx]
head(indiv)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1), respect=TRUE)
for (m in c("m1","m2","m3")){  
  h = hist(indiv[[m]], nclass=20, main=m, xlab="", ylab="", col="grey")
  for (all in unique(indiv$all)) {
    h = hist(indiv[indiv$all==all,][[m]], plot=FALSE, breaks=h$breaks)    
    lines(h, col=adjustcolor(all*2, alpha=0.3))
  }
}

## ------------------------------------------------------------------------
genodata[96:106, 1:7]

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
plot(0,0, col=0, xlim=c(-1.5,27), ylim=c(0,0.3), 
  main="distribution of phenotypes", xlab="single cell value", ylab="density")
foo = sapply(rownames(indiv), function(i) {
  lines(density(cells[[i]], bw=1), col=indiv[i,"all"]*2)
})

## ------------------------------------------------------------------------
library(ptlmapper)

## ------------------------------------------------------------------------
pheno_hists = build_pheno_hists(cells, bin_width=1)

## ------------------------------------------------------------------------
kd_matrix = build_kd_matrix(pheno_hists)
kd_matrix[1:6,1:6]

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
kanto_analysis = ptl_scan(kd_matrix, genodata_ptl, method="kanto")
mmoments_analysis = ptl_scan(mm_matrix, genodata_ptl, method="mmoments")

## ------------------------------------------------------------------------
rqtl_analysis = rqtl_launch(genodata_ptl, pheno_hists, kanto_analysis, mmoments_analysis)

## ---- results='hide'-----------------------------------------------------
ptl_mapping_result = ptl_mapping(genodata, cells, bckg, nb_perm=20, bin_width=1)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
plot_rqtl(ptl_mapping_result, which_pheno=1)
plot_rqtl(ptl_mapping_result, which_pheno=2)
plot_rqtl(ptl_mapping_result, which_pheno=3)

## ----fig.height=3, fig.width=9-------------------------------------------
best_marker_kanto = get_best_markers_rptl(ptl_mapping_result, method="kanto")
print(best_marker_kanto)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
best_marker_name = rownames(best_marker_kanto)[1]
col = marker2col(ptl_mapping_result, best_marker_name)
plot_orth_trans(ptl_mapping_result, col=col, method="kanto")
plot_wilks(ptl_mapping_result, method="kanto")
plot_can(ptl_mapping_result, best_marker_name, col=col, method="kanto")

## ----fig.height=3, fig.width=9-------------------------------------------
best_marker_mmoments = get_best_markers_rptl(ptl_mapping_result, method="mmoments")
print(best_marker_mmoments)

## ----fig.height=3, fig.width=9-------------------------------------------
layout(matrix(1:3, 1, byrow=TRUE), respect=TRUE)
best_marker_name = rownames(best_marker_mmoments)[1]
col = marker2col(ptl_mapping_result, best_marker_name)
plot_orth_trans(ptl_mapping_result, col=col, method="mmoments")
plot_wilks(ptl_mapping_result, method="mmoments")
plot_can(ptl_mapping_result, best_marker_name, col=col, method="mmoments")

