## ------------------------------------------------------------------------
library(ptldata)

## ------------------------------------------------------------------------
err = 0.2
nb_indiv = 80
idx = sample(which(indiv$err==err), nb_indiv)
cells = cells[idx]
indiv = indiv[idx, ]
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

## ----fig.width=6, fig.height=6-------------------------------------------
# plot_dist(ptl_mapping_result, main=main, col=col)
library(ptlmapper)

## ------------------------------------------------------------------------
pheno_hists = build_pheno_hists(cells, bin_width=1)

## ------------------------------------------------------------------------
kd_matrix = build_kd_matrix(pheno_hists)

## ------------------------------------------------------------------------
mm_matrix = build_mmoments_matrix(pheno_hists)

## ------------------------------------------------------------------------
bckg = names(genodata)
genodata$chromosome = rep(1, nrow(genodata))
genodata$position = 1:nrow(genodata)
genodata$prob_name = paste("marker", 1:nrow(genodata), sep="_")
genodata$rec_fractions = 0.05

genodata_ptl = preprocess_genodata(genodata, bckg)

## ------------------------------------------------------------------------
kanto_analysis = ptl_scan(kd_matrix, genodata_ptl, nb_perm=3)

## ------------------------------------------------------------------------
mmoments_analysis = ptl_scan(mm_matrix, genodata_ptl, nb_perm=3, method="mmoment")

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------

## ----fig.height=6, fig.width=12------------------------------------------


###################################################
### code chunk number 7: simdata_pres.Rnw:553-556
###################################################

main=paste("err=", err, sep="")
# ptl_mapping_filename = paste("cache/simdata", err, nb_perm, "ptl_mapping.Rds", sep="_")
# genodata = g
# cells = c
# bckg = bckg
# nb_perm=nb_perm
# nb_dim=0
# bin_width=1
# ptl_mapping_filename=ptl_mapping_filename
# DO_MMOMENTS=TRUE
# errs = c(0.05, 0.01, 0.005)
# minimal_proportion=0.1
# COMPUTE_KD_MATRIX=TRUE
# DO_KANTO=TRUE
# DO_RQTL=TRUE
# nb_moments=3
# LIGHTWEIGHT=FALSE
ptl_mapping_result = ptl_mapping(genodata, cells, bckg, nb_perm=3, nb_dim=0, bin_width=1, 
# ptl_mapping_filename=ptl_mapping_filename, 
DO_MMOMENTS=TRUE)

layout(matrix(1:8, 2, byrow=TRUE), respect=TRUE)
plot_rqtl(ptl_mapping_result, main=paste(main, dim(genodata)[1], "markers and", dim(genodata)[2], "indiviuals"), which_pheno=1)
plot_rqtl(ptl_mapping_result, main=main, which_pheno=2)
marker_names = rownames(get_best_markers_rptl(ptl_mapping_result))
col = mn_2_col(ptl_mapping_result, marker_names[1])
plot_mds(ptl_mapping_result, main=main, col=col)
plot_wilks(ptl_mapping_result, main=main, method="kanto")
plot_wilks(ptl_mapping_result=ptl_mapping_result, pa=ptl_mapping_result$mmoments_analysis, main=paste(nb_indiv, "indiv (mm)"), method="mmoment")
plot_can(ptl_mapping_result, marker_names[1], main=main, col=col)
# plot_noise(ptl_mapping_result, main=main, col=col)
plot_dist(ptl_mapping_result, main=main, col=col)
plot_empirical_test(ptl_mapping_result,main=main)

