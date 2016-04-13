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

## ----fig.height=6, fig.width=12------------------------------------------


###################################################
### code chunk number 7: simdata_pres.Rnw:553-556
###################################################
nb_perm = 0
layout(matrix(1:8, 2, byrow=TRUE), respect=TRUE)
g = genodata
bckg = names(g)
nb_allele = 1
g$chromosome = as.vector(sapply(1:nb_allele, function(i){rep(i, nrow(g))}))
g$RQTL.name = paste("marker", 1:nrow(g), sep="_")
g$rec.fractions = 0.01
g$position = 1:nrow(g)
# "position"
c = cells
main=paste("err=", err, sep="")
ptl_mapping_filename = paste("cache/simdata", err, nb_perm, "ptl_mapping.Rds", sep="_")
# genodata = g
# cells = c
# bckg = bckg
# nb_perm=nb_perm
# types="raw"
# nb_dim=0
# bin_width=1
# ptl_mapping_filename=ptl_mapping_filename
# USE_BOT_FOR_EXPLORATION=FALSE
# DO_MMOMENTS=TRUE
# errs = c(0.05, 0.01, 0.005)
# minimal_proportion=0.1
# COMPUTE_KD_MATRIX=TRUE
# DO_KANTO=TRUE
# DO_RQTL=TRUE
# nb_moments=3
# LIGHTWEIGHT=FALSE
ptl_mapping_result = ptl_mapping(g, c, bckg, nb_perm=nb_perm, types="raw", nb_dim=0, bin_width=1, ptl_mapping_filename=ptl_mapping_filename, USE_BOT_FOR_EXPLORATION=FALSE, DO_MMOMENTS=TRUE)

plot_rqtl(ptl_mapping_result, main=paste(main, dim(g)[1], "markers and", dim(g)[2], "indiviuals"), which_pheno=1)
plot_rqtl(ptl_mapping_result, main=main, which_pheno=2)
marker_names = rownames(get_best_markers_rptl(ptl_mapping_result))
col = mn_2_col(ptl_mapping_result, marker_names[1])
plot_mds(ptl_mapping_result, main=main, col=col)
plot_wilks(ptl_mapping_result,main=main)
plot_can(ptl_mapping_result, marker_names[1], main=main, col=col)
plot_noise(ptl_mapping_result, main=main, col=col)
plot_dist(ptl_mapping_result, main=main, col=col)
plot_empirical_test(ptl_mapping_result,main=main)



# ###################################################
# ### code chunk number 18: simdata_pres.Rnw:653-674
# ###################################################
# nb_perm = 20
# layout(matrix(1:44, 4, byrow=FALSE), respect=TRUE)
# for (err in c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20)) {
#   idx = which(indiv$err==err)
#   nb_indiv = length(idx)
#   g = genodata[,idx]
#   bckg = names(g)
#   nb_allele = 1
#   g$chromosome = as.vector(sapply(1:nb_allele, function(i){rep(i, nrow(g))}))
#   g$RQTL.name = paste("marker", 1:nrow(g), sep="_")
#   g$rec.fractions = 0.01
#   # "position"
#   c = cells[idx]
#   main=paste("err=", err, sep="")
#   ptl_mapping_filename = paste("cache/simdata", err, nb_perm, "ptl_mapping.rds", sep="_")
#   ptl_mapping_result = ptl_mapping(g, c, bckg, nb_perm=nb_perm, types="raw", nb_dim=0, bin_width=1, ptl_mapping_filename=ptl_mapping_filename)
#   plot_rqtl(ptl_mapping_result, main=main, which_pheno=1, ylim=c(0, 6))
#   plot_rqtl(ptl_mapping_result, main=main, which_pheno=2, ylim=c(0, 10))
#   plot_wilks(ptl_mapping_result,main=main, ylim=c(0, 40))
#   plot_wilks(ptl_mapping_result=NULL, pa=ptl_mapping_result$mmoments_analysis[["raw"]],main=paste(main, nb_indiv, "(mm)"), ylim=c(0, 40))
# }
#
#
# ###################################################
# ### code chunk number 19: simdata_pres.Rnw:683-709
# ###################################################
# nb_perm = 20
# layout(matrix(1:66, 6, byrow=FALSE), respect=TRUE)
# for (err in c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20)) {
#   for (nb_indiv in c(125, 100, 75)) {
#     set.seed(seed=0)
#     # over_idx = runif(nrow(indiv)) <= prop_indiv
#     # idx = which((indiv$err==err)[over_idx])
#     # nb_indiv = length(idx)
#     idx = sample(which((indiv$err==err)), nb_indiv)
#     g = genodata[,idx]
#     bckg = names(g)
#     nb_allele = 1
#     g$chromosome = as.vector(sapply(1:nb_allele, function(i){rep(i, nrow(g))}))
#     g$RQTL.name = paste("marker", 1:nrow(g), sep="_")
#     g$rec.fractions = 0.01
#     # "position"
#     c = cells[idx]
#     main=paste("err=", err, sep="")
#     ptl_mapping_filename = paste("cache/simdata", err, nb_perm, nb_indiv, "_less_indiv_ptl_mapping.rds", sep="_")
#     ptl_mapping_result = ptl_mapping(g, c, bckg, nb_perm=nb_perm, types="raw", nb_dim=0, bin_width=1, ptl_mapping_filename=ptl_mapping_filename, DO_MMOMENTS=TRUE)
#     # plot_rqtl(ptl_mapping_result, main=main, which_pheno=1, ylim=c(0, 6))
#     # plot_rqtl(ptl_mapping_result, main=main, which_pheno=2, ylim=c(0, 10))
#     plot_wilks(ptl_mapping_result, ylim=c(0, 10), main=paste(main, "(kanto)"))
#     plot_wilks(ptl_mapping_result=NULL, pa=ptl_mapping_result$mmoments_analysis[["raw"]],main=paste(nb_indiv, "indiv (mm)"), ylim=c(0, 10))
#   }
# }
#
#
# ###################################################
# ### code chunk number 20: simdata_pres.Rnw:716-742
# ###################################################
# nb_perm = 20
# layout(matrix(1:66, 6, byrow=FALSE), respect=TRUE)
# for (err in c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20)) {
#   for (nb_indiv in c(50, 25, 15)) {
#     set.seed(seed=0)
#     # over_idx = runif(nrow(indiv)) <= prop_indiv
#     # idx = which((indiv$err==err)[over_idx])
#     # nb_indiv = length(idx)
#     idx = sample(which((indiv$err==err)), nb_indiv)
#     g = genodata[,idx]
#     bckg = names(g)
#     nb_allele = 1
#     g$chromosome = as.vector(sapply(1:nb_allele, function(i){rep(i, nrow(g))}))
#     g$RQTL.name = paste("marker", 1:nrow(g), sep="_")
#     g$rec.fractions = 0.01
#     # "position"
#     c = cells[idx]
#     main=paste("err=", err, sep="")
#     ptl_mapping_filename = paste("cache/simdata", err, nb_perm, nb_indiv, "_less_indiv_ptl_mapping.rds", sep="_")
#     ptl_mapping_result = ptl_mapping(g, c, bckg, nb_perm=nb_perm, types="raw", nb_dim=0, bin_width=1, ptl_mapping_filename=ptl_mapping_filename, DO_MMOMENTS=TRUE,)
#     # plot_rqtl(ptl_mapping_result, main=main, which_pheno=1, ylim=c(0, 6))
#     # plot_rqtl(ptl_mapping_result, main=main, which_pheno=2, ylim=c(0, 10))
#     plot_wilks(ptl_mapping_result, ylim=c(0, 10), main=paste(main, "(kanto)"))
#     plot_wilks(ptl_mapping_result=NULL, pa=ptl_mapping_result$mmoments_analysis[["raw"]],main=paste(nb_indiv, "indiv (mm)"), ylim=c(0, 10))
#   }
# }
#
# ```
#

