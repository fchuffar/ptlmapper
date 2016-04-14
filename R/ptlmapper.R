#' A Function That Extracts Axis Information from a 'ptl_mapping' Data Structure.
#'
#' This function extracts axis information from a 'ptl_mapping' data structure. 
#' @param ptl_mapping_result A `ptl_mapping` data structure.
#' @param delta A numeric proportional to the space lets between chromosomes on the x axis.
extract_axis_info = function(ptl_mapping_result, delta = 0.6) {
  sub_genodata = genodata[ptl_mapping_result$genodata$prob_name %in% names(ptl_mapping_result$genodata_ptl),]
  chrs = sub_genodata$chromosome        
  xs = unlist(
  lapply(1:length(unique(chrs)), function(i) {
    chr = unique(chrs)[i]
    which(sub_genodata$chromosome == chr)
    nb_ticks = sum(sub_genodata$chromosome == chr)
    if (nb_ticks>1) {
      start = i - (delta/2)
      ret = start + ((((1:nb_ticks) - 1) / (nb_ticks-1)) * delta)
    } else {
      ret = i
    }
    return(ret)
  }))
  return(list(chrs=chrs, xs=xs))
}

#' A Function That Scans Genome to Detect QTLs (based on R package `qtl`).
#'
#' This function embed R package `qtl`. It scans genome to detect QTLs using `ptlmapper` data structires. 
#' @param rqtl_data Preprocessed genotypes (typically outpout of preprocess_data_rqtl function).
#' @param nb_perm An integer that specifies the number of permuation to do.
#' @param errs A vector of integer (error) that will be used to compute threshold from the permutation test.
#' @export
rqtl_launch = function(rqtl_data, nb_perm, errs=0.05) {
  rqtl_data$scan = lapply(1:length(rqtl_data$pheno), function(i) {
    tmp_scan_output = scanone(rqtl_data, method="imp", pheno.col=i)
    nb_perm = max(nb_perm, 1000)
    permtest = scanone(rqtl_data, method="hk", n.perm=nb_perm, pheno.col=i)
    thres = summary(permtest, alpha=errs)
    # summary(ptl_mapping_result$rqtl_data$scan[[3]]$scan_output , perms=ptl_mapping_result$rqtl_data$scan[[3]]$permtest, alpha=0.2, pvalues=TRUE)
    return(list(scan_output = tmp_scan_output, thres=thres, errs=errs, permtest=permtest))
  })
  return(rqtl_data)
}


#' A Function That Scans Genome to Detect PTLs.
#'
#' This function scans genome to detect PTLs. 
#' It determines with the Kaiser criterion how many dimensions of the MDS are relevants for canonical analysis (`eig_to_scan`).
#' If nb_dim=0, a z-core criterion id used to determine with how many dimensions the canonical analysis gives the best result (to fixe 2 <= nb_dim <= eig_to_scan).
#' If nb_dim=NULL, `nb_dim` is automatically set to eig_to_scan.
#' If nb_dim>0 the canonical analysis is performed used provided `nb_dim`. In this case, keep attention to the fact that a high value of `nb_dim` introduce a lot of degree of freedom in your canonical analysis. The number of dimension `nb_dim` needs to be __substantially__ smaller than the size of your population. 
#' The permutations batch uses `ptl_scan` function with the `nb_dim` value as the one provided by the initial ptl_scan call.
#' @param kd_matrix The Kantorivitch distance matrix.
#' @param genodata_ptl The preprocessed for ptl scanning genome.
#' @param nb_perm An integer that specifies the number of permuation to do.
#' @param nb_dim An integer that specifies the number of dimension of the MDS space to explore.
#' @param SHOW_SCAN_PROG=FALSE A boolean specifying if genome scan progression need to to report on console.
#' @param SHOW_PERM_PROG=TRUE A boolean specifying if permutation progression need to to report on console.
#' @param scan_prog_freq=50 An integer specifying the frequency of the scan progression reporting.
#' @param perm_prog_freq=5 An integer specifying the frequency of the permution progression reporting.
#' @param method A character string in c("kanto", "mmoment") that specify the method to use to characterize phenotypic distribution. Default is our favourite: "kanto"!
#' @export
ptl_scan = function(pheno_matrix, genodata_ptl, nb_perm, nb_dim=0, SHOW_SCAN_PROG=FALSE, SHOW_PERM_PROG=TRUE, scan_prog_freq=50, perm_prog_freq=5, method="kanto") {
  # function used to compute the number of significants eigen values
  get_nb_eig_sign = function(pheno_matrix, method=method) {
    if (method == "kanto") {
      mds = cmdscale(pheno_matrix,2,eig=TRUE)
      mds_eig = mds$eig[mds$eig > 0]
      nb_eig_sign = sum(mds_eig/sum(mds_eig) >= 1/length(mds$eig))      
    } else {
      pca = prcomp(t(pheno_matrix), scale=TRUE)
      pca_var = pca$sdev * pca$sdev
      nb_eig_sign = sum(pca_var/sum(pca_var) >= 1/length(pca_var))
    }
    return(nb_eig_sign)
  }
  #######
  # GO! #
  #######
  nb_dim_orig = nb_dim
  # How many Dimension ?
  NEED_TO_COMPUTE_MDS_AND_CAN = TRUE
  dim_scan = dim_scan_zscores = eig_to_scan = NULL
  if (is.null(nb_dim)) {
    nb_dim = max(2, get_nb_eig_sign(pheno_matrix, method=method))
    print(paste("number of dimension for can (nb_eig_sign):", nb_dim))
  } else if (nb_dim==0) {
    eig_to_scan = 2:max(2, get_nb_eig_sign(pheno_matrix, method=method))
    dim_scan = apply(t(eig_to_scan), 2, function(test_dim){
      print(paste("Testing nb_dim=", test_dim, "/", max(eig_to_scan), "...", sep=""))
      ptl_scan(pheno_matrix, genodata_ptl, nb_perm=0, nb_dim=test_dim, method=method)
    })
    dim_scan_zscores = sapply(dim_scan, function(d){
    return(zscore(-log10(d$Ws)))
    })
    best_zscore_idx = which(dim_scan_zscores == max(dim_scan_zscores))[1]
    nb_dim = (eig_to_scan)[best_zscore_idx]
    print(paste("number of dimension for can (best zscore):", nb_dim))
    # Do not recalculated MDS nor CAN
    mds = dim_scan[[best_zscore_idx]]$mds
    Ws = dim_scan[[best_zscore_idx]]$Ws
    cans = dim_scan[[best_zscore_idx]]$cans
    NEED_TO_COMPUTE_MDS_AND_CAN = FALSE
    # Clean objects ?
    dim_scan = lapply(dim_scan, function(p_scan){
      p_scan$genodata_ptl = NULL
      p_scan$mds = NULL
      p_scan$cans = NULL
      p_scan
    })
  } else {
    print(paste("number of dimension for can (forced):", nb_dim))
  }

  # Wilks score for our genodata_ptl (collect mds WS and cans if needed)
  if (NEED_TO_COMPUTE_MDS_AND_CAN) {      
    if (method == "kanto") {
      mds = cmdscale(pheno_matrix,nb_dim,eig=TRUE)
      data = data.frame(mds$points)      
    } else {
      mds = prcomp(t(pheno_matrix), scale=TRUE)
      data = data.frame(mds$rotation[,1:nb_dim])      
    }
    foo = apply(t(1:length(genodata_ptl)), 2, function(i) {
      if (i %% scan_prog_freq == 1 & SHOW_SCAN_PROG) {print(paste("Marker ", i, "/", length(genodata_ptl), sep=""))}
      all = genodata_ptl[,i]
      get_wilks_score(data, all)
    })
    foo = do.call(rbind,foo)
    Ws = unlist(foo[,"W"])
    cans = foo[,"can"]
  }

  # Perform permutation if needed.
  if (nb_perm > 0) {
    perm_Ws = sapply(1:nb_perm, function(p){
      if (p %% perm_prog_freq == 1 & SHOW_PERM_PROG) {print(paste("permutation", p , "/", nb_perm))}
        perm_index = sample(1:nrow(genodata_ptl))
        perm_geno = genodata_ptl[perm_index,]
        kanto_analysis = ptl_scan(pheno_matrix, perm_geno, nb_perm=0, nb_dim=nb_dim_orig, SHOW_SCAN_PROG=SHOW_SCAN_PROG, scan_prog_freq=scan_prog_freq, method=method)
        ret = min(kanto_analysis$Ws)
        return(ret)            
    })
  } else {
    perm_Ws = NULL
  }
  
  ptl_analysis = list(genodata=genodata_ptl, mds=mds, Ws=Ws, cans=cans, 
      eig_to_scan=eig_to_scan, dim_scan=dim_scan, dim_scan_zscores=dim_scan_zscores, nb_dim=nb_dim, 
      nb_perm=nb_perm, perm_Ws=perm_Ws)
  return(ptl_analysis)
}

#' A Function That Compute the Wilks Score.
#'
#' This function compute the Wilks score. The score formula is extracted from source code of the candisc package.
#' @param data Multivariate coordinates of individuals.
#' @param all A vector of factors that describes groups of individuals.
#' @importFrom candisc candisc
#' @export
get_wilks_score = function(data, all) {
  # On the web:
  #  * http://www.philender.com/courses/multivariate/notes2/can1.html
  #  * http://www.philender.com/courses/multivariate/notes3/manova.html
  # NOTE: The score formula is extracted from source code of the candisc package.
  # seqWilks <- function (eig, p, df.h, df.e) {
  #   p.full <- length(eig)
  #   result <- matrix(0, p.full, 4)
  #   m <- df.e + df.h - ( p.full + df.h + 1)/2
  #   for (i in seq(p.full)) {
  #     test <- prod(1/(1 + eig[i:p.full]))
  #     p <- p.full + 1 - i
  #     q <- df.h + 1 - i
  #     s <- p^2 + q^2 - 5
  #     if (s > 0) {
  #       s = sqrt(((p * q)^2 - 4)/s)
  #     } else {
  #       s = 1
  #     }
  #     df1 <- p * q
  #     df2 <- m * s - (p * q)/2 + 1
  #     result[i,] <- c(test, ((test^(-1/s) - 1) * (df2/df1)),
  #       df1, df2)
  #   }
  #   result <- cbind(result, pf(result[,2], result[,3], result[,4], lower.tail = FALSE))
  #   colnames(result) <- c("LR test stat", "approx F", "num Df", "den Df", "Pr(> F)")
  #   rownames(result) <- 1:p.full
  #   return(result)
  # }
  data$allele = as.factor(all)
  data = data[!is.na(data$allele),]
  model_formula = paste("cbind(", paste(names(data)[-length(names(data))], collapse=", "), ") ~ allele")
  mod = lm(model_formula, data=data)
  can = candisc::candisc(mod, data=data)
  p  = can$rank
  eig = can$eigenvalues[1:p]
  # df.h =  can$dfh
  # df.e = can$dfe
  # tests = seqWilks(eig, p, df.h, df.e)
  # return(list(W=tests[1,5], can=can, tests=tests))
  W = prod(1/(1 + eig[1:length(eig)]))
  return(list(W=W, can=can))
}


#' A Function That Preprocesses Genotype of Individuals.
#'
#' This function preprocesses genotype of individuals. 
#' In geneodata, "2" stands for undetermined and is replaced by NA.
#' @param genodata A matrix describing the genotype of individuals.
#' @param bckg Vector of character describing the individuals as they are describe in genodata.
#' @param minimal_proportion A numeric specifying the minimal proportion of each parental allele under which a marker is discard from the analysis.
#' @export
preprocess_genodata = function(genodata, bckg, minimal_proportion=0.1) {
  g = genodata
  colnames =  as.character(g$prob_name)
  g =  data.frame(cbind(bckg, t(g[, bckg])), stringsAsFactors=FALSE)
  colnames(g) = c("spores", colnames)
  # Check genodata
  g[g=="2"] = NA
  kept = sapply(1:ncol(g), function(i) {
    col = g[bckg, i]
    nb_ind = length(col)
    inds = as.vector(na.omit(col))
    u_inds = unique(inds)
    for (i in u_inds) {
      if (sum(inds==i)/nb_ind < minimal_proportion | sum(inds==i)/nb_ind == 1) {
        # print("Removing something")
        return(FALSE)
      }
    }
    return(TRUE)
  })
  return(g[bckg, kept])
}

#' A Function That Preprocesses Genotypes and Phenotypes of Individuals for R Package QTL.
#'
#' This function preprocesses genotypes and phenotypes of individuals for R Package QTL. 
#' @inheritParams preprocess_phenodata_rqtl
#' @inheritParams preprocess_genodata_rqtl
#' @export
preprocess_data_rqtl = function(genodata, bckg, pheno_hists, kanto_analysis=NULL, mmoments_analysis=NULL) {
  genodata_rqtl = preprocess_genodata_rqtl(genodata, bckg)
  phenodata_rqtl = preprocess_phenodata_rqtl(pheno_hists, kanto_analysis=kanto_analysis, mmoments_analysis=mmoments_analysis)
  rqtl_data = list()
  class(rqtl_data) = c("bc", "cross")
  rqtl_data$geno = list()
  for (chr in unique(genodata_rqtl$chromosome)) {
    chr = as.character(chr)
    rqtl_data$geno[[chr]] = list()
    rqtl_data$geno[[chr]]$data = t(genodata_rqtl[ genodata_rqtl$chromosome == chr, bckg])
    colnames(rqtl_data$geno[[chr]]$data) = genodata_rqtl[ genodata_rqtl$chromosome == chr, ]$prob_name
    rqtl_data$geno[[chr]]$map = as.list(cumsum(genodata_rqtl[genodata_rqtl$chromosome == chr, ]$rec_fractions * 100))
    names(rqtl_data$geno[[chr]]$map) = genodata_rqtl[genodata_rqtl$chromosome == chr, ]$prob_name
    rqtl_data$geno[[chr]]$map = unlist(rqtl_data$geno[[chr]]$map)
    class(rqtl_data$geno[[chr]]) = "A"
  }
  rqtl_data$pheno = data.frame(phenodata_rqtl)
  rqtl_data = jittermap(rqtl_data)
  rqtl_data = calc.genoprob(rqtl_data, step=10, err=0.01)
  rqtl_data = sim.geno(rqtl_data, step=10, n.draws=8, err=0.01)
  return(rqtl_data)
}

#' A Function That Preprocesses Phenotypes of Individuals for R Package QTL.
#'
#' This function preprocesses phenotypes of individuals for R Package QTL. 
#' @param pheno_hists A list of object containing mean and var attribute. Typically outpout of the `build_pheno_hists` function.
#' @param kanto_analysis Output of `ptl_scan` function used with the method "kanto".
#' @param mmoments_analysis Output of `ptl_scan` function used with the method "mmoment". 
#' @export
preprocess_phenodata_rqtl = function(pheno_hists, kanto_analysis=NULL, mmoments_analysis=NULL) {
  phenodata_rqtl = sapply(pheno_hists, function(h) {
    c(mean=h$mean, var=h$var)
  })
  phenodata_rqtl = data.frame(t(phenodata_rqtl))
  phenodata_rqtl$noise = sqrt(phenodata_rqtl$var) / phenodata_rqtl$mean
  if (!is.null(kanto_analysis)) {
    phenodata_rqtl$mds1_k = kanto_analysis$mds$points[,1]
  }
  if (!is.null(mmoments_analysis)) {
    phenodata_rqtl$acp1_mm = mmoments_analysis$mds$rotation[,1]
  }
  return(phenodata_rqtl)
}



#' A Function That Preprocesses Genotypes of Individuals for R Package QTL.
#'
#' This function preprocesses genotypes of individuals for R package QTL. 
#' In geneodata, "2" stands for undetermined and is replaced by NA.
#' @param genodata A matrix describing the genotype of individuals.
#' @param bckg Vector of character describing the individuals as they are describe in genodata.
#' @export
preprocess_genodata_rqtl = function(genodata, bckg) {
  g = genodata
  prob_name =  as.character(g$prob_name)
  chromosome =  g$chromosome
  rec_fractions = genodata$rec_fractions
  g[g==2] = NA
  g = g[, bckg] + 1
  g = data.frame(g)
  g$chromosome = chromosome
  g$prob_name = prob_name
  g$rec_fractions = rec_fractions
  return(g)
}

#' A Function That Builds Kantorivitch Distance Matrix.
#'
#' This function builds Kantorivitch distance matrix from a list of histograms with all the same `breaks` value.
#' @param hists A list of histograms.
#' @export
build_kd_matrix = function(hists) {
  l = length(hists)
  bar = sapply(1:l, function(i) {
    foo = sapply(1:l, function(j) {
      if (i<j) {
        h1 = hists[[i]]
        h2 = hists[[j]]
        kantorovich(x=0, y=0, h1=h1, h2=h2)
      } else{
        0
      }
    })
    foo
  })
  bar = bar + t(bar)
  return(bar)
}
#' A Function That Builds Multivariate Moment Matrix.
#'
#' This function builds multivariate moment matrix from a list of histograms.
#' @param hists A list of histograms.
#' @param nb_moments An integer specifying the number of moments to compute.
#' @importFrom moments moment
#' @export
build_mmoments_matrix = function(hists, nb_moments=4) {
  l = length(hists)
  print(paste("Building mmoments... ", "(", nb_moments, ")", sep=""))
  bar = lapply(1:l, function(i) {
    c = hists[[i]]$cells
    moments = sapply(1:nb_moments, function(i){
      moment(c, order=i, central=i!=1)      
    })
    return(moments)
  })
  bar = do.call(rbind,bar)
  return(as.matrix(bar))
}

#' A Function That Builds Histograms From Single Cell Data.
#'
#' This function builds histograms form single cell data.
#' @param cells A list of vector of integer (single cell values).
#' @param bin_width An integer specifying the width of the bin to use.
#' @param nb_bin An integer specifying the number of to use.
#' @export
build_pheno_hists = function(cells, bin_width=NULL, nb_bin=100) {
  r = range(unlist(cells))
  if (is.null(bin_width)) {
    breaks = seq(r[1], r[2], length.out=nb_bin)
  } else {
    breaks = seq(r[1], r[2], by=bin_width)    
  }
  pheno_hists = lapply(cells, function(c){
    fl1 = c
    x = hist(fl1, plot = FALSE, breaks = breaks)
    x$mean = mean(fl1)
    x$var = var(fl1)
    x$nb_cells = length(fl1)
    x$cells = c
    x
  })
  return(pheno_hists)  
}

#' A Function That Computes Z-Score.
#'
#' This function returns z-score of serie.
#' @param dWs A vector of numeric.
#' @param p A numeric that specifies the quantile considered a noise.
#' @return the z-score 
#' @export
zscore = function(dWs, p=0.95) {
  best = max(dWs)
  q = quantile(dWs, probs=c(p))
  noise = dWs[dWs<q]
  zscore = (best - mean(noise))/sd(noise)
}


#' A Function That Computes Kantorovich Distance.
#'
#' This function computes Kantorovich distance betwenn two vector of numerics, or two histograms.
#' @param x A vector of numerics.
#' @param y A vector of numerics.
#' @param nbreaks The number of breaks used to build histograms.
#' @param lims vector of numerics of length 2 that specifies the range on which the histograms need to be built.
#' @param h1 An histogram.
#' @param h2 An histogram with the same breaks value as h1.
#' @return The kantorovich distance between two samples.
#' @examples
#' layout(matrix(1:3, 3))
#' x = rnorm(100000,0,1)
#' y = rnorm(100000,0,2)
#' lims = c(min(x, y),max(x, y))
#' H = hist(lims, plot = FALSE, breaks = 200);
#' tmp_breaks = H$breaks
#' h1 <-  hist(x, plot = FALSE, breaks = tmp_breaks);
#' h2 <-  hist(y, plot = FALSE, breaks = tmp_breaks);
#' diff = h1$density - h2$density
#' cumdiff = cumsum(diff)
#'
#' plot(h1$density, type="l", main="Distributions", col=2, lwd=3, xlab="", ylab="density", xaxt="n", yaxt="n")
#' lines(h2$density, type="l", col=4, lwd=3)
#' legend("topright", col=c(2,4),  legend=c("P", "Q"), lwd=2)
#' abline(h=0, lty=2)
#' axis(2, at=0)
#'
#' plot(diff, type="l", main="Difference", lwd=3, xlab="", ylab="diff", xaxt="n", yaxt="n")
#' abline(h=0, lty=2)
#' axis(2, at=0)
#'
#' plot(cumdiff, type="l", main="Cumulative Sum", lwd=3, xlab="", ylab="cumsum", xaxt="n", yaxt="n")
#' polygon(x= 1:length(cumdiff), y = cumdiff, col = "grey")
#' lines(cumdiff, type="l", lwd=3)
#' abline(h=0, lty=2)
#' axis(2, at=0)
#' @export
kantorovich = function(x, y, nbreaks=100, lims=NULL, h1=NULL, h2=NULL){
  if (is.null(h1)) {
    if (is.null(lims)) {
      lims = c(min(x, y),max(x, y))
    }
    H = hist(lims, plot = FALSE, breaks = nbreaks);
    tmp_breaks = H$breaks
    h1 <-  hist(x, plot = FALSE, breaks = tmp_breaks);
    h2 <-  hist(y, plot = FALSE, breaks = tmp_breaks);
  }
  binsize = h1$breaks[2] - h1$breaks[1]
  diff = h1$density - h2$density
  cumdiff = cumsum(diff);
  KD = binsize*binsize*sum(abs(cumdiff));
  return(KD)
}


#' A Function That Computes Kantorovich Distance between two 2D samples.
#'
#' This function computes Kantorovich distance betwenn two 2D samples.
#' @param x A vector of numerics.
#' @param y A vector of numerics.
#' @param nbreaks The number of breaks used to build histograms.
#' @param lims vector of numerics of length 2 that specifies the range on which the histograms need to be built.
#' @param h1 An histogram.
#' @param h2 An histogram with the same breaks value as h1.
#' @return The kantorovich distance between two 2D samples.
#' @examples
#' # # dataset
#' # x = cbind(rnorm(100000,0,1),rnorm(100000,0,1))
#' # y = cbind(rnorm(100000,0,2),rnorm(100000,0,2))
#' #
#' # # kantorovich2D algo
#' # nbreaks=32
#' # lims=c(min(x[,1], y[,1]),max(x[,1], y[,1]),min(x[,2], y[,2]), max(x[,2],y[,2]))
#' # k1 = MASS::kde2d(x[,1], x[,2], n=nbreaks, lims=lims)
#' # k2 = MASS::kde2d (y[,1], y[,2], n=nbreaks, lims=lims)
#' # binsizex = k1$x[2] - k1$x[1]
#' # binsizey = k1$y[2] - k1$y[1]
#' # diff = k1$z - k2$z
#' # cumsum2D = function(A){
#' #   t(apply(apply(A, 2, cumsum), 1, cumsum))
#' # }
#' # cumdiff = cumsum2D(diff)
#' # KD = binsizex*binsizey*binsizex*binsizey*sum(abs(cumdiff))
#' #
#' # # illustration
#' # layout(matrix(1:4, 2), respect=TRUE)
#' # persp(k1, phi = 30, theta = 20, d = 5)
#' # persp(k2, phi = 30, theta = 20, d = 5)
#' # kdiff = k1
#' # kdiff$z = diff
#' # persp(kdiff, phi = 30, theta = 20, d = 5)
#' # kcumdiff = k1
#' # kcumdiff$z = cumdiff
#' # persp(kcumdiff, phi = 30, theta = 20, d = 5)
#' @export
kantorovich2D = function(x, y, nbreaks=32, lims=NULL, k1=NULL, k2=NULL) {
  if (is.null(k1)) {
    if (is.null(lims)) {
      lims=c(min(x[,1], y[,1]),max(x[,1], y[,1]),min(x[,2], y[,2]), max(x[,2],y[,2]))
    }
    k1 = kde2d(x[,1], x[,2], n=nbreaks, lims=lims)
    k2 = kde2d(y[,1], y[,2], n=nbreaks, lims=lims)
  }  
  binsizex = k1$x[2] - k1$x[1]
  binsizey = k1$y[2] - k1$y[1]
  diff = k1$z - k2$z
  cumsum2D = function(A){
    t(apply(apply(A, 2, cumsum), 1, cumsum))
  }
  cumdiff = cumsum2D(diff)
  KD = binsizex*binsizey*binsizex*binsizey*sum(abs(cumdiff))
  return(KD)
}

#' A Workflow That Maps PTLs, QTLs, Returns involved Objects, Embed Caching Fetures.
#'
#' This function is a workflow that encapsulated many functions of the R package `ptlmapper`. 
#' It aggregates parameter values, inputs and outputs of the call to the R package `ptlmapper` function.
#' The resulting data structure could be cache on the file system. 
#' It is useful to massivelly save multiple call to `plt_scan` and `rqtl_launch` function, for example in a complexe design.  
#' @inheritParams preprocess_genodata
## #' @param genodata ...
## #' @param bckg ...
## #' @param minimal_proportion ...
#' @inheritParams build_pheno_hists
## #' @param cells ...
## #' @param bin_width ...
## #' @param nb_bin ...
#' @param nb_perm An integer that specifies the number of permuation to do.
#' @param nb_dim An integer that specifies the number of dimension of the MDS space to explore.
#' @param nb_moments An integer specifying the number of moments to compute.
#' @param errs A vector of integer (error) that will be used to compute threshold from the permutation test.
#' @param ptl_mapping_filename A character string that specifies the file to save the `ptl_mapping` results.
#' @param COMPUTE_KD_MATRIX A boolean that specifies if `kd_matrix` needs to be computed.
#' @param COMPUTE_MM_MATRIX A boolean that specifies if `mm_matrix` needs to be computed.
#' @param DO_KANTO A boolean that specifies if  `ptl_scan` with method="kanto" needs to be performed.
#' @param DO_MMOMENTS A boolean that specifies if `ptl_scan` with method="mmoment" needs to be performed.
#' @param DO_RQTL A boolean that specifies if `rqtl_launch` needs to be called.
#' @param CLEAN_OBJECT A boolean that specifies if `ptl_mapping` results needs to be cleaned.
#' @importFrom qtl jittermap
#' @importFrom qtl calc.genoprob
#' @importFrom qtl sim.geno
#' @importFrom qtl scanone
#' @export
ptl_mapping = function(
genodata, 
cells, 
bckg, 
nb_perm=20, 
bin_width=NULL, 
nb_dim=NULL, 
errs = 0.05, 
nb_bin=100, 
nb_moments=4,  
minimal_proportion=0.1, 
ptl_mapping_filename=NULL, 
COMPUTE_KD_MATRIX=TRUE, 
COMPUTE_MM_MATRIX=TRUE, 
DO_KANTO=TRUE, 
DO_MMOMENTS=TRUE, 
DO_RQTL=TRUE, 
CLEAN_OBJECT=FALSE
) {
  
  #################
  # cache feature # 
  #################
  if (!is.null(ptl_mapping_filename)) {
    if (file.exists(ptl_mapping_filename)) {
      ptl_mapping = readRDS(ptl_mapping_filename)
      return(ptl_mapping)
    } else {
      print(paste("Computing ", ptl_mapping_filename, "...", sep =""))
    }
  }
  ###############
  # pheno_hists #
  ###############
  pheno_hists = build_pheno_hists(cells, bin_width=bin_width, nb_bin=nb_bin)

  #############
  # kd_matrix #
  #############
  if (COMPUTE_KD_MATRIX | DO_KANTO) {
    kd_matrix = build_kd_matrix(pheno_hists)
  } else {
    kd_matrix = NULL    
  }

  #############
  # mm_matrix #
  #############
  if (COMPUTE_MM_MATRIX | DO_MMOMENTS) {
    mm_matrix = build_mmoments_matrix(pheno_hists, nb_moments=nb_moments)
  } else {
    mm_matrix = NULL
  }
  
  if (DO_KANTO | DO_MMOMENTS) {
    ################
    # genodata_ptl #
    ################
    genodata_ptl = preprocess_genodata(genodata, bckg)

    ##################
    # kanto_analysis #
    ##################
    if (DO_KANTO) {
      kanto_analysis = ptl_scan(kd_matrix, genodata_ptl, nb_perm=nb_perm, nb_dim=nb_dim, method="kanto")
    } else {
      kanto_analysis=NULL
    }      

    #####################
    # mmoments_analysis #
    #####################
    if (DO_MMOMENTS) {
      mmoments_analysis = ptl_scan(mm_matrix, genodata_ptl, nb_perm=nb_perm, nb_dim=nb_dim, method="mmoment")
    } else {
      mmoments_analysis=NULL
    }
    
  } else {
    kanto_analysis=NULL
    mmoments_analysis=NULL
    genodata_ptl=NULL
  }
  
  if (DO_RQTL) {
    #############
    # rqtl_data #
    #############
    rqtl_data = preprocess_data_rqtl(genodata, bckg, pheno_hists, kanto_analysis, mmoments_analysis)

    #################
    # rqtl_analysis #
    #################
    rqtl_analysis = rqtl_launch(rqtl_data, nb_perm=nb_perm, errs=errs)
  } else {
    rqtl_analysis=NULL
  }

  ###########
  # results # 
  ###########
  ptl_mapping = list(genodata=genodata, genodata_ptl=genodata_ptl, pheno_hists=pheno_hists, bckg=bckg, rqtl_analysis=rqtl_analysis, kanto_analysis=kanto_analysis, mmoments_analysis=mmoments_analysis, kd_matrix=kd_matrix, mm_matrix=mm_matrix)

  #################
  # cache feature # 
  #################
  if (!is.null(ptl_mapping_filename)) {
    if (!file.exists(dirname(ptl_mapping_filename))) {
      dir.create(path=dirname(ptl_mapping_filename), recursive=TRUE)
    }
    saveRDS(ptl_mapping, file=ptl_mapping_filename)
  }
    
  #################
  # Clean objects #
  #################
  if (CLEAN_OBJECT) {
    ptl_mapping$genodata = NULL
    ptl_mapping$genodata_ptl = NULL
    ptl_mapping$pheno_hists = NULL
    ptl_mapping$bckg = NULL    
    ptl_mapping$rqtl_analysis$pheno = NULL
    ptl_mapping$rqtl_analysis$geno = NULL
    ptl_mapping$pheno = NULL
  }

  return(ptl_mapping)
}





































#' A Function That Plots Encapsulated RQTL Outputs.
#'
#' This function plots encapsulated RQTL outputs.
#' @param ptl_mapping_result ...
#' @param main ...
#' @param ylim ...
#' @param which_pheno ...
#' @export
plot_rqtl = function(ptl_mapping_result, main="", ylim=NULL, which_pheno=NULL) {
  if (is.null(which_pheno)) {
    which_pheno = 1:length(ptl_mapping_result$rqtl_analysis$scan)
  }
  for (i in which_pheno) {
    rqtl_res = ptl_mapping_result$rqtl_analysis$scan[[i]]
    if (is.null(ylim)) {
      cur_ylim=c(0,max(rqtl_res$scan_output$lod))
    } else {
      cur_ylim=ylim
    }
    main = paste(names(ptl_mapping_result$rqtl_analysis$pheno)[i], main)
    # plot(rqtl_res$scan_output, main=main, ylim=c(0, max(c(rqtl_res$thres, rqtl_res$scan_output$lod))))
    plot(rqtl_res$scan_output, main=main, ylim=cur_ylim)
    abline(h=rqtl_res$thres, lty=2:4, col="grey")
    legend("topright", legend=rqtl_res$errs, lty=2:4, col="grey")
  }
}

#' A Function That Plots Encapsulated Wilks Scores.
#'
#' This function plots encapsulated Wilks scores.
#' @param ptl_mapping_result ...
#' @param main ... 
#' @param pa ... 
#' @param errs ... 
#' @param ylim ...
#' @param y_thres ...
#' @export
plot_wilks = function(ptl_mapping_result, main="", pa=NULL, errs = c(0.05, 0.01, 0.005), ylim=NULL, y_thres=NULL, method="kanto") {
  if (!is.null(ptl_mapping_result)) {
    axis_info = extract_axis_info(ptl_mapping_result)    
  } else {
    stop("ptl_mapping_result is NULL, can't get axis information.")
  }
  if (is.null(pa)) {
    if (method=="kanto") {
      pa = ptl_mapping_result$kanto_analysis
    } else {
      pa = ptl_mapping_result$mmoments_analysis
    }
  }
  pa$xs = axis_info$xs
  pa$chrs = axis_info$chrs
  if (is.null(ylim)) {
    cur_ylim=c(0, max(-log10(pa$Ws)[!is.infinite(-log10(pa$Ws))]))
  } else {
    cur_ylim=ylim
  }
  main = paste(main, " (", method, ")", sep="")
  plot(0,0, main=main, xlab="markers", ylab="-log10(F)", xaxt='n', col=0,
    xlim=c(min(pa$xs), max(pa$xs)), ylim=cur_ylim)
  for (chr in unique(pa$chrs)) {
    idx = which(pa$chrs == chr)
    if (length(idx) > 1) {
      lines(pa$xs[idx], -log10(pa$Ws)[idx], type="l", lw=2)
    } else {
      points(pa$xs[idx], -log10(pa$Ws)[idx], lw=2)
    }
  }
  axis(1, at=1:length(unique(pa$chrs)), labels=unique(pa$chrs))
  abline(v=(0:length(unique(pa$chrs))) + 0.5, lty=3, col="grey")
  if (length(pa$perm_Ws) > 1) {
    q = quantile(-log10(pa$perm_Ws),1-errs)
    abline(h=q, lty=2:length(errs), col="grey")
    legend("topright", legend=errs, lty=2:length(errs), col="grey")
  }
}

#' A Function That Plots Encapsulated Wilks Scores for MultiVariate Moment Approach.
#'
#' This function plots encapsulated Wilks scores for multiVariate moment approach.
#' @param ptl_mapping_result ...
#' @param main ... 
#' @export
plot_wilks_mm = function(ptl_mapping_result, main="", ...) {
  main = paste("Wilks on MM", main)
  plot_wilks(ptl_mapping_result=NULL, pa=ptl_mapping_result$mmoments_analysis, main=main, ...)
}

#' A Function That Plots Encapsulated Wilks Scores for Kantoritch Approach.
#'
#' This function plots encapsulated Wilks scores for Kantoritch approach.
#' @param ptl_mapping_result ...
#' @param main ... 
#' @export
plot_wilks_k = function(ptl_mapping_result, main="", ...) {
  main = paste("Wilks on Kanto", main)
  plot_wilks(ptl_mapping_result=ptl_mapping_result, main=main, ...)
}


#' A Function That Plots Empirical P-Values.
#'
#' This function plots empirical p-values.
#' @param ptl_mapping_result ...
#' @param main ... 
#' @param errs ... 
#' @export
plot_empirical_test = function(ptl_mapping_result, main="", errs = c(0.05, 0.01, 0.005)) {
  if (!is.null(ptl_mapping_result$kanto_analysis$perm_Ws)) {
    pa = ptl_mapping_result$kanto_analysis
    if (!is.null(ptl_mapping_result)) {
      axis_info = extract_axis_info(ptl_mapping_result)    
    } else {
      stop("ptl_mapping_result is NULL, can't get axis information.")
    }
    pa$xs = axis_info$xs
    pa$chrs = axis_info$chrs
    emp = sapply(pa$Ws, function(W){
      max(sum(pa$perm_Ws < W), 1) / length(pa$perm_Ws)
    })
    main = paste("Empical test", main)
    plot(0,0, main=main, xlab="markers", ylab="-log10(empirical)", xaxt='n', col=0,
      xlim=c(min(pa$xs), max(pa$xs)), ylim=c(0, max(-log10(emp)[!is.infinite(-log10(emp))])))
    for (chr in unique(pa$chrs)) {
      idx = which(pa$chrs == chr)
      if (length(idx) > 1) {
        lines(pa$xs[idx], -log10(emp)[idx], type="l", lw=2)
      } else {
        points(pa$xs[idx], -log10(emp)[idx], lw=2)
      }
    }
    axis(1, at=1:length(unique(pa$chrs)), labels=unique(pa$chrs))
    # abline(v=0:length(unique(pa$chrs)) + 0.5, col="grey", lty=3)
    abline(h=-log10(errs), lty=2:4, col="grey")
    legend("topright", legend=errs, lty=2:4, col="grey")
  }
}


#' A Function That Plots MDS.
#'
#' This function plots MDS.
#' @param ptl_mapping_result ...
#' @param col ... 
#' @param main ... 
#' @param mds ... 
#' @export
plot_mds = function(ptl_mapping_result, main="", col=NULL, mds=NULL, ...){
  if (is.null(col)) {
    col = 1
  }
  if (is.null(mds)) {
    mds = ptl_mapping_result$kanto_analysis$mds
  }
  p = mds$points
  plot(p[,1], p[,2], cex = 2, pch=16, col=adjustcolor(col,alpha=0.3), main=paste("MDS", main),
    xlab=paste("Mds1 (", round(mds$eig[1]/sum(mds$eig[1:(dim(mds$points)[2])]) * 100,2), "%)", sep=""),
    ylab=paste("Mds2 (", round(mds$eig[2]/sum(mds$eig[1:(dim(mds$points)[2])]) * 100,2), "%)", sep=""), ...)
}
#' A Function That Plots the Linear Disciminant Analysis.
#'
#' This function plots the linear disciminant analysis.
#' @param ptl_mapping_result ...
#' @param marker_names ...
#' @param main ... 
#' @param col ... 
#' @export
plot_can = function(ptl_mapping_result, marker_names, main="", col=NULL){
  col = na.omit(col)
  idx = which(names(ptl_mapping_result$kanto_analysis$genodata) == marker_names)
  can = ptl_mapping_result$kanto_analysis$cans[[idx]]
  nb_dim_mds = ptl_mapping_result$kanto_analysis$nb_dim
  # all_col = (can$scores$allele * 2) + 2
  main = paste("CAN@", marker_names, " #dim_can=", can$ndim, " #dim_mds=", nb_dim_mds, " ", main, sep="")
  if (ncol(can$scores) == 2 ){
    # plot(can$scores$Can1, pch=16, cex = 2, col=adjustcolor(col,alpha=0.3), main=main)
    boxplot(can$scores$Can1 ~ can$scores$allele, main=main)
    # ylab=paste("Can1 (", round(can$eigenvalues[1]/sum([1:(dim(mds$points)[2])]) * 100,2), "%)", sep=""),
  } else {
      plot(can$scores$Can1, can$scores$Can2, pch=16, cex = 2, col=adjustcolor(col,alpha=0.3), main=main)
  }
    # xlab=paste("Can1 (", round(mds$eig[1]/sum(mds$eig[1:(dim(mds$points)[2])]) * 100,2), "%)", sep=""),
    # ylab=paste("Can2 (", round(mds$eig[2]/sum(mds$eig[1:(dim(mds$points)[2])]) * 100,2), "%)", sep=""))
}

#' A Function That Plots Empirical P-Values.
#'
#' This function plots empirical p-values.
#' @param ptl_mapping_result ...
#' @param col ... 
#' @param main ... 
#' @export
plot_noise = function(ptl_mapping_result, main="", col=NULL){
  if (is.null(col)) {
    col = 1
  }
  phenodata_rqtl = ptl_mapping_result$rqtl_analysis$pheno
  plot(phenodata_rqtl$mean, sqrt(phenodata_rqtl$var) / phenodata_rqtl$mean, xlab = "mean", ylab = "noise", cex = 2, pch=16, col=adjustcolor(col,alpha=0.3), main=paste("Noise", main))
}












#' A Function Distribution of the Phenotypes.
#'
#' This function plots distribution of the phenotypes.
#' @param ptl_mapping_result ...
#' @param col ... 
#' @param main ... 
#' @param xlim ... 
#' @param ylim ... 
#' @param only_n_first ... 
#' @param KANTO_DENSITY ... 
#' @param MEAN_DENSITY ... 
#' @param ... ...
#' @export
plot_dist = function(ptl_mapping_result, col=NULL, main="", xlim=NULL, ylim=NULL, only_n_first=0, KANTO_DENSITY=TRUE, MEAN_DENSITY=FALSE, ...){
  h = ptl_mapping_result$pheno_hists
  if (is.null(col)) {
    col = rep(1,length(h))
  }
  if (is.null(ylim)) {
    ylim = c(0,max(sapply(h, function(i){max(i$density)})))
  }
  if (is.null(xlim)) {
    xlim=c(min(h[[1]]$breaks[-1]),max(h[[1]]$breaks[-1]))
  }
  plot(0,0, col=0, xlim=xlim, ylim=ylim, main=paste("Distrib.", main), xlab="", ylab="")
  if (only_n_first == 0) {
    only_n_first = length(h)
  }
  for (i in 1:only_n_first)   {
    if (KANTO_DENSITY) {
      lines(h[[i]]$breaks[-1], h[[i]]$density, col=adjustcolor(col[i],alpha=0.3))
    } else {
      lines(density(h[[i]]$cells), col=adjustcolor(col[i],alpha=0.3))
    }
  }
  if (!KANTO_DENSITY) {
    mean_density = function(h) {
      nb_c_tot = sum(sapply(1:length(h), function(i){
        length(h[[i]]$cells)
      }))
      nb_c_per_sample = ceiling(nb_c_tot / length(h))
      print(paste(nb_c_tot, "cells", nb_c_per_sample, "per sample"))
      cells = lapply(1:length(h), function(i) {
        c = h[[i]]$cells
        print(paste(length(c), "cells in", i))
        sample(rep(c, ceiling(nb_c_per_sample / length(c))), nb_c_per_sample)
      })
      d = density(unlist(cells))
      return(d)
    }
    if (MEAN_DENSITY) {
      for (c in na.omit(unique(col))) {
        # cells1 = sapply(which(col == c), function(i){
        #   h[[i]]$cells
        # })
        # d = density(unlist(cells1)
        d = mean_density(h[which(col == c)])
        lines(d, col=c, lwd=3, lty=2)
      }
    }
  }
}





#' A Function That Returns Best Markers Obtained in the Classical QTL Search.
#'
#' This function returns best markers obtained in the classical QTL search.
#' @param ptl_mapping_result ...
#' @param which_pheno ... 
#' @export
get_best_markers_rqtl = function(ptl_mapping_result, which_pheno) {
  tmp_scan_output = ptl_mapping_result$rqtl_analysis$scan[[which_pheno]]$scan_output
  tmp_scan_output = tmp_scan_output[!grepl("c*loc*",rownames(tmp_scan_output)),]
  best_marker_names = rownames(tmp_scan_output[which(tmp_scan_output$lod == max(tmp_scan_output$lod)),])
  lods = tmp_scan_output$lod
  return(tmp_scan_output[which (lods == max(lods)), ])
}


#' A Function That Returns Best Markers Obtained in the PTL Search.
#'
#' This function returns best markers obtained in the PTL search.
#' @param ptl_mapping_result ...
#' @param chr ... 
#' @export
get_best_markers_rptl = function(ptl_mapping_result, chr=NULL) {
  pa = ptl_mapping_result$kanto_analysis
  if (!is.null(ptl_mapping_result)) {
    axis_info = extract_axis_info(ptl_mapping_result)    
  } else {
    stop("ptl_mapping_result is NULL, can't get axis information.")
  }
  pa$xs = axis_info$xs
  pa$chrs = axis_info$chrs
  if (!is.null(chr)) {
    min_Ws = min(pa$Ws[pa$chrs==chr])
  } else {
    min_Ws = min(pa$Ws)
  }
  best_marker_names = names(pa$genodata[which(pa$Ws == min_Ws)])
  # Retro compatibility for previously cached ptl_mapping_result (sindataset)
  if (!("prob_name" %in% names(ptl_mapping_result$genodata))) {
      ptl_mapping_result$genodata$prob_name = ptl_mapping_result$genodata$RQTL.name
  }
  ret = ptl_mapping_result$genodata[which(ptl_mapping_result$genodata$prob_name %in% best_marker_names),c("chromosome","position")]
  ret$W = -log10(min_Ws)
  rownames(ret) = best_marker_names
  return(ret)
}

#' A Function That Converts Marker Names Into Color Index.
#'
#' This function converts marker names Into color index.
#' @param ptl_mapping_result ...
#' @param marker_names ...
#' @export
mn_2_col = function(ptl_mapping_result, marker_names) {
  as.numeric(ptl_mapping_result$kanto_analysis$genodata[,marker_names]) * 2 + 2
}




# get_confidence_interval_rqtl = function(ptl_mapping_result, which_pheno=1, chr=1, prob=0.95) {
#   tmp_scan_output = ptl_mapping_result$rqtl_analysis$scan[[which_pheno]]$scan_output
#   # scanoneboot(ptl_mapping_result$rqtl_analysis, chr=5, pheno.col=2)
#   get_pos_bp = function(li) {
#       sapply(rownames(li), function(n){
#       p = ptl_mapping_result$genodata[which(ptl_mapping_result$genodata$prob_name == n),"position"]
#       if (length(p) == 0) {
#         p=""
#       }
#       return(p)
#     })
#   }
#   li = lodint(tmp_scan_output, chr=chr, expandtomarkers=TRUE)
#   li = cbind(data.frame(li), get_pos_bp(li))
#   bi = bayesint(tmp_scan_output, chr=chr, prob=prob, expandtomarkers=TRUE)
#   bi = cbind(data.frame(bi), get_pos_bp(bi))
#   return(list(lodint = li, bayesint = bi))
# }

# plot_dim_scan = function(ptl_mapping_result, main="") {
#   main=paste("CAN #dim scan", main)
#   if (!is.null(ptl_mapping_result$kanto_analysis$eig_to_scan)) {
#     pa = ptl_mapping_result$kanto_analysis
#     z = pa$dim_scan_zscores
#     plot(pa$eig_to_scan, z, type="b", main=main, xlab="#dim", ylab="z-score")
#     abline(v=ptl_mapping_result$kanto_analysis$nb_dim, lty=2, col=2)
#   } else {
#     print("WARNING! CAN dimension have not been scanned." )
#   }
# }
#
#
#
# plot_genodata_map = function(ptl_mapping_result, main="") {
#   genome_matrix = t(as.matrix(ptl_mapping_result$genodata[ptl_mapping_result$bckg]))
#   library(matlab)
#   imagesc(genome_matrix, xlab="chromosome I", ylab="individuals", main=paste("Genodata Map", main))
# }
#
#
# plot_recomb_fract = function(ptl_mapping_result, main="") {
#   genome_matrix = t(as.matrix(ptl_mapping_result$genodata[ptl_mapping_result$bckg]))
#   plot(density(apply(genome_matrix[,-1] != genome_matrix[,-ncol(genome_matrix)], 2, sum) / nrow(genome_matrix)), main=paste("Recombination Fraction (teta)", main))
# }
