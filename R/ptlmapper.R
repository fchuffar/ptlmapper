#' A Function That Computes Z-Score.
#'
#' This function computes z-score.
#' @param dWs ...
#' @param p ...
#' @export
zscore = function(dWs, p=0.95) {
  best = max(dWs)
  q = quantile(dWs, probs=c(p))
  noise = dWs[dWs<q]
  zscore = (best - mean(noise))/sd(noise)
}


#' A Function That Computes Kantorovich Distance.
#'
#' This function computes Kantorovich distance.
#' @param x ...
#' @param y ...
#' @param nbreaks ...
#' @param lims ...
#' @param h1 ...
#' @param h2 ...
#' @export
kantorovich = structure(function#  Kantorovitch Distance
### Compute the kantorovich distance between two 1D samples
##author<< Gael Yvert,
( 
x, ##<<
y, ##<<
nbreaks=1000, ##<<
lims=NULL, ##<<
h1=NULL, ##<<
h2=NULL  ##<<
){
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
### Returns the kantorovich distance between two samples
  }, ex=function(){
    print("example")
})




#' A Function That Maps PTLs.
#'
#' This function maps PTLs.
#'
#' @param genodata ...
#' @param cells ...
#' @param bckg ...
#' @param nb_perm ...
#' @param bin_width ...
#' @param ptl_mapping_filename ...
#' @param types ...
#' @param nb_dim ...
#' @param errs ...
#' @param minimal_proportion ...
#' @param USE_BOT_FOR_EXPLORATION ...
#' @param COMPUTE_KD_MATRIX ...
#' @param DO_KANTO ...
#' @param DO_RQTL ...
#' @param DO_MMOMENTS ...
#' @param nb_bin ...
#' @param nb_moments ...
#' @param LIGHTWEIGHT ...
#' @importFrom candisc candisc
#' @importFrom qtl jittermap
#' @importFrom qtl calc.genoprob
#' @importFrom qtl sim.geno
#' @importFrom qtl scanone
#' @export
ptl_mapping = function(
genodata, 
cells, 
bckg, 
nb_perm=0, 
bin_width=NULL, 
ptl_mapping_filename=NULL, 
types=c("raw"), 
nb_dim=NULL, 
errs = c(0.05, 0.01, 0.005), 
minimal_proportion=0.1, 
USE_BOT_FOR_EXPLORATION=FALSE, 
COMPUTE_KD_MATRIX=TRUE, 
DO_KANTO=TRUE, 
DO_RQTL=TRUE, 
DO_MMOMENTS=TRUE, 
nb_bin=100, 
nb_moments=3,  
LIGHTWEIGHT=FALSE
) {
  if ("rec.fractions" %in% names(genodata)){
    names(genodata)[which(names(genodata) == "rec.fractions")] = "rec_fractions"
  }
  if ("RQTL.name" %in% names(genodata)){
    names(genodata)[which(names(genodata) == "RQTL.name")] = "prob_name"
  }
  
  if (!is.null(ptl_mapping_filename)) {
    if (file.exists(ptl_mapping_filename)) {
      ptl_mapping = readRDS(ptl_mapping_filename)
      return(ptl_mapping)
    } else {
      print(paste("Computing ", ptl_mapping_filename, "...", sep =""))
    }
  }
  ##
  rqtl_launch = function(rqtldata, nb_perm, errs) {
    ret_scan_output = list()
    for (i in 1:length(rqtldata$pheno)) {
      thres=NULL
      tmp_scan_output = scanone(rqtldata, method="imp", pheno.col=i)
      nb_perm = max(nb_perm, 1000)
      permtest = scanone(rqtldata, method="hk", n.perm=nb_perm, pheno.col=i)
      thres = summary(permtest, alpha=errs)
      # summary(ptl_mapping_result$rqtldata$scan[[3]]$scan_output , perms=ptl_mapping_result$rqtldata$scan[[3]]$permtest, alpha=0.2, pvalues=TRUE)
      ret_scan_output[[length(ret_scan_output) + 1]] = list(scan_output = tmp_scan_output, thres=thres, errs=errs, permtest=permtest)
    }
    return(ret_scan_output)
  }
  ##
  rqtl_builddata = function(genodata_rqtl, phenodata_rqtl, bckg) {
    rqtldata = list()
    class(rqtldata) = c("bc", "cross")
    rqtldata$geno = list()
    for (chr in unique(genodata_rqtl$chromosome)) {
      # print(chr)
      chr = as.character(chr)
      rqtldata$geno[[chr]] = list()
      rqtldata$geno[[chr]]$data = t(genodata_rqtl[ genodata_rqtl$chromosome == chr, bckg])
      colnames(rqtldata$geno[[chr]]$data) = genodata_rqtl[ genodata_rqtl$chromosome == chr, ]$prob_name
      # rqtldata$geno[[chr]]$data[rqtldata$geno[[chr]]$data == 2] = NA
      # rqtldata$geno[[chr]]$data = rqtldata$geno[[chr]]$data + 1
      rqtldata$geno[[chr]]$map = as.list(cumsum(genodata_rqtl[genodata_rqtl$chromosome == chr, ]$rec_fractions * 100))
      names(rqtldata$geno[[chr]]$map) = genodata_rqtl[genodata_rqtl$chromosome == chr, ]$prob_name
      rqtldata$geno[[chr]]$map = unlist(rqtldata$geno[[chr]]$map)
      class(rqtldata$geno[[chr]]) = "A"
    }
    rqtldata$pheno = data.frame(phenodata_rqtl)
    rqtldata = jittermap(rqtldata)
    # names(rqtldata$geno[[3]])
    rqtldata = calc.genoprob(rqtldata, step=10, err=0.01)
    # names(rqtldata$geno[[3]])
    rqtldata = sim.geno(rqtldata, step=10, n.draws=8, err=0.01)
    # names(rqtldata$geno[[3]])
    # rqtldata = argmax.geno(rqtldata, step=10, err=0.01)
    # names(rqtldata$geno[[3]])
    # rqtldata = calc.errorlod(rqtldata, err=0.01)
    return(rqtldata)
  }
  ## PTL and MMOMENTS
  get_wilks_score = function(data, all) {
    # Fucntion extract from candisc package to obtain return pvalue.
    ## http://www.philender.com/courses/multivariate/notes2/can1.html
    ## http://www.philender.com/courses/multivariate/notes3/manova.html
    # ptlWilks = function (eig) {
    #   prod(1/(1 + eig[1:length(eig)]))
    # }
    # W = ptlWilks(can$eigenvalues[1:(can$rank)])
    seqWilks <- function (eig, p, df.h, df.e) {
      p.full <- length(eig)
      result <- matrix(0, p.full, 4)
      m <- df.e + df.h - ( p.full + df.h + 1)/2
      for (i in seq(p.full)) {
        test <- prod(1/(1 + eig[i:p.full]))
        p <- p.full + 1 - i
        q <- df.h + 1 - i
        s <- p^2 + q^2 - 5
        if (s > 0) {
          s = sqrt(((p * q)^2 - 4)/s)
        } else {
          s = 1
        }
        df1 <- p * q
        df2 <- m * s - (p * q)/2 + 1
        result[i,] <- c(test, ((test^(-1/s) - 1) * (df2/df1)),
          df1, df2)
      }
      result <- cbind(result, pf(result[,2], result[,3], result[,4], lower.tail = FALSE))
      colnames(result) <- c("LR test stat", "approx F", "num Df", "den Df", "Pr(> F)")
      rownames(result) <- 1:p.full
      return(result)
    }
    data$allele = as.factor(all)
    data = data[!is.na(data$allele),]
    model_formula = paste("cbind(", paste(names(data)[-length(names(data))], collapse=", "), ") ~ allele")
    mod = lm(model_formula, data=data)
    can = candisc(mod, data=data)
    p  = can$rank
    eig = can$eigenvalues[1:p]
    df.h =  can$dfh
    df.e = can$dfe
    tests = seqWilks(eig, p, df.h, df.e)
    return(list(W=tests[1,5], can=can, tests=tests))
  }
  # common function for both conditions get_min_perm_ws
  get_min_perm_ws = function(genodata, data){
    nb_indiv = nrow(genodata)
    perm_index = sample(1:nb_indiv)
    perm_geno = genodata[perm_index,]
    perm_Ws = apply(perm_geno,2, function(all) {
      get_wilks_score(data, all)$W
    })
    ret = min(perm_Ws)
    return(ret)            
  }      

  ptl_scan = function(pheno_KD, genodata, nb_perm, nb_dim, USE_BOT_FOR_EXPLORATION) {
    # function used to compute the number of significants eigen values
    get_nb_eig_sign = function(pheno_KD) {
      mds =   cmdscale(pheno_KD,2,eig=TRUE)
      mds_eig = mds$eig[mds$eig > 0]
      nb_eig_sign = sum(mds_eig/sum(mds_eig) >= 1/length(mds$eig))
      return(nb_eig_sign)
    }
    # 
    get_min_perm_ws_with_dim_scan = function(pheno_KD, genodata){
      nb_indiv = nrow(genodata)
      perm_index = sample(1:nb_indiv)
      perm_geno = genodata[perm_index,]
      ptl_analysis = ptl_scan(pheno_KD, perm_geno, nb_perm=0, nb_dim=0, USE_BOT_FOR_EXPLORATION=FALSE)
      ret = min(ptl_analysis$Ws)
      return(ret)            
    }      
    #######
    # GO! #
    #######
    # How many Dimension ?
    NEED_TO_COMPUTE_MDS_AND_CAN = TRUE
    dim_scan = dim_scan_zscores = eig_to_scan = NULL
    if (is.null(nb_dim)) {
      nb_dim = max(2, get_nb_eig_sign(pheno_KD))
      print(paste("number of dimesion for can (nb_eig_sign):", nb_dim))
    } else if (nb_dim==0) {
      eig_to_scan = 2:max(2, get_nb_eig_sign(pheno_KD))
      if (USE_BOT_FOR_EXPLORATION) {
        library(bot)
        dim_scan = botapply(eig_to_scan, function(test_dim){
          test_dim = test_dim[[1]]
          print(paste("Testing nb_dim=", test_dim, "/", max(eig_to_scan), "...", sep=""))
          ptl_scan(pheno_KD, genodata, nb_perm=0, test_dim)
        })
      } else {
        dim_scan = apply(t(eig_to_scan), 2, function(test_dim){
          print(paste("Testing nb_dim=", test_dim, "/", max(eig_to_scan), "...", sep=""))
          ptl_scan(pheno_KD, genodata, nb_perm=0, test_dim)
        })
      }      
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
        p_scan$genodata = NULL
        p_scan$mds = NULL
        p_scan$cans = NULL
        p_scan
      })
    } else {
      print(paste("number of dimension for can (forced):", nb_dim))
    }

    # Wilks score for our genodata
    if (NEED_TO_COMPUTE_MDS_AND_CAN, SHOW_SCAN_PROGRESSION=TRUE) {      
      mds = cmdscale(pheno_KD,nb_dim,eig=TRUE)
      data = data.frame(mds$points)
      foo = apply(t(1:length(genodata)), 2, function(i) {
        if (i %% 50 == 1 & SHOW_SCAN_PROGRESSION) {print(paste("Marker ", i, "/", length(genodata), sep=""))}
        all = genodata[,i]
        get_wilks_score(data, all)
      })
      foo = do.call(rbind,foo)
      Ws = unlist(foo[,"W"])
      cans = foo[,"can"]
    }

    # Perform permutation
    if (nb_perm > 0) {
      if (USE_BOT_FOR_EXPLORATION) {
        library(bot)
        perm_Ws = botapply(1:nb_perm, function(p){
          p = p[[1]]
          print(paste("permutation", p , "/", nb_perm))
          return(get_min_perm_ws_with_dim_scan(pheno_KD, genodata))
        })
        perm_Ws = unlist(perm_Ws)
      } else {
        perm_Ws = sapply(1:nb_perm, function(p){
          if (p %% 10 == 1) {print(paste("permutation", p , "/", nb_perm))}
            return(get_min_perm_ws_with_dim_scan(pheno_KD, genodata))
        })
      }
    } else {
      perm_Ws = NULL
    }
    
    return(list(
                genodata=genodata, 
                mds=mds, 
                Ws=Ws, 
                cans=cans, 

                eig_to_scan=eig_to_scan,
                dim_scan=dim_scan, 
                dim_scan_zscores=dim_scan_zscores, 
                nb_dim=nb_dim, 

                nb_perm=nb_perm, 
                perm_Ws=perm_Ws
                ))
  }
  mmoments_scan = function(pheno_moments, genodata, nb_perm, nb_dim, USE_BOT_FOR_EXPLORATION) {
    # function used to compute the number of significants eigen values
    get_nb_eig_sign = function(pheno_moments) {
      # pheno_moments = as.matrix(indiv[601:800,9:18])
      pca = prcomp(t(pheno_moments), scale=TRUE)
      pcavar = pca$sdev * pca$sdev
      # round(pcavar/sum(pcavar), 3) * 100
      # plot(pcavar/sum(pcavar), type = "h", lwd = 8, ylab = "Proportion of variance explained", xlab = "Principal component");
      # abline(h=1/length(pcavar), col = 2, lty = 1, lwd = 2);
      nb_eig_sign = sum(pcavar/sum(pcavar) >= 1/length(pcavar))
      return(nb_eig_sign)
    }
    # 
    get_min_perm_ws_with_dim_scan = function(pheno_moments, genodata){
      nb_indiv = nrow(genodata)
      perm_index = sample(1:nb_indiv)
      perm_geno = genodata[perm_index,]
      ptl_analysis = mmoments_scan(pheno_moments, perm_geno, nb_perm=0, nb_dim=0, USE_BOT_FOR_EXPLORATION=FALSE)
      ret = min(ptl_analysis$Ws)
      return(ret)            
    }      
    #######
    # GO! #
    #######
    # How many Dimension ?
    NEED_TO_COMPUTE_MDS_AND_CAN = TRUE
    dim_scan = dim_scan_zscores = eig_to_scan = NULL
    if (is.null(nb_dim)) {
      nb_dim = max(2, get_nb_eig_sign(pheno_moments))
      print(paste("number of dimension for can (nb_eig_sign):", nb_dim))
    } else if (nb_dim==0) {
      eig_to_scan = 2:max(2, get_nb_eig_sign(pheno_moments))
      if (USE_BOT_FOR_EXPLORATION) {
        library(bot)
        dim_scan = botapply(eig_to_scan, function(test_dim){
          test_dim = test_dim[[1]]
          print(paste("Testing nb_dim=", test_dim, "/", max(eig_to_scan), "...", sep=""))
          mmoments_scan(pheno_moments, genodata, nb_perm=0, test_dim)
        })
      } else {
        dim_scan = apply(t(eig_to_scan), 2, function(test_dim){
          print(paste("Testing nb_dim=", test_dim, "/", max(eig_to_scan), "...", sep=""))
          mmoments_scan(pheno_moments, genodata, nb_perm=0, test_dim)
        })
      }
      dim_scan_zscores = sapply(dim_scan, function(d){
        return(zscore(-log10(d$Ws)))
      })
      best_zscore_idx = which(dim_scan_zscores == max(dim_scan_zscores))[1]
      nb_dim = (eig_to_scan)[best_zscore_idx]
      print(paste("number of dimension for can (best zscore):", nb_dim))
      # Do not recalculated MDS nor CAN
      pca = dim_scan[[best_zscore_idx]]$mds
      Ws = dim_scan[[best_zscore_idx]]$Ws
      cans = dim_scan[[best_zscore_idx]]$cans
      NEED_TO_COMPUTE_MDS_AND_CAN = FALSE
      # Clean objects ?
      dim_scan = lapply(dim_scan, function(p_scan){
        p_scan$genodata = NULL
        p_scan$mds = NULL
        p_scan$cans = NULL
        p_scan
      })
    } else {
      print(paste("number of dimension for can (forced):", nb_dim))
    }

    # Wilks score for our genodata
    if (NEED_TO_COMPUTE_MDS_AND_CAN) {      
      pca = prcomp(t(pheno_moments), scale=TRUE)
      data = data.frame(pca$rotation[,1:nb_dim])
      foo = apply(t(1:length(genodata)), 2, function(i) {
        if (i %% 50 == 1) {print(paste("Marker ", i, "/", length(genodata), sep=""))}
        all = genodata[,i]
        get_wilks_score(data, all)
      })
      foo = do.call(rbind,foo)
      Ws = unlist(foo[,"W"])
      cans = foo[,"can"]
    }

    # Perform permutation
    if (nb_perm > 0) {
      if (USE_BOT_FOR_EXPLORATION) {
        library(bot)
        perm_Ws = botapply(1:nb_perm, function(p){
          p = p[[1]]
          print(paste("permutation", p , "/", nb_perm))
          # return(get_min_perm_ws(genodata, data))
          return(get_min_perm_ws_with_dim_scan(pheno_moments, genodata))
        })
        perm_Ws = unlist(perm_Ws)
      } else {
        perm_Ws = sapply(1:nb_perm, function(p){
          if (p %% 10 == 1) {print(paste("permutation", p , "/", nb_perm))}
            # return(get_min_perm_ws(genodata, data))
            return(get_min_perm_ws_with_dim_scan(pheno_moments, genodata))
        })
      }
    } else {
      perm_Ws = NULL
    }
    
    return(list(Ws=Ws, perm_Ws=perm_Ws, nb_perm=nb_perm, genodata=genodata, mds=pca, cans=cans, dim_scan=dim_scan, nb_dim=nb_dim, dim_scan_zscores=dim_scan_zscores, eig_to_scan=eig_to_scan))
  }
  ##
  get_KD_matrix = function(cur_hists) {
    l = length(cur_hists)
    print("Building KD...")
    bar = lapply(1:l, function(i) {
      foo = lapply(1:l, function(j) {
        if (i<j) {
          h1 = cur_hists[[i]]
          h2 = cur_hists[[j]]
          kantorovich(x=0, y=0, h1=h1, h2=h2)
        } else{
          0
        }
      })
      foo = do.call(rbind,foo)
      foo
    })
    bar = do.call(cbind,bar)
    bar = bar + t(bar)
    return(as.matrix(bar))
  }
  ##
  get_mmoments_matrix = function(cur_hists, nb_moments=nb_moments) {
    l = length(cur_hists)
    print(paste("Building mmoments... ", "(", nb_moments, ")", sep=""))
    bar = lapply(1:l, function(i) {
      c = cur_hists[[i]]$cells
      moments = sapply(1:nb_moments, function(i){
        moment(c, order=i, central=i!=1)      
      })
      return(moments)
    })
    bar = do.call(rbind,bar)
    return(as.matrix(bar))
  }
  ###############
  # pheno_hists #
  ###############
  r = lapply(cells, function(c) {
    range(c)
  })
  r = do.call(rbind,r)
  if (!is.null(bin_width)) {
    tmp_breaks = floor(min(r[,1])/bin_width):ceiling(max(r[,2])/bin_width) * bin_width
  # } else if (!any(unique(unlist(cells)) != as.integer(unique(unlist(cells))))) {
  #   tmp_breaks = min(r[,1]):max(r[,2])
  } else {
    tmp_breaks = (0:nb_bin) * (max(r[,2]) - min(r[,1])) / nb_bin + min(r[,1])
  }

  pheno_hists_raw = lapply(cells, function(c){
    fl1 = c
    x = hist(fl1, plot = FALSE, breaks = tmp_breaks)
    x$mean = mean(fl1)
    x$var = var(fl1)
    x$nb_cells = length(fl1)
    x$cells = c
    x
  })

  pheno_hists = list(raw=pheno_hists_raw)

  # if ("centered" %in% types) {
  #   cells_centered = lapply(cells, function(c) {
  #     fl1 = c
  #     m = round(mean(fl1))
  #     c = fl1 - m
  #     c
  #   })
  #   # do.call()
  #
  #   r_centered = lapply(cells_centered, function(c) {
  #     range(c)
  #   })
  #   r_centered = do.call(rbind,r_centered)
  #   if (!is.null(bin_width)) {
  #     tmp_breaks_centered = floor(min(r_centered[,1])/bin_width):ceiling(max(r_centered[,2])/bin_width) * bin_width
  #   } else {
  #     tmp_breaks_centered = (0:nb_bin) * (max(r_centered[,2]) - min(r_centered[,1])) / nb_bin + min(r_centered[,1])
  #   }
  #
  #   pheno_hists_centered =  lapply(cells_centered, function(c){
  #     fl1 = c
  #     m = mean(fl1)
  #     x = hist(fl1, plot = FALSE, breaks = tmp_breaks_centered)
  #     x$mean = mean(fl1)
  #     x$var = var(fl1)
  #     x$nb_cells = length(fl1)
  #     x$cells = c
  #     x
  #   })
  #   pheno_hists[["centered"]] = pheno_hists_centered
  # }

  if (COMPUTE_KD_MATRIX) {
    kd_matrix = list()
    for (type in types) {
      kd_matrix[[type]] = get_KD_matrix(pheno_hists[[type]])
    }  
  } else {
    kd_matrix = NULL    
  }
  
  if (DO_KANTO | DO_MMOMENTS) {
    # ################
    # # genodata_ptl #
    # ################
    genodata_ptl = genodata
    colnames =  as.character(genodata_ptl$prob_name)
    genodata_ptl =  data.frame(cbind(bckg, t(genodata_ptl[, bckg])), stringsAsFactors=FALSE)
    names(genodata_ptl) = c("spores", colnames)
    # Check genodata
    genodata_ptl[genodata_ptl=="2"] = NA
    keep = sapply(1:ncol(genodata_ptl), function(i) {
      col = genodata_ptl[bckg, i]
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
    genodata_ptl = genodata_ptl[, keep]
    gen_only = genodata_ptl[bckg, ]
    # 

    get_axis_info = function(analysis, delta = 0.6) {
      sub_genodata = genodata[genodata$prob_name %in% names(genodata_ptl),]
      analysis$chrs = sub_genodata$chromosome        
      analysis$xs = unlist(
      lapply(1:length(unique(analysis$chrs)), function(i) {
        chr = unique(analysis$chrs)[i]
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
      return(analysis)
    }

    if (DO_KANTO) {
      ptl_analysis = list()
      for (type in types) {
        ptl_analysis[[type]] = ptl_scan(kd_matrix[[type]], gen_only, nb_perm=nb_perm, nb_dim=nb_dim, USE_BOT_FOR_EXPLORATION=USE_BOT_FOR_EXPLORATION)
        ptl_analysis[[type]] = get_axis_info(ptl_analysis[[type]])
      }
    } else {
      ptl_analysis=NULL
    }      
    if (DO_MMOMENTS) {
      library(moments)
      mmoments_analysis = list()
      for (type in types) {
        pheno_moments = get_mmoments_matrix(pheno_hists[[type]], nb_moments=nb_moments)
        mmoments_analysis[[type]] = mmoments_scan(pheno_moments, gen_only, nb_perm=nb_perm, nb_dim=nb_dim, USE_BOT_FOR_EXPLORATION=USE_BOT_FOR_EXPLORATION)
        mmoments_analysis[[type]] = get_axis_info(mmoments_analysis[[type]])
        mmoments_analysis[[type]]$pheno_moments = pheno_moments
      }
    } else {
      mmoments_analysis=NULL
    }
  } else {
    ptl_analysis=NULL
    mmoments_analysis=NULL
    genodata_ptl=NULL
  }

  if (DO_RQTL) {
    #################
    # genodata_rqtl #
    #################
    genodata_rqtl = genodata

    prob_name =  as.character(genodata_rqtl$prob_name)
    chromosome =  genodata_rqtl$chromosome
    rec_fractions = genodata$rec_fractions

    genodata_rqtl = genodata_rqtl[, bckg] + 1
    genodata_rqtl[genodata_rqtl==3] = NA
    genodata_rqtl = data.frame(genodata_rqtl)

    genodata_rqtl$chromosome = chromosome
    genodata_rqtl$prob_name = prob_name
    genodata_rqtl$rec_fractions = rec_fractions

    # phenodata_rqtl
    phenodata_rqtl = lapply(pheno_hists[["raw"]], function(h) {
      c(h$mean, h$var)
    })
    phenodata_rqtl = data.frame(do.call(rbind, phenodata_rqtl))
    names(phenodata_rqtl) = c("mean", "var")
    phenodata_rqtl$noise = sqrt(phenodata_rqtl$var) / phenodata_rqtl$mean

    if (DO_KANTO) {
      phenodata_rqtl$mds1_k = ptl_analysis[[type]]$mds$points[,1]
    }
    if (DO_MMOMENTS) {
      phenodata_rqtl$acp1_mm = mmoments_analysis[[type]]$mds$rotation[,1]
    }

    rqtldata = rqtl_builddata(genodata_rqtl, phenodata_rqtl, bckg)
    rqtldata$scan = rqtl_launch(rqtldata, nb_perm = nb_perm, errs=errs)
  } else {
    rqtldata=NULL
  }

  #
  ptl_mapping = list(genodata=genodata, genodata_ptl=genodata_ptl, pheno_hists=pheno_hists, bckg=bckg, rqtldata=rqtldata, ptl_analysis=ptl_analysis, mmoments_analysis=mmoments_analysis, kd_matrix=kd_matrix)
  if (!is.null(ptl_mapping_filename)) {
    if (!file.exists(dirname(ptl_mapping_filename))) {
      dir.create(path=dirname(ptl_mapping_filename), recursive=TRUE)
    }
    saveRDS(ptl_mapping, file=ptl_mapping_filename)
  }
    
  # Clean objects ?
  if (LIGHTWEIGHT) {
    ptl_mapping$genodata = NULL
    ptl_mapping$genodata_ptl = NULL
    ptl_mapping$pheno_hists = NULL
    ptl_mapping$bckg = NULL    
    ptl_mapping$rqtldata$pheno = NULL
    ptl_mapping$rqtldata$geno = NULL
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
    which_pheno = 1:length(ptl_mapping_result$rqtldata$scan)
  }
  for (i in which_pheno) {
    rqtl_res = ptl_mapping_result$rqtldata$scan[[i]]
    if (is.null(ylim)) {
      cur_ylim=c(0,max(rqtl_res$scan_output$lod))
    } else {
      cur_ylim=ylim
    }
    main = paste(names(ptl_mapping_result$rqtldata$pheno)[i], main)
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
#' @param type ...  
#' @param pa ... 
#' @param errs ... 
#' @param ylim ...
#' @param y_thres ...
#' @export
plot_wilks = function(ptl_mapping_result, main="", type="raw", pa=NULL, errs = c(0.05, 0.01, 0.005), ylim=NULL, y_thres=NULL) {
  if (is.null(pa)) {
    pa = ptl_mapping_result$ptl_analysis[[type]]
  } else if (!is.null(ptl_mapping_result)){
    pa$xs = ptl_mapping_result$ptl_analysis$raw$xs
    pa$chrs = ptl_mapping_result$ptl_analysis$raw$chrs
  }
  if (is.null(ylim)) {
    cur_ylim=c(0, max(-log10(pa$Ws)[!is.infinite(-log10(pa$Ws))]))
  } else {
    cur_ylim=ylim
  }
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
  plot_wilks(ptl_mapping_result=NULL, pa=ptl_mapping_result$mmoments_analysis[["raw"]], main=main, ...)
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
#' @param type ...  
#' @param errs ... 
#' @export
plot_empirical_test = function(ptl_mapping_result, main="", type="raw", errs = c(0.05, 0.01, 0.005)) {
  if (!is.null(ptl_mapping_result$ptl_analysis[[type]]$perm_Ws)) {
    pa = ptl_mapping_result$ptl_analysis[[type]]
    emp = sapply(pa$Ws, function(W){
      max(sum(pa$perm_Ws < W), 1) / length(pa$perm_Ws)
    })
    main = paste("Empical test", type, main)
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
#' @param type ...  
#' @param mds ... 
#' @export
plot_mds = function(ptl_mapping_result, main="", type="raw", col=NULL, mds=NULL, ...){
  if (is.null(col)) {
    col = 1
  }
  if (is.null(mds)) {
    mds = ptl_mapping_result$ptl_analysis[[type]]$mds
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
#' @param type ...  
#' @param col ... 
#' @export
plot_can = function(ptl_mapping_result, marker_names, main="", type="raw", col=NULL){
  col = na.omit(col)
  idx = which(names(ptl_mapping_result$ptl_analysis[[type]]$genodata) == marker_names)
  can = ptl_mapping_result$ptl_analysis[[type]]$cans[[idx]]
  nb_dim_mds = ptl_mapping_result$ptl_analysis[[type]]$nb_dim
  # all_col = (can$scores$allele * 2) + 2
  main = paste("CAN@", marker_names, " ", type, " #dim_can=", can$ndim, " #dim_mds=", nb_dim_mds, " ", main, sep="")
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
  phenodata_rqtl = ptl_mapping_result$rqtldata$pheno
  plot(phenodata_rqtl$mean, sqrt(phenodata_rqtl$var) / phenodata_rqtl$mean, xlab = "mean", ylab = "noise", cex = 2, pch=16, col=adjustcolor(col,alpha=0.3), main=paste("Noise", main))
}












#' A Function Distribution of the Phenotypes.
#'
#' This function plots distribution of the phenotypes.
#' @param ptl_mapping_result ...
#' @param col ... 
#' @param main ... 
#' @param type ...  
#' @param xlim ... 
#' @param ylim ... 
#' @param only_n_first ... 
#' @param KANTO_DENSITY ... 
#' @param MEAN_DENSITY ... 
#' @param ... ...
#' @export
plot_dist = function(ptl_mapping_result, col=NULL, main="", type="raw", xlim=NULL, ylim=NULL, only_n_first=0, KANTO_DENSITY=TRUE, MEAN_DENSITY=FALSE, ...){
  h = ptl_mapping_result$pheno_hists[[type]]
  if (is.null(col)) {
    col = rep(1,length(h))
  }
  if (is.null(ylim)) {
    ylim = c(0,max(sapply(h, function(i){max(i$density)})))
  }
  if (is.null(xlim)) {
    xlim=c(min(h[[1]]$breaks[-1]),max(h[[1]]$breaks[-1]))
  }
  plot(0,0, col=0, xlim=xlim, ylim=ylim, main=paste("Distrib.", type, main), xlab="", ylab="")
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
  tmp_scan_output = ptl_mapping_result$rqtldata$scan[[which_pheno]]$scan_output
  tmp_scan_output = tmp_scan_output[!grepl("c*loc*",rownames(tmp_scan_output)),]
  best_marker_names = rownames(tmp_scan_output[which(tmp_scan_output$lod == max(tmp_scan_output$lod)),])
  lods = tmp_scan_output$lod
  return(tmp_scan_output[which (lods == max(lods)), ])
}


#' A Function That Returns Best Markers Obtained in the PTL Search.
#'
#' This function returns best markers obtained in the PTL search.
#' @param ptl_mapping_result ...
#' @param type ...  
#' @param chr ... 
#' @export
get_best_markers_rptl = function(ptl_mapping_result, type="raw", chr=NULL) {
  p = ptl_mapping_result$ptl_analysis[[type]]
  if (!is.null(chr)) {
    min_Ws = min(p$Ws[p$chrs==chr])
  } else {
    min_Ws = min(p$Ws)
  }
  best_marker_names = names(p$genodata[which(p$Ws == min_Ws)])
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
  as.numeric(ptl_mapping_result$ptl_analysis[["raw"]]$genodata[,marker_names]) * 2 + 2
}




# get_confidence_interval_rqtl = function(ptl_mapping_result, which_pheno=1, chr=1, prob=0.95) {
#   tmp_scan_output = ptl_mapping_result$rqtldata$scan[[which_pheno]]$scan_output
#   # scanoneboot(ptl_mapping_result$rqtldata, chr=5, pheno.col=2)
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
#   if (!is.null(ptl_mapping_result$ptl_analysis$raw$eig_to_scan)) {
#     pa = ptl_mapping_result$ptl_analysis$raw
#     z = pa$dim_scan_zscores
#     plot(pa$eig_to_scan, z, type="b", main=main, xlab="#dim", ylab="z-score")
#     abline(v=ptl_mapping_result$ptl_analysis$raw$nb_dim, lty=2, col=2)
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
