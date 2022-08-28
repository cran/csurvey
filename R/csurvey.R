#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
csvy <- function(formula, design, subset = NULL, family = stats::gaussian(),
                 nD = NULL, level = 0.95, n.mix = 100L, test = TRUE,...){
  cl <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized \n")
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match("formula", names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, model.frame(design), parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")

  if (length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nm))
      names(y) <- nm
  }
  shapes_add <- NULL
  xnms_add <- NULL
  nums_add <- NULL
  block.ave.lst <- vector("list", length = (ncol(mf) - 1))
  block.ord.lst <- vector("list", length = (ncol(mf) - 1))
  cnms <- colnames(mf)
  for (i in 2:ncol(mf)) {
    block.ord.lst[[i - 1]] <- -1
    block.ave.lst[[i - 1]] <- -1
    if (is.numeric(attributes(mf[, i])$shape)) {
      shapes_add <- c(shapes_add, attributes(mf[, i])$shape)
      xnms_add <- c(xnms_add, attributes(mf[, i])$nm)
      if ((attributes(mf[, i])$categ == "block.ord")) {
        block.ord.lst[[i - 1]] <- attributes(mf[, i])$order
      }
      if ((attributes(mf[, i])$categ == "block.ave")) {
        block.ave.lst[[i - 1]] <- attributes(mf[, i])$order
      }
    } else if (is.null(attributes(mf[, i])$shape)) {
      shapes_add <- c(shapes_add, 0)
      xnms_add <- c(xnms_add, cnms[i])
    }
  }
  xid_add <- 2:ncol(mf)
  #check more
  xmat_add <- dplyr::select(mf, all_of(xid_add))
  #xmat_add <- mf[, -1, drop = FALSE]
  #how to handle ...?
  fit <- csvy.fit(design = design, fo = formula, subset = subset, family = family, M = nD, xm = xmat_add, xnms_add = xnms_add,
                     sh = shapes_add, ynm = ynm, block.ave.lst = block.ave.lst, block.ord.lst = block.ord.lst,
                     level = level, n.mix = n.mix, test = test,...)
  #classes <- if (inherits(design, "svyrep.design")) {
  #  union(class(fit$ans.unc), c("svrepglm", "svyglm"))
  #} else {
  #  union(class(fit$ans.unc), "svyglm")
  #}
  classes <- if (inherits(design, "svyrep.design")) {
    c("svrepglm", "svyglm")
  } else {
    "svyglm"
  }
  structure(c(fit, list(call = cl, n.mix = n.mix, xmat_add = xmat_add, xnms_add = xnms_add,
                        shapes_add = shapes_add, ynm = ynm, family = family, terms = mt)),
            class = c("csvy", "cgam", classes))
}


#---------------------------------------------------------------------------------------------------------------------------------------
#main method of csvy
#---------------------------------------------------------------------------------------------------------------------------------------
csvy.fit <- function(design = NULL, fo = NULL, subset = NULL, family = stats::gaussian(), M = NULL, xm = NULL, xnms_add = NULL,
                     sh = 1, ynm = NULL, block.ave.lst = NULL, block.ord.lst = NULL,
                     level = 0.95, n.mix = 100L, test = TRUE,...) {
  #----------------------------------------------------------------------------------
  #create parameters for coneA, ID, w, nd, etc
  #----------------------------------------------------------------------------------
  #print (list(...))
  #print (sh)
  #print (head(xm))
  #learn from surveyrep
  subset <- substitute(subset)
  subset <- eval(subset, design$variables, parent.frame())
  if (any(is.na(subset))) {
    stop("subset must not contain NA values \n")
  }
  if (!is.null(subset)) {
    design <- design[subset, ]
  }

  #----------------------------------------------------------------------------------
  #local helper function: transfor x's into factors and code them with integers
  #----------------------------------------------------------------------------------
  # make_factor <- function(vec) {
  #   #if(!is.factor(vec)){
  #   vec_fac <- as.factor(vec)
  #   vec_fac_lvls <- suppressWarnings(as.numeric(levels(vec_fac))[vec_fac])
  #   if(anyNA(vec_fac_lvls)){
  #     rslt <- as.factor(as.integer(vec_fac))
  #   } else {
  #     rslt <- as.factor(vec_fac_lvls)
  #   }
  #   return(rslt)
  #   #}
  # }
  #xm <- apply(xm, 2, function(e) make_factor(e))
  Ds <- map_dbl(xm, .f = function(.x) length(unique(.x)))
  Nd <- length(Ds)
  xm.red <- unique(xm)
  xvs2 <- map(xm.red, .f = function(.x) sort(unique(.x)))
  grid2 <- expand.grid(xvs2)
  colnames(grid2) <- xnms_add
  df <- design$variables
  ds <- design
  xm.df <- xm
  if (Nd == 1) {
    ID <- apply(xm, 1, function(elem) grid2[apply(grid2, 1, function(gi) all(gi == elem)), ])
  }
  if (Nd > 1) {
    ID <- apply(xm, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
  }
  ID <- as.ordered(ID)
  uds <- unique(ID)
  domain.ord <- order(uds)
  uds <- sort(uds) #avoid sort for factors not coded as numbers? require users specify x's by ordered numbers
  obs_cells <- uds

  #print (M)
  #Max_obs <- max(obs_cells) #Max_obs is factor -> gives NA when compared with M
  Tot_obs <- length(obs_cells)
  #print (obs_cells)
  #print (Max_obs)
  if (is.null(M)) {
    stop("nD (total number of domains) is not provided! \n")
    #M <- max(obs_cells)
  } else if (M %% 1 != 0 | M <= 1) {
    stop("nD (total number of domains) must be a positive integer > 1! \n")
    #M <- max(obs_cells)
  } else if (M < Tot_obs) {
    stop(paste("There are at least", Tot_obs, "domains! nD (total number of domains) should be >= ",
            Tot_obs, sep = " ", "\n"))
  }
  M_0 <- M

  weights <- stats::weights(design, "sampling")
  # Nhat and stratawt will follow the order: 1 to M
  if(family$family == "gaussian"){
    y <- df[[ynm]]
  } else {
    fo <- update(fo, ~. -1) #handle the case when x is ordered; remove the intercept; do nothing if fo already includes -1
    ans.unc <- suppressWarnings(svyglm(formula = fo, design = design, family = family,...))
    y <- ans.unc$y
  }

  # n <- nrow(grid2)
  # ys <- vector("list", length = M)
  # stratawt <- vector("list", length = M)
  # Nhat <- 1:M * 0
  # nd <- 1:M * 0
  # uds_num <- as.numeric(levels(uds))[uds] #some factors are not coded as numbers....
  # for (udi in uds_num) {
  #   ps <- which(ID %in% udi)
  #   ysi <- y[ps]
  #   wi <- weights[ps]
  #   Nhi <- sum(wi)
  #   Nhat[udi] <- Nhi
  #   nd[udi] <- length(wi)
  #   ys[[udi]] <- ysi
  #   stratawt[[udi]] <- wi
  # }

  n <- nrow(grid2)
  ys <- vector("list", length = M)
  stratawt <- vector("list", length = M)
  Nhat <- 1:M*0
  nd <- 1:M*0
  uds_num <- as.numeric(levels(uds))[uds] #some factors don't have numeric levels...
  empty_cells <- which(!1:M %in% uds)

  ne <- length(empty_cells)
  zeros_ps <- empty_cells
  iter <- 1
  for(i in 1:M){
    if(!i %in% empty_cells){
      udi <- uds[iter]
      ps <- which(ID %in% udi)
      #ps_lst[[i]] = ps
      #print (ps)
      ysi <- y[ps]
      wi <- weights[ps]
      Nhi <- sum(wi)
      Nhat[i] <- Nhi
      nd[i] <- length(wi)
      ys[[i]] <- ysi
      stratawt[[i]] <- wi
      iter <- iter + 1
    }
  }

  w <- Nhat / sum(Nhat)
  w_0 <- w
  obs_cells <- which(nd > 0)
  empty_cells <- if (any(nd == 0)) which(nd == 0) else NULL
  ne <- length(empty_cells)
  zeros_ps <- empty_cells
  zeros <- ne

  nd_ne <- nd[empty_cells]
  small_cells <- which(nd >= 1 & nd <= 10)
  nsm <- length(small_cells)

  if (!is.null(empty_cells)) {
    M <- M - length(empty_cells)
    w <- w[-empty_cells]
  }
  imat <- diag(M)
  # define weight mat and its inverse here to avoid repeating it
  # need to consider zeros
  wt <- diag(w)
  wtinv <- diag(1/w)

  fo_y <- as.formula(paste0("~", ynm))
  fo_null <- formula(paste(ynm, "~", 1))
  fo_by <- if (length(xnms_add) == 1) formula(paste("~", xnms_add)) else formula(paste("~", paste(xnms_add, collapse = "+")))
  ans.unc_cp_0 <- NULL

  if (family$family == "gaussian") {
    ans.unc <- suppressWarnings(svyby(formula = fo_y, by = fo_by, design = ds, FUN = svymean, keep.var = TRUE,
                                      keep.names = TRUE, verbose = FALSE, vartype = "se",
                                      drop.empty.groups = FALSE, covmat = TRUE,
                                      na.rm.by = FALSE, na.rm.all = FALSE,...))
    v1 <- vcov(ans.unc)
    etahatu <- yvecu <- yvec <- ans.unc[[ynm]]

    #check more...
    #predictors in fo_by must be factors; otherwise, na will be removed
    if(length(zeros_ps) > 0) {
      v1 <- v1[-zeros_ps, -zeros_ps]
      etahatu <- etahatu[-zeros_ps]
      yvecu <- yvecu[-zeros_ps]
      yvec <- yvec[-zeros_ps]
    }

    vsu <- round(diag(v1), 5)
    ans.unc_null <- svyglm(formula = fo_null, design = ds, family = family)
    null.deviance <- ans.unc_null$deviance
    prior.weights <- ans.unc_null$prior.weights
    #new:
    ans.unc_cp_0 <- ans.unc
  } else {
    ans.unc_cp_0 <- suppressWarnings(svyby(formula = fo_y, by = fo_by, design = ds, FUN = svymean, keep.var = TRUE,
                                           keep.names = TRUE, verbose = FALSE, vartype = "se",
                                           drop.empty.groups = FALSE, covmat = TRUE,
                                           na.rm.by = FALSE, na.rm.all = FALSE,...))
    null.deviance <- ans.unc$null.deviance
    prior.weights <- ans.unc$prior.weights
    #newdata <- grid2
    #if Nd == 1, then grid2 already deletes empty cells
    #if (Nd > 1 & ne  > 0) {
    #  newdata <- newdata[-empty_cells, , drop = FALSE]
    #}
    #----------------------------------------------------------------------------------
    #use predict to get unconstrained estimate and variance and ci
    #----------------------------------------------------------------------------------
    newdata <- grid2
    #if Nd == 1, then grid2 already deletes empty cells
    if (Nd > 1 & ne  > 0) {
      newdata <- newdata[-empty_cells, , drop = FALSE]
    }
    #if Nd > 1 and ne > 0, predict won't work -> separate ne == 0 and ne > 0 as below
    if (ne == 0 | (Nd == 1 & ne > 0)) {
      p.unc <- predict(ans.unc, newdata, se.fit = TRUE, type = "link", vcov = TRUE)
      v1 <- vcov(p.unc)
      #when some nd == 1, as.data.frame will give the warning:In sqrt(diag(as.matrix(vv))) : NaNs produced
      #SE for nd == 1 will be NA: won't affect yvec or vsu
      etahatu <- yvecu <- yvec <- suppressWarnings(as.data.frame(p.unc)$link)
      vsu <- round(diag(v1), 5)
    }
    #if (ne > 0) {
    if (Nd > 1 & ne > 0) {
      # if (Nd > 1 & length(empty_cells) > 0) {
      # newdata <- newdata[-empty_cells, , drop = FALSE]
      # }
      tt <- delete.response(terms(formula(ans.unc)))
      mf <- model.frame(tt, data = newdata)
      mm <- model.matrix(tt, mf)
      # modified from predict.svrepglm
      # check more...when predictors are not coded as factor, then mm will delete the empty cells; otherwise not
      mm <- mm[, -empty_cells, drop = FALSE]
      # print (dim(mm))
      # print (dim(vcov(ans.unc)))
      vv <- mm %*% vcov(ans.unc) %*% t(mm)
      v1 <- vv
      #etahatu <- yvecu <- yvec <- NULL
      etahatu <- drop(mm %*% coef(ans.unc))
      yvec <- yvecu <- etahatu
      vsu <- round(diag(v1), 5)
    }
  }

  muhat.un <- c(family$linkinv(etahatu))
  z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
  hlu <- z.mult * sqrt(vsu)
  lwru <- yvecu - hlu
  uppu <- yvecu + hlu

  ps <- which(round(diag(v1), 7) == 0)
  #check: impute zero variance with average of non-zero variance
  if (length(ps) > 0) {
    #diag(v1)[ps] <- mean(diag(v1)[-ps])
    diag(v1)[ps] <- 0
  }

  #set hkeep = NULL for now
  #----------------------------------------------------------------------------------
  #create amat by .makeamat
  #----------------------------------------------------------------------------------
  hkeep <- NULL
  noord_cells <- NULL
  bool_ave <- map_dbl(block.ave.lst, .f = function(.x) all(.x == -1))
  bool_ord <- map_dbl(block.ord.lst, .f = function(.x) all(.x == -1))
  if (any(bool_ave == 0)) {
    noord_cells <- map(block.ave.lst, .f = function(.x) which(.x == 0))
    if (length(Ds) == 1) {
      noord_cells <- noord_cells[[1]]
    }
  }
  if (any(bool_ord == 0)) {
    noord_cells <- map(block.ord.lst, .f = function(.x) which(.x == 0))
    if (length(Ds) == 1) {
      noord_cells <- noord_cells[[1]]
    }
  }

  ans_amat <- .MakeAmat(sh, Nd, block.ave.lst, block.ord.lst, zeros_ps=empty_cells, M,
                        M_0, noord_cells, xvs2, grid2, Ds)
  amat <- ans_amat$amat
  amat_0 <- ans_amat$amat_0
  #print (dim(amat))
  #----------------------------------------------------------------------------------
  #project unconstrained estimate to get constrained estimate
  #----------------------------------------------------------------------------------
  ans.polar <- coneA(yvec, amat, w = w, msg = FALSE)
  etahat <- ans.polar$thetahat
  #etahat <- round(ans.polar$thetahat, 10)
  #print (etahat)
  face <- ans.polar$face
  # new: to be used in summary and anova
  edf <- 1.5 * ans.polar$df
  if (length(face) == 0) {
    mat <- wt %*% v1
  } else {
    dd <- amat[face, , drop = FALSE]
    wtinvdd <- wtinv %*% t(dd)
    pmat_is <- wtinvdd %*% solve(dd %*% wtinvdd) %*% dd
    mat <- wt %*% (imat - pmat_is) %*% v1
  }
  evec <- c(family$linkinv(yvec)) - c(family$linkinv(etahat))
  CIC <- t(evec) %*% wt %*% evec + 2 * (sum(diag(mat)))
  CIC <- as.numeric(CIC)
  # new: add CIC.un
  CIC.un <- 2 * (sum(diag(wt %*% v1)))
  CIC.un <- as.numeric(CIC.un)

  #----------------------------------------------------------------------------------
  #mixture variance-covariance matrix
  #----------------------------------------------------------------------------------

  acov <- v1
  # add a warning for n.mix
  if (n.mix %% 1 != 0 | n.mix <= 0) {
    message("n.mix must be a positive integer for simulation of constrained covariance! The unconstrained covariance will be used! \n")
    z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
    vsc <- diag(v1)
    hl <- z.mult * sqrt(vsc)
    lwr <- etahat - hl
    upp <- etahat + hl
  } else {
    # get the constrained variance-covariance matrix
    lwr <- NULL
    upp <- NULL
    #acov <- v1
    dp <- -t(amat)
    dp <- apply(dp, 2, function(e) e / (sum(e^2))^(.5))
    m_acc <- M
    sector <- NULL
    times <- NULL
    df.face <- NULL
    pmat_is <- NULL
    iter <- 1
    obs <- 1:M
    # non-gaussian can use mvrnorm by CLT
    ysims <- MASS::mvrnorm(n.mix, mu = etahat, Sigma = v1)
    for (iloop in 1:n.mix) {
      ysim <- ysims[iloop, ]
      ansi <- coneA(ysim, amat, w = w, msg = FALSE)
      etahati <- round(ansi$thetahat, 10)
      facei <- ansi$face
      if (length(facei) == 0) {
        next
      } else {
        sec <- 1:nrow(amat) * 0
        sec[facei] <- 1
        r <- .makebin(sec) + 1
        if (iter == 1) {
          df.face <- rbind(df.face, c(r, 1))
          sector <- rbind(sector, sec)
        } else {
          if (r %in% df.face[, 1]) {
            ps <- which(df.face[, 1] %in% r)
            df.face[ps, 2] <- df.face[ps, 2] + 1
          } else {
            df.face <- rbind(df.face, c(r, 1))
            sector <- rbind(sector, sec)
          }
        }
        iter <- iter + 1
      }
    }

    if (!is.null(df.face)) {
      imat <- diag(M)
      sm_id <- which((df.face[, 2] / n.mix) < 1e-3)
      if (any(sm_id)) {
        df.face <- df.face[-sm_id, , drop = FALSE]
        sector <- sector[-sm_id, , drop = FALSE]
      }
      nsec <- nrow(df.face)
      bsec <- df.face
      bsec[, 2] <- bsec[, 2] / sum(bsec[, 2])

      acov <- matrix(0, nrow = M, ncol = M)
      wtinv <- diag(1 / w)
      for (is in 1:nsec) {
        jvec <- sector[is, ]
        smat <- dp[, which(jvec == 1), drop = FALSE]
        wtinvs <- wtinv %*% smat
        pmat_is_p <- wtinv %*% smat %*% solve(t(smat) %*% wtinv %*% smat) %*% t(smat)
        pmat_is <- (imat - pmat_is_p)
        acov <- acov + bsec[is, 2] * pmat_is %*% v1 %*% t(pmat_is)
      }
    } else {
      pmat_is <- imat
      acov <- v1
    }
    z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
    vsc <- diag(acov)
    hl <- z.mult * sqrt(vsc)
    lwr <- etahat - hl
    upp <- etahat + hl
  }

  #----------------------------------------------------------------------------------
  #impute
  #----------------------------------------------------------------------------------
  #new
  acov_cp <- acov
  #print (acov_cp)
  if (ne >= 1) {
    # new: if there's empty cells, just augment muhatu to include NA's
    etahatu_all <- 1:(M + zeros) * 0
    etahatu_all[zeros_ps] <- NA
    etahatu_all[obs_cells] <- etahatu
    etahatu <- etahatu_all

    lwru_all <- 1:(M + zeros) * 0
    lwru_all[zeros_ps] <- NA
    lwru_all[obs_cells] <- lwru
    lwru <- lwru_all

    uppu_all <- 1:(M + zeros) * 0
    uppu_all[zeros_ps] <- NA
    uppu_all[obs_cells] <- uppu
    uppu <- uppu_all

    # for monotone only: to find the c.i. for empty/1 cells by using smallest L and largest U of neighbors
    # new: create grid2 again because the routines to do imputation need to have each x start from 1
    xvs2 <- apply(xm.red, 2, function(x) 1:length(unique(x)))
    if (is.matrix(xvs2)) {
      grid3 <- expand.grid(as.data.frame(xvs2))
    } else if (is.list(xvs2)) {
      grid3 <- expand.grid(xvs2)
    }
    # revise later: replace lwru with lwr etc
    ans_im <- .impute_em2(empty_cells, obs_cells, M = (M + ne), etahatu, Nd, nd, Nhat, w,
                          domain.ord, grid3, Ds, sh, etahat, lwr = lwr, upp = upp, vsc = vsc, acov=acov, amat_0)
    etahat <- ans_im$muhat
    lwr <- ans_im$lwr
    upp <- ans_im$upp
    vsc <- ans_im$vsc
    domain.ord <- ans_im$domain.ord
    acov_cp <- ans_im$acov_cp
  }

  sig1 <- vsc
  sig2 <- NULL
  if (Nd >= 1 & nsm >= 1) {
    new_obs_cells <- sort(unique(c(empty_cells, obs_cells)))
    if (Nd == 1) {
      ans_im2 <- .impute_em2(small_cells, new_obs_cells,
                             M = (M + ne), etahatu,
                             Nd, nd, Nhat, w, domain.ord, grid2, Ds, sh, etahat, lwr, upp,
                             vsc, acov=NULL, amat_0)
      sig2 <- ans_im2$vsc
    }
    if (Nd >= 2) {
      ans_im2 <- .impute_em3(small_cells, M = (M + ne), etahatu, Nd, nd, Nhat, w, domain.ord,
                             grid2, Ds, sh, etahat, lwr, upp, vsc, new_obs_cells, amat_0)
      sig2 <- ans_im2$vsc
    }
  }

  vsc_mix <- sig1
  vsc_mix_imp <- NULL
  if (Nd >= 1 & nsm >= 1) {
    imp_ps <- sort(c(zeros_ps, small_cells))
    nd_imp_ps <- nd[imp_ps]
    l_imp_ps <- length(nd_imp_ps)

    sig1_imp <- sig1[imp_ps]
    sig2_imp <- sig2[imp_ps]

    vsc_mix_imp <- 1:l_imp_ps * 0
    # tmp:
    if (Nd > 1) {
      nbors_lst_0 <- .fn_nbors_2(small_cells, amat_0, etahat)
    }

    k_sm <- 1
    for (k in 1:l_imp_ps) {
      ndek <- nd_imp_ps[k]
      s1k <- sig1_imp[k]
      s2k <- sig2_imp[k]
      # tmp
      if (Nd > 1 & ndek > 0) {
        max_ps <- nbors_lst_0$nb_lst_maxs[[k_sm]]
        min_ps <- nbors_lst_0$nb_lst_mins[[k_sm]]
        # max_ps = min_ps = NULL
        if (length(max_ps) > 0) {
          mu_max <- etahat[max_ps]
          # mu_max = muhat[max_ps]
        }
        if (length(min_ps) > 0) {
          mu_min <- etahat[min_ps]
          # mu_min = muhat[min_ps]
        }
        k_sm <- k_sm + 1
      } else if (Nd == 1 | ndek == 0) {
        min_ps <- max_ps <- 1
      }

      # when min_ps or max_ps = NULL: kth cell is in a 2D case, small cell, at the edge of the
      # grid and there is only one constraint in amat
      if (length(min_ps) == 0 | length(max_ps) == 0) {
        varek <- s1k
      } else {
        if (ndek == 0) {
          # varek = s1k
          varek <- s2k
        } else if (ndek == 1) {
          varek <- s1k * .2 + s2k * .8
        } else if (ndek == 2 | ndek == 3) {
          varek <- s1k * .3 + s2k * .7
        } else if (ndek == 4 | ndek == 5)  {
          varek <- s1k * .6 + s2k * .4
        } else if (ndek == 6 | ndek == 7) {
          varek <- s1k * .6 + s2k * .4
        } else if (ndek == 8 | ndek == 9 | ndek == 10) {
          varek <- s1k * .8 + s2k * .2
        } # else if(ndek == 10){
        #  varek = s1k*.8 + s2k*.2
        # }
      }
      vsc_mix_imp[k] <- varek
    }
    vsc_mix[imp_ps] <- vsc_mix_imp
  }

  hl <- z.mult * sqrt(vsc_mix)
  lwr <- etahat - hl
  upp <- etahat + hl
  # new: test of H_0: theta in V and H_1: theta in C
  # do lsq fit with Atheta = 0 to get fit under H0
  # use the 1st option to compute T1 and T2

  #----------------------------------------------------------------------------------
  #one-sided test
  #----------------------------------------------------------------------------------

  bval <- NULL
  pval <- NULL
  if (test) {
    m <- qr(amat)$rank
    vec <- amat %*% yvecu
    T1_hat <- t(vec) %*% ginv(amat %*% v1 %*% t(amat)) %*% vec
    eans <- eigen(v1, symmetric = TRUE)
    evecs <- eans$vectors
    evecs <- round(evecs, 8)
    evals <- eans$values
    # new: change negative eigenvalues to be a small positive value
    sm <- 1e-7
    neg_ps <- which(evals < sm)
    if(length(neg_ps) > 0){
      evals[neg_ps] <- sm
    }

    L <- evecs %*% diag(sqrt(evals))
    Linv <- solve(L) # solve will give singular error message when n is small
    #Linv <- ginv(L)
    atil <- amat %*% L
    Z_s <- Linv %*% yvecu

    theta_hat <- coneA(Z_s, atil, msg = FALSE)$thetahat
    T2_hat <- t(Z_s - theta_hat) %*% (Z_s - theta_hat)

    bval <- (T1_hat - T2_hat) / T1_hat

    if (bval > sm) {
      nloop <- 100
      zerovec <- rep(0, M)
      zsims <- MASS::mvrnorm(nloop, mu = zerovec, Sigma = imat)

      J_length <- rep(0, nloop)
      for (i in 1:nloop) {
        J_length[i] <- length(coneA(zsims[i, ], atil, msg = FALSE)$face)
      }

      mdist <- 0:m * 0
      for (i in 1:(m + 1)) {
        mdist[i] <- length(which(J_length == (i - 1)))
      }
      mdist <- mdist / nloop
      pval <- 1 - sum(pbeta(bval, m:0 / 2, 0:m / 2) * mdist)
    } else {
      pval <- 1
    }
  }

  #----------------------------------------------------------------------------------
  #return results
  #----------------------------------------------------------------------------------
  # set h_grid = NULL for now
  h_grid <- NULL
  muhat <- c(family$linkinv(etahat))
  muhat.un <- c(family$linkinv(etahatu))
  muhat_all <- 1:n * 0
  etahat_all <- 1:n * 0
  for (udi in uds) {
    udi <- as.numeric(udi)
    ps <- which(ID %in% udi)
    # should use muhat, not etahat
    muhat_all[ps] <- muhat[which(uds %in% udi)]
    etahat_all[ps] <- etahat[which(uds %in% udi)]
    # ybar_all[ps] = ybar_vec[which(uds%in%udi)]
  }
  deviance <- sum(family$dev.resids(y, muhat_all, prior.weights))
  logLike_CIC <- deviance + 2 * (sum(diag(mat)))
  #new: inherit attributes if svyby is used
  if (family$family == "gaussian") {
    ans.unc_cp <- data.frame(matrix(nrow = M_0, ncol = ncol(ans.unc)))
    colnames(ans.unc_cp) <- colnames(ans.unc)
    #print (ans.unc_cp)
    #ans.unc_cp$ID <- factor(1:M_0)
    #print (which(colnames(ans.unc_cp) %in% xnms_add))
    ans.unc_cp <- ans.unc_cp_0
    #ans.unc_cp[,which(colnames(ans.unc_cp) %in% xnms_add)] <- ans.unc_cp_0[,which(colnames(ans.unc_cp) %in% xnms_add)]
    ans.unc_cp[[ynm]] <- as.numeric(etahat)
    ans.unc_cp$se <- vsc_mix^(.5)
    #print (ans.unc_cp)
    #print (ans.unc_cp_0)
    attr(ans.unc_cp, "class") <- c("svyby", "data.frame")
    attr(ans.unc_cp, "svyby") <- attr(ans.unc, "svyby")
    attr(ans.unc_cp, "var") <- acov_cp
  } else {
    ans.unc_cp <- data.frame(matrix(nrow = M_0, ncol = ncol(ans.unc_cp_0)))
    colnames(ans.unc_cp) <- colnames(ans.unc_cp_0)
    #ans.unc_cp$ID <- factor(1:M_0)
    ans.unc_cp <- ans.unc_cp_0
    #ans.unc_cp[,which(colnames(ans.unc_cp) %in% xnms_add)] <- ans.unc_cp_0[,which(colnames(ans.unc_cp) %in% xnms_add)]
    ans.unc_cp[[ynm]] <- as.numeric(etahat)
    ans.unc_cp$se <- vsc_mix^(.5)
    attr(ans.unc_cp, "class") <- c("svyby", "data.frame")
    attr(ans.unc_cp, "svyby") <- attr(ans.unc_cp_0, "svyby")
    attr(ans.unc_cp, "var") <- acov_cp
  }
  if ((n - edf) <= 0) {
    df.residual <- edf
  } else {
    df.residual <- n - edf
  }
  #set pskeep = NULL for now
  pskeep <- NULL
  ans <- list(
    linear.predictors = etahat_all, etahat = c(etahat), linear.predictors.un = ans.unc$linear.predictors,
    etahatu = c(yvecu), fitted.values = muhat_all, muhat = muhat, fitted.values.un = ans.unc$fitted.values,
    muhatu = muhat.un, cov.un = v1, # unconstrained cov
    cov.unscaled = vsc_mix, n.mix = n.mix,# constrained cov with imputation; diagonal only
    null.deviance = null.deviance, df.null = n - 1, deviance = deviance, df.residual = df.residual,
    edf = edf, rank = edf, lwr = c(lwr), upp = c(upp), lwru = c(lwru), uppu = c(uppu), y = y,
    domain = as.numeric(levels(ID))[ID], acov = acov, acov_cp = acov_cp,
    grid_ps = pskeep, xvs2 = xvs2, uds = uds,
    amat = amat, Ds = Ds, domain.ord = domain.ord, nd = nd, grid = grid2, ec = empty_cells,
    hkeep = hkeep, h_grid = h_grid, zeros_ps = zeros_ps, pval = pval, bval = bval,
    CIC = CIC, CIC.un = CIC.un, logLike_CIC = logLike_CIC, w = w, xm.red = xm.red, xvs2 = xvs2,
    ne = ne, empty_cells = empty_cells, small_cells = small_cells, obs_cells = obs_cells,
    zeros = zeros, Nd = Nd, Nhat = Nhat, amat_0 = amat_0, weights = w,
    prior.weights = prior.weights, data = df, #terms = ans.unc$terms,
    offset = ans.unc$offset, contrasts = ans.unc$contrasts, survey.design = design, ans.unc = ans.unc, ans.unc_cp = ans.unc_cp
  )
  return(ans)
}

#---------------------------------------------------------------------------------------------------------------------------------------
# csvy.control <- function(epsilon = 1e-05, maxit = 10, deff = FALSE, multicore = FALSE){
#   if (!is.numeric(epsilon) | epsilon <= 0) {
#     stop("value of 'epsilon' must be > 0")
#   }
#   if (!is.numeric(maxit) | maxit <= 0) {
#     stop("maximum number of iterations must be > 0")
#   }
#   list(epsilon = epsilon, maxit = maxit, deff = deff, multicore)
# }


#############################################
#block ordering and block average
#other shape functions are from cgam
#############################################
block.Ord <- function(x, order = NULL, numknots = 0, knots = 0, space = "E") {
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 9
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "block.ord"
  attr(x, "order") <- order
  # class(x) <- "additive"
  x
}

# block.Ave <- function(x, order = NULL, numknots = 0, knots = 0, space = "E") {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 10
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "block.ave"
#   attr(x, "order") <- order
#   # class(x) <- "additive"
#   x
# }


########################################################
# confint.csvy
# inherit from svyby?
# use type = c("link", "response") as ...
# set parm = NULL
########################################################
confint.csvy <- function(object, parm, level = 0.95, type = c("link", "response"),...) {
  type <- match.arg(type)
  n.mix <- object$n.mix
  #if (n.mix %% 1 != 0 | n.mix <= 0){
  #  stop("confint function only works when n.mix is a positive integer!")
  #} else {
    z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
    # vs <- vcov(object)
    vs <- object$cov.unscaled
    lwr <- object$etahat - z.mult * sqrt(vs)
    upp <- object$etahat + z.mult * sqrt(vs)

    lwr <- switch(type,
                  link = lwr,
                  response = object$family$linkinv(lwr)
    )
    upp <- switch(type,
                  link = upp,
                  response = object$family$linkinv(upp)
    )

    ci.bands <- as.data.frame(do.call(cbind, list(lwr, upp)))
    # learn from confint.svyglm:
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    pct <- paste0(format(100 * a,
                         trim = TRUE,
                         scientific = FALSE, digits = 3
    ), "%")
    colnames(ci.bands) <- pct
    return(ci.bands)
  #}
}


#############################################
# summary for csvy
# inherit from cgam, use NextMethod()?
#############################################
summary.csvy <- function(object, ...) {
  family <- object$family
  CIC <- object$CIC
  CIC.un <- object$CIC.un
  deviance <- object$deviance
  null.deviance <- object$null.deviance
  bval <- object$bval
  pval <- object$pval
  # edf
  df.residual <- object$df.residual
  df.null <- object$df.null
  edf <- object$edf
  tms <- object$terms
  rslt <- NULL
  # learn from summary.svyglm
  # not work for gaussian; svyby
  #dispersion <- svyvar(resid(object,"pearson"), object$survey.design, na.rm = TRUE)
  #disperson <- survey:::summary.svrepglm(object)$disperson
  if(inherits(object, "svrepglm")){
    presid <- resid(object,"pearson")
    dispersion <- sum(object$survey.design$pweights*presid^2, na.rm = TRUE) / sum(object$survey.design$pweights)
  } else if (inherits(object, "svyglm")){
    dispersion <- svyvar(resid(object, "pearson"), object$survey.design, na.rm = TRUE)
  }
  if (!is.null(pval)) {
    # rslt = data.frame("edf" = round(resid_df_obs, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
    # same as cgam:
    rslt <- data.frame("edf" = round(edf, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
    # check?
    rownames(rslt) <- rev((attributes(tms)$term.labels))[1]
    structure(list(
      call = object$call, coefficients = rslt, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
      null.deviance = null.deviance, deviance = deviance, family = family, dispersion = dispersion
    ), class = "summary.csvy")
  } else {
    structure(list(
      call = object$call, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
      null.deviance = null.deviance, deviance = deviance, family = family, dispersion = dispersion
    ), class = "summary.csvy")
    # warning("summary function for a csvy object only works when test = TRUE!")
  }
}


##############################
# print.summary.csvy
# inherit from cgam??
##############################
print.summary.csvy <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  # cat(x$call, "\n")
  cat("Null deviance: ", round(x$null.deviance, 4), "", "on", x$df.null, "", "degrees of freedom", "\n")
  cat("Residual deviance: ", round(x$deviance, 4), "", "on", x$df.residual, "", "observed degrees of freedom", "\n")
  if (!is.null(x$coefficients)) {
    # cat("\n")
    cat("Approximate significance of constrained fit: \n")
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  }
  # if (!is.null(x$CIC)) {
  cat("CIC (constrained estimator): ", round(x$CIC, 4), "\n")
  cat("CIC (unconstrained estimator): ", round(x$CIC.un, 4), "\n")
  # }
  # cat("(Dispersion parameter for",  x$family$family, "taken to be", round(x$dispersion, 4),  ")", "\n")
}


#---------------------------------------------------------------------------------------------------------------------------------------
#use fitted to return etahat
#add an ans.unc_cp to each family option
coef.csvyby <- function(object,...) {
  x <- object$ans.unc_cp
  coef(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
fitted.csvy <- function(object, ...) {
  ans <- object$muhat
  ans
}


##########################################################
#3d persp plot when the object has at least two predictors
##########################################################
plotpersp.csvy <- function(object, x1 = NULL, x2 = NULL, x1nm = NULL, x2nm = NULL, data = NULL, ci = "none",
                           transpose = FALSE, main = NULL, categ = NULL, categnm = NULL,
                           surface = c("C", "U"), type = c("link", "response"),
                           col = "white", cex.main = .8, xlab = NULL,
                           ylab = NULL, zlab = NULL, zlim = NULL, box = TRUE,
                           axes = TRUE, th = NULL, ltheta = NULL,
                           ticktype = "detailed", nticks = 5, palette = NULL, NCOL = NULL, ...) {
  x1nm <- deparse(substitute(x1))
  x2nm <- deparse(substitute(x2))
  type <- match.arg(type)
  surface <- match.arg(surface)
  if (!inherits(object, "csvy")) {
    warning("calling plotpersp(<fake-csvy-object>) ... \n")
  }
  t_col <- function(color, percent = 80, name = NULL) {
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 maxColorValue = 255,
                 alpha = (100 - percent) * 255 / 100, names = name
    )
    ## Save the color
    invisible(t.col)
  }
  col.upp <- t_col("green", perc = 90, name = "lt.green")
  col.lwr <- t_col("pink", perc = 80, name = "lt.pink")
  Ds <- object$Ds
  # muhat = object$muhat.s
  # lwr = object$lwr.s
  # upp = object$upp.s
  xnms <- object$xnms_add
  xmat <- object$xmat_add
  #bool <- apply(xmat, 2, function(x) is.numeric(x))
  #if (any(!bool)) {
  #  xmat <- apply(xmat, 2, function(x) as.numeric(x))
  #}
  # xmat0 <- object$xmat_add0
  xmat0 <- object$xmat_add
  #bool <- apply(xmat0, 2, function(x) is.numeric(x))
  #if (any(!bool)) {
  #  xmat0 <- apply(xmat0, 2, function(x) as.numeric(x))
  #}
  shp <- object$shapes_add
  ynm <- object$ynm
  dm <- object$domain
  #new: to avoid empty cell
  #dm = as.numeric(levels(dm))[dm]
  nd <- object$nd
  ec <- object$ec
  xvs2 <- object$xvs2
  uds <- object$uds
  data <- object$data
  # new: for more families
  family <- object$family
  # muhat = object$muhat
  etahat <- object$etahat
  lwr <- object$lwr # for eta
  upp <- object$upp # for eta
  main <- ifelse(is.null(main), "Constrained Fit", main)

  switch(surface,
         U = {
           etahat <- object$etahatu
           lwr <- object$lwru
           upp <- object$uppu
           main <- "Unconstrained Fit"
         },
         C =
  )

  switch(type,
         response = {
           etahat <- family$linkinv(etahat)
           lwr <- family$linkinv(lwr)
           upp <- family$linkinv(upp)
         },
         link =
  )

  knms <- length(xnms)
  obs <- 1:knms
  if (x1nm == "NULL" | x2nm == "NULL") {
    if (length(xnms) >= 2) {
      if (is.null(categ)) {
        x1nm <- xnms[1]
        x2nm <- xnms[2]
        x1id <- 1
        x2id <- 2
      } else {
        if (!is.character(categ)) {
          stop("categ must be a character argument! \n")
        } else if (any(grepl(categ, xnms))) {
          id <- which(grepl(categ, xnms))
          x12nms <- xnms[-id]
          x12ids <- obs[-id]
          x1nm <- x12nms[1]
          x2nm <- x12nms[2]
          x1id <- x12ids[1]
          x2id <- x12ids[2]
        }
      }
    } else {
      stop("Number of non-parametric predictors must >= 2! \n")
    }
  } else {
    # new: use the data as the environment for x1 and x2
    x1 <- data[[x1nm]]
    x2 <- data[[x2nm]]
    if (all(xnms != x1nm)) {
      if (length(x1) != nrow(xmat)) {
        stop("Number of observations in the data set is not the same as the number of elements in x1! \n")
      }
      bool <- apply(xmat, 2, function(x) all(x1 == x))
      if (any(bool)) {
        x1id <- obs[bool]
        # change x1nm to be the one in formula
        x1nm <- xnms[bool]
      } else {
        stop(paste(paste("'", x1nm, "'", sep = ""), "is not a predictor defined in the model!", "\n"))
      }
    } else {
      x1id <- obs[xnms == x1nm]
    }

    if (all(xnms != x2nm)) {
      if (length(x2) != nrow(xmat)) {
        stop("Number of observations in the data set is not the same as the number of elements in x2! \n")
      }
      bool <- apply(xmat, 2, function(x) all(x2 == x))
      if (any(bool)) {
        x2id <- obs[bool]
        x2nm <- xnms[bool]
      } else {
        stop(paste(paste("'", x2nm, "'", sep = ""), "is not a predictor defined in the model! \n"))
      }
    } else {
      x2id <- obs[xnms == x2nm]
    }
  }
  # xmat keeps the order of x's in the original data frame, but x1 and x2 may change depending on x1id and x2id
  # new: not use
  x1 <- xmat[, x1id]
  x2 <- xmat[, x2id]

  if (is.null(xlab)) {
    # xlab = deparse(x1nm)
    xlab <- x1nm
  }
  if (is.null(ylab)) {
    # ylab = deparse(x2nm)
    ylab <- x2nm
  }
  if (is.null(zlab)) {
    if ("response" %in% type) {
      zlab <- paste("Est mean of", ynm)
    }
    if ("link" %in% type) {
      zlab <- paste("Est systematic component of", ynm)
    }
  }

  D1 <- Ds[x1id]
  D2 <- Ds[x2id]

  #grid2 <- cbind(x1, x2, dm)
  #ord <- order(dm)
  # order x1 and x2 in ascending order as how we estimate muhat
  #grid2 <- grid2[ord, , drop = FALSE]
  #dm <- sort(dm)

  #new
  grid2 <- object$grid

  #new: augment dm to include zeros_ps
  # zeros_ps <- object$zeros_ps
  # dm2 <- 1:length(object$nd)
  # dm2[which(object$nd %in% table(dm))] <- dm

  ps <- NULL
  ps_id <- 1:nrow(grid2)
  knms <- length(xnms)
  obs <- 1:knms
  if (knms >= 3) {
    if (is.null(categ)) {
      x3id <- obs[-c(x1id, x2id)]
      kx3 <- length(x3id)
      # for (i in 1:kx3) {
      #   x3i <- xmat[, x3id[i]]
      #   x3i <- x3i[ord]
      #   x3i0 <- xmat0[, x3id[i]]
      #   x3i0 <- x3i0[ord]
      #   x3i_use <- max(x3i[x3i <= median(x3i)])
      #   ps <- c(ps, x3i_use)
      #   ps_idi <- which(x3i == x3i_use)
      #   ps_id <- intersect(ps_id, ps_idi)
      # }
      # domain_id <- unique(dm[ps_id])
      # etahat_use <- etahat[domain_id]

      #new: to handle zero_ps
      for (i in 1:kx3) {
        x3i <- grid2[, x3id[i]]
        if(is.factor(x3i)){
          #new: choose the middle level if x3i is a factor
          nlvs <- length(levels(x3i))
          ps_lvl <- round(nlvs/2) + 1
          x3i_use <- levels(x3i)[ps_lvl]
        } else {
          #choose the median value if x3i is numeric
          x3i_use <- max(x3i[x3i <= median(x3i)])
        }
        #ps_idi <- which(x3i == x3i_use)
        ps_idi <- which(x3i %in% x3i_use)
        ps_id <- intersect(ps_id, ps_idi)
      }
      domain_id <- ps_id
      etahat_use <- etahat[domain_id]
      upp_use <- upp[domain_id]
      lwr_use <- lwr[domain_id]
    } else {
      #need to test more
      x3nms <- xnms[-c(x1id, x2id)]
      x3mat <- xmat[, -c(x1id, x2id), drop = FALSE]
      x3mat0 <- xmat0[, -c(x1id, x2id), drop = FALSE]
      # print (head(x3mat))
      if (!is.character(categ)) {
        stop("categ must be a character argument! \n")
      } else if (any(grepl(categ, x3nms))) {
        id <- which(grepl(categ, x3nms))
        x3nmi <- x3nms[id]
        if (x3nmi %in% c(x1nm, x2nm)) {
          stop("categ must be different than x1 and x2! \n")
        }
        x3i <- x3mat[, id]
        #x3i <- x3i[ord]
        x3i0 <- x3mat0[, id]
        #x3i0 <- x3i0[ord]

        ps_id4 <- NULL
        #ps_id is the same as 1:M
        ps_id <- 1:nrow(grid2)
        if (ncol(x3mat) > 1) {
          x4s <- x3mat[, -id, drop = FALSE]
          #x4s_use <- apply(x4s, 2, min) #x4s may not be integers
          #new:
          x4s_use <- NULL
          for(i in 1:ncol(x4s)){
            x4i <- x4s[, i]
            if (is.numeric(x4i) | is.integer(x4i)) {
              x4s_use <- c(x4s_use, min(x4i))
            } else if (is.factor(x4i)) {
              nlvs <- length(levels(x4i))
              ps_lvl <- 1
              x4i_use <- levels(x4i)[ps_lvl]
              x4s_use <- c(x4s_use, x4i_use)
            }
          }

          if(is.factor(x3i)){
            #new: choose the middle level if x3i is a factor
            nlvs <- length(levels(x3i))
            ps_lvl <- round(nlvs/2) + 1
            x3i_use <- levels(x3i)[ps_lvl]
          } else {
            #choose the median value if x3i is numeric
            x3i_use <- max(x3i[x3i <= median(x3i)])
          }

          #new:
          x4nms <- x3nms[-id]
          for (i in 1:ncol(x4s)) {
            # print (i)
            # print (x4s)
            #x4i <- x4s[, i, drop = FALSE]
            #x4i <- x4i[ord]
            #x4i_use <- x4s_use[i]
            #ps_id4 <- which(x4i == x4i_use)
            #ps_id <- intersect(ps_id, ps_id4)
            x4i_use <- x4s_use[i]
            x4i_col_id <- which(grepl(x4nms[i], object$xnms_add))
            ps_id4 <- which(grid2[,  x4i_col_id] %in% x4i_use)
            ps_id <- intersect(ps_id, ps_id4)
          }
        }

        surfs <- list()
        dms <- list()

        #new: to be used in grid2
        x3_col_id <- which(grepl(categ, object$xnms_add))

        #ux3i <- unique(x3i)
        #ux3i0 <- unique(x3i0)
        ux3i <- xvs2[[x3_col_id]]
        ux3i0 <- ux3i

        kz_add <- length(ux3i)
        mins <- maxs <- NULL
        dm <- object$domain
        dm <- sort(dm)
        for (iz in 1:kz_add) {
          x3i_use <- ux3i[iz]
          #ps_idi <- which(x3i == x3i_use)
          #ps_idi <- which(x3i %in% x3i_use)
          #ps_id_iz <- intersect(ps_id, ps_idi)
          #domain_id <- unique(dm[ps_id_iz])
          ps_idi <- which(grid2[,  x3_col_id] %in% x3i_use)
          ps_id_iz <- intersect(ps_id, ps_idi)
          domain_id <- ps_id_iz
          dms[[iz]] <- domain_id
          etahat_use <- etahat[domain_id]
          surfs[[iz]] <- etahat_use
          mins <- c(mins, min(etahat_use))
          maxs <- c(maxs, max(etahat_use))
        }
      } else {
        # print(paste(categ, "is not an exact character name defined in the csurvey fit!"), "\n")
        cat(categ, "is not an exact character name defined in the csurvey fit!", "\n")
      }
    }
  } else {
    #domain_id <- unique(dm)
    domain_id <- object$uds
    etahat_use <- etahat
    upp_use <- upp
    lwr_use <- lwr
  }

  if (is.null(categ)) {
    if (ci != "none") {
      upp_use <- upp
      lwr_use <- lwr
      rg <- max(upp_use, na.rm = TRUE) - min(lwr_use, na.rm = TRUE)
      z.lwr <- min(etahat_use, na.rm = TRUE) - rg / 3
      z.upp <- max(etahat_use, na.rm = TRUE) + rg / 3
    } else if (ci == "none") {
      rg <- max(etahat_use, na.rm = TRUE) - min(etahat_use, na.rm = TRUE)
      z.lwr <- min(etahat_use, na.rm = TRUE) - rg / 3
      z.upp <- max(etahat_use, na.rm = TRUE) + rg / 3
    }
  } else {
    z.lwr <- min(mins, na.rm = TRUE)
    z.upp <- max(maxs, na.rm = TRUE)
  }
  if (is.null(zlim)) {
    zlim0 <- c(z.lwr, z.upp)
  } else {
    zlim0 <- zlim
  }

  ang <- NULL
  if (is.null(th)) {
    if (all(shp == 1)) {
      ang <- -40
    } else if (all(shp == 2)) {
      ang <- 40
    } else if (shp[1] == 1 & shp[2] == 2) {
      ang <- -50
    } else if (shp[1] == 2 & shp[2] == 1) {
      ang <- -230
    } else {
      ang <- -40
    }
  } else {
    ang <- th
  }
  if (is.null(categ)) {
    x1p <- 1:D1
    x2p <- 1:D2
    surf1 <- surf2 <- surf3 <- matrix(0, nrow = D1, ncol = D2)
    # new: some x won't start from 1
    #ux1 <- sort(unique(x1))
    #ux2 <- sort(unique(x2))
    ux1 <- sort(xvs2[[x1id]])
    ux2 <- sort(xvs2[[x2id]])
    if (knms == 1 | knms == 2) {
      # for (di in domain_id) {
      #   rid <- (grid2[grid2[, 3] == di, , drop = FALSE])[1, ]
      #   # new: some x won't start from 1
      #   ps1 <- which(ux1 %in% rid[1])
      #   ps2 <- which(ux2 %in% rid[2])
      #   surf1[ps1, ps2] <- etahat_use[di]
      #   # surf1[rid[1], rid[2]] = muhat_use[di]
      #   if (!is.null(upp) & !is.null(lwr)) {
      #     surf2[ps1, ps2] <- upp_use[di]
      #     surf3[ps1, ps2] <- lwr_use[di]
      #     # surf2[rid[1], rid[2]] = upp_use[di]
      #     # surf3[rid[1], rid[2]] = lwr_use[di]
      #   }
      # }
      surf1 <- matrix(etahat_use, nrow = D1, ncol = D2)
      if (!is.null(upp) & !is.null(lwr)) {
        surf2 <- matrix(upp_use, nrow = D1, ncol = D2)
        surf3 <- matrix(lwr_use, nrow = D1, ncol = D2)
      }
    } else if (knms >= 3) {
      # for (di in domain_id) {
      #   rid <- (grid2[grid2[, 3] == di, , drop = FALSE])[1, ]
      #   # new: some x won't start from 1
      #   ps1 <- which(ux1 %in% rid[1])
      #   ps2 <- which(ux2 %in% rid[2])
      #   surf1[ps1, ps2] <- etahat_use[which(domain_id %in% di)]
      #   # surf1[rid[1], rid[2]] = muhat_use[which(domain_id %in% di)]
      #   if (!is.null(upp) & !is.null(lwr)) {
      #     surf2[ps1, ps2] <- upp_use[which(domain_id %in% di)]
      #     surf3[ps1, ps2] <- lwr_use[which(domain_id %in% di)]
      #     # surf2[rid[1], rid[2]] = upp_use[which(domain_id %in% di)]
      #     # surf3[rid[1], rid[2]] = lwr_use[which(domain_id %in% di)]
      #   }
      # }
      surf1 <- matrix(etahat_use, nrow = D1, ncol = D2)
      if (!is.null(upp) & !is.null(lwr)) {
        surf2 <- matrix(upp_use, nrow = D1, ncol = D2)
        surf3 <- matrix(lwr_use, nrow = D1, ncol = D2)
      }
    }

    # empty cell: test 3D
    # if (length(ec) >= 1) {
    #   if (knms == 1 | knms == 2) {
    #     surf1[ec] <- etahat_use[ec]
    #     if (!is.null(upp) & !is.null(lwr)) {
    #       surf2[ec] <- upp_use[ec]
    #       surf3[ec] <- lwr_use[ec]
    #     }
    #   } else if (knms >= 3) {
    #     surf1[which(domain_id %in% ec)] <- etahat_use[which(domain_id %in% ec)]
    #     if (!is.null(upp) & !is.null(lwr)) {
    #       surf2[which(domain_id %in% ec)] <- upp_use[which(domain_id %in% ec)]
    #       surf3[which(domain_id %in% ec)] <- lwr_use[which(domain_id %in% ec)]
    #     }
    #   }
    # }
    # not to change users' setting
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (ci == "none") {
      persp(x1p, x2p, surf1, zlim = zlim0, col = col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang,
            ltheta = -135, ticktype = ticktype, box = box, axes = axes, nticks = nticks,...)
      par(new = FALSE)
    } else if (ci == "up") {
      persp(x1p, x2p, surf1, zlim = zlim0, col = col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang, ltheta = -135, ticktype = ticktype,
            box = box, axes = axes, nticks = nticks,...)
      par(new = TRUE)
      persp(x1p, x2p, surf2, zlim = zlim0, col = col.upp, theta = ang, box = FALSE, axes = FALSE,...)
      par(new = FALSE)
    } else if (ci == "lwr") {
      persp(x1p, x2p, surf1, zlim = zlim0, col = col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang, ltheta = -135, ticktype = ticktype,
            box = box, axes = axes, nticks = nticks,...)
      par(new = TRUE)
      persp(x1p, x2p, surf3, zlim = zlim0, col = col.lwr, theta = ang, box = FALSE, axes = FALSE,...)
      par(new = FALSE)
    }
  } else {
    # new: some x won't start from 1
    #ux1 <- sort(unique(x1))
    #ux2 <- sort(unique(x2))
    ux1 <- sort(xvs2[[x1id]])
    ux2 <- sort(xvs2[[x2id]])
    if (is.null(palette)) {
      palette <- c(
        "peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen",
        "seagreen1", "palegreen", "azure", "whitesmoke"
      )
    }
    # palette = scales::hue_pal(h = c(275, 360))(40)
    ks <- length(surfs)
    if (ks > 1 & ks < 11) {
      col <- palette[1:ks]
    } else {
      # new: use rainbow
      col <- topo.colors(ks)
    }
    width <- ks^.5
    wd <- round(width)
    if (width > wd) {
      wd <- wd + 1
    }
    if ((wd^2 - ks) >= wd) {
      fm <- c(wd, wd - 1)
    } else {
      fm <- rep(wd, 2)
    }
    if (wd > 3) {
      cex <- .7
      cex.main <- .8
    }
    if (transpose) {
      fm <- rev(fm)
    }
    # not to change users' setting
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (is.null(NCOL)) {
      par(mfrow = c(fm[1], fm[2]))
    } else {
      nr <- fm[1] * fm[2] / NCOL
      par(mfrow = c(nr, NCOL))
    }

    par(mar = c(4, 1, 1, 1))
    par(cex.main = cex.main)

    # print (ks)
    x1p <- 1:D1
    x2p <- 1:D2
    for (i in 1:ks) {
      surfi <- surfs[[i]]
      dmi <- dms[[i]]
      # print (surfi)
      #surf1 <- matrix(0, nrow = D1, ncol = D2)
      # for (j in 1:(D1 * D2)) {
      #   di <- dmi[j]
      #   rid <- (grid2[grid2[, 3] == di, , drop = FALSE])[1, ]
      #   # new: some x won't start from 1
      #   ps1 <- which(ux1 %in% rid[1])
      #   ps2 <- which(ux2 %in% rid[2])
      #   surf1[ps1, ps2] <- surfi[j]
      #   # surf1[rid[1], rid[2]] = surfi[j]
      # }
      surf1 <- matrix(surfi, nrow = D1, ncol = D2)
      if (is.null(categnm)) {
        persp(x1p, x2p, surf1,
              zlim = zlim0, col = col[i], xlab = xlab,
              ylab = ylab, zlab = zlab, main = paste0(x3nmi, " = ", ux3i0[i]),
              theta = ang, ltheta = -135, ticktype = ticktype, box = box, axes = axes, nticks = nticks,...
        )
      } else {
        if (length(categnm) == ks) {
          persp(x1p, x2p, surf1,
                zlim = zlim0, col = col[i], xlab = xlab,
                ylab = ylab, zlab = zlab, main = categnm[i],
                theta = ang, ltheta = -135, ticktype = ticktype, box = box, axes = axes, nticks = nticks,...
          )
        } else if (length(categnm) == 1) {
          persp(x1p, x2p, surf1,
                zlim = zlim0, col = col[i], xlab = xlab,
                ylab = ylab, zlab = zlab, main = paste0(categnm, " = ", ux3i0[i]),
                theta = ang, ltheta = -135, ticktype = ticktype, box = box, axes = axes, nticks = nticks,...
          )
        }
      }
      # par(new = TRUE)
    }
    par(new = FALSE)
  }
  # par(new=FALSE)
}

##############
# predict.csvy#
##############
predict.csvy <- function(object, newdata, type = c("link", "response"),
                         se.fit = FALSE, level = 0.95, ...) {
  family <- object$family
  type <- match.arg(type)
  if (!inherits(object, "csvy")) {
    warning("calling predict.csvy(<fake-csvy-object>) ... \n")
  }
  if (missing(newdata)) {
    fit <- object$etahat
    newdata <- object$xmat_add
    #print ( newdata )
    #ans <- switch(type, link = fit, response = object$family$linkinv(fit))
    #return (ans)
  } else if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame! \n")
  }
  # } else if (!all(colnames(newdata) %in% object$xnms_add) | !object$xnms_add %in% all(colnames(newdata))){
  #   stop("variables in newdata much be the same variables used in the fit!")
  # }
  #else {
  tt <- object$terms
  Terms <- delete.response(tt)
  newdata <- model.frame(Terms, newdata)
  #print (head(newdata))
  #}
  # match newdata with grid
  xnms_add <- object$xnms_add #names without symbolic functions
  #grid2 will have colnames as incr + x; same as newdata
  grid <- object$grid
  grid2 <- model.frame(Terms, grid)
  xvs2 <- object$xvs2

  #if (!is.data.frame(newdata)) {
  #  stop("newdata must be a data frame!")
  #}
  # learnt from predict.svyglm:
  # if (type=="terms")
  #  return(predterms(object,se=se.fit,...))
  # make it simpler?
  new_nmes <- colnames(newdata)
  if (ncol(newdata) == 1) {
    defined_nm <- map_lgl(new_nmes, function(elem) any(colnames(grid2) %in% elem))
    if(any(!defined_nm)){
      ps_not_def_nm <- new_nmes[!defined_nm]
      stop(cat(ps_not_def_nm, "cannot be found among predictors in the fit!", "\n"))
    } else {
      defined <- apply(newdata, 1, function(elem) any(apply(grid, 1, function(gi) all(gi == elem))))
      ps_not_def <- which(!defined)
      lvl_not_def <- unique(newdata[ps_not_def, ])
      if(all(!defined)){
        stop(paste("new factor level", paste(shQuote(lvl_not_def), collapse = ","), "cannot be found among factor levels in the fit!", sep = " ", "\n"))
      } else {
        if(any(!defined)){
          #print (as.character(lvl_not_def))
          warning(cat("new factor level:",  paste(shQuote(lvl_not_def), collapse = ","),
                      "cannot be found among factor levels in the fit!", "\n"))
        }
        #ps <- apply(newdata, 1, function(elem) grid2[which(apply(grid2, 1, function(gi) all(gi == elem))), ])
        ps <- apply(newdata, 1, function(elem) grid[apply(grid, 1, function(gi) all(gi == elem)), ])
        ps <- as.numeric(unlist(ps))
      }
    }
    #print (ps)
  }
  if (ncol(newdata) > 1) {
    #check more
    defined_nm <- map_lgl(new_nmes, function(elem) any(colnames(grid2) %in% elem))
    if(any(!defined_nm)){
      ps_not_def_nm <- new_nmes[!defined_nm]
      stop(cat(ps_not_def_nm, "cannot be found among predictors in the fit!", "\n"))
    } else {
      # order the column in case the user doesn't define the variables as the order used in the formula
      #print (newdata)
      #newdata <- newdata[, ps0]
      #print (newdata)
      newd.red <- unique(newdata)
      defined <- matrix(NA, nrow = nrow(newd.red), ncol = ncol(newd.red))
      for(ic in 1:ncol(newd.red)){
        vi <- newd.red[,ic]
        defined[,ic] <- vi %in% grid[,ic]
      }

      defined_all <- apply(defined, 1, all)
      ps_not_def <- which(!defined, arr.ind = TRUE) #return row and colunm of not defined
      if(all(!defined_all)){
        #stop(paste("factor level of", paste(new_nmes[ps_not_def[,2]], collapse = ","),
        #           "cannot be found among factor levels in the fit!", sep = " "))
        stop(paste("factor level of", paste(unique(xnms_add[ps_not_def[,2]]), collapse = ","),
                   "cannot be found among factor levels in the fit!", sep = " ", "\n"))
        #stop(paste("new factor level", paste(shQuote(lvl_not_def), collapse = ","), "cannot be found among factor levels in the fit!", sep = " "))
      } else {
        if(any(!defined)){
          #print (as.character(lvl_not_def))
          #warning(cat("new factor level:",  paste(shQuote(lvl_not_def), collapse = ","),
          #            "cannot be found among factor levels in the fit!", "\n"))
          warning(cat("factor level of", paste(unique(xnms_add[ps_not_def[,2]]), collapse = ","), "cannot be found among factor levels in the fit!", "\n"))
          ps <- apply(newdata[-ps_not_def[1],], 1, function(elem) which(apply(grid, 1, function(gi) all(gi == elem))))
        } else {
          ps <- apply(newdata, 1, function(elem) which(apply(grid, 1, function(gi) all(gi == elem))))
          ps <- as.numeric(unlist(ps))
        }
      }
    }
  }

  etahat <- object$etahat
  n.mix <- object$n.mix
  vsc_mix <- object$cov.unscaled
  # interval = match.arg(interval)
  if (!se.fit) {
    fit <- etahat[ps]
    fit <- switch(type,
                  link = fit,
                  response = object$family$linkinv(fit)
    )
    return(fit)
  } else {
    # when n.mix = 0, then there is no vsc_mix
    #if (n.mix == 0) {
    if(n.mix %% 1 != 0 | n.mix <= 0){
      message("mixture of variance-covariance matrix is not found and is being simulated...\n")
      #ynm <- object$ynm
      amat <- object$amat
      nd <- object$nd
      M <- length(nd)
      #muhat <- object$muhat
      etahat <- object$etahat
      # add in csvy
      # revise later: etahat.s
      etahatu <- object$etahatu
      w <- object$w
      v1 <- object$cov.un
      domain.ord <- object$domain.ord
      # move the imputation code to be in predict.csvy
      xm.red <- object$xm.red
      zeros <- ne <- object$ne
      xvs2 <- apply(xm.red, 2, function(x) 1:length(unique(x)))
      zeros_ps <- empty_cells <- object$empty_cells
      small_cells <- object$small_cells
      obs_cells <- object$obs_cells
      nsm <- length(small_cells)
      #zeros <- object$zeros
      Nhat <- object$Nhat
      Ds <- object$Ds
      sh <- object$shapes_add
      amat_0 <- object$amat_0
      Nd <- object$Nd

      lwr <- upp <- NULL
      acov <- v1
      dp <- -t(amat)
      dp <- apply(dp, 2, function(e) e / (sum(e^2))^(.5))

      M <- M - zeros
      m_acc <- M
      sector <- NULL
      times <- NULL
      df.face <- NULL
      iter <- 1
      obs <- 1:M
      n.mix <- 100
      if (length(empty_cells) > 0) {
        etahat <- etahat[-empty_cells]
      }
      ysims <- MASS::mvrnorm(n.mix, mu = etahat, Sigma = v1)
      for (iloop in 1:n.mix) {
        ysim <- ysims[iloop, ]
        ansi <- coneA(ysim, amat, w = w, msg = FALSE)
        etahati <- round(ansi$thetahat, 10)
        facei <- ansi$face
        if (length(facei) == 0) {
          next
        } else {
          sec <- 1:nrow(amat) * 0
          sec[facei] <- 1

          r <- .makebin(sec) + 1
          if (iter == 1) {
            df.face <- rbind(df.face, c(r, 1))
            sector <- rbind(sector, sec)
          } else {
            if (r %in% df.face[, 1]) {
              psi <- which(df.face[, 1] %in% r)
              df.face[psi, 2] <- df.face[psi, 2] + 1
            } else {
              df.face <- rbind(df.face, c(r, 1))
              sector <- rbind(sector, sec)
            }
          }
          iter <- iter + 1
        }
      }

      if (!is.null(df.face)) {
        imat <- diag(M)
        sm_id <- which((df.face[, 2] / n.mix) < 1e-3)
        if (any(sm_id)) {
          df.face <- df.face[-sm_id, , drop = FALSE]
          sector <- sector[-sm_id, , drop = FALSE]
        }
        nsec <- nrow(df.face)
        bsec <- df.face
        bsec[, 2] <- bsec[, 2] / sum(bsec[, 2])

        acov <- matrix(0, nrow = M, ncol = M)
        wtinv <- diag(1 / w)
        for (is in 1:nsec) {
          jvec <- sector[is, ]
          smat <- dp[, which(jvec == 1), drop = FALSE]
          wtinvs <- wtinv %*% smat
          pmat_is_p <- wtinv %*% smat %*% solve(t(smat) %*% wtinv %*% smat) %*% t(smat)
          pmat_is <- (imat - pmat_is_p)
          acov <- acov + bsec[is, 2] * pmat_is %*% v1 %*% t(pmat_is)
        }
      } else {
        acov <- v1
      }
      z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
      # z.mult = 2
      vsc <- diag(acov)
      hl <- z.mult * sqrt(vsc)
      lwr <- loweta <- etahat - hl
      upp <- uppeta <- etahat + hl

      vscu <- diag(v1)
      hlu <- z.mult * sqrt(vscu)
      lwru <- lowetau <- etahatu - hlu
      uppu <- uppetau <- etahatu + hlu

      if (is.matrix(xvs2)) {
        grid2 <- expand.grid(as.data.frame(xvs2))
      } else if (is.list(xvs2)) {
        grid2 <- expand.grid(xvs2)
      }

      acov_cp <- acov
      if (ne >= 1) {
        # muhatu = yvecu
        # new: if there's empty cells, just augment muhatu to include NA's
        etahatu_all <- 1:(M + zeros) * 0
        etahatu_all[zeros_ps] <- NA
        etahatu_all[obs_cells] <- etahatu
        etahatu <- etahatu_all

        ans_im <- .impute_em2(empty_cells, obs_cells, M = (M + zeros), etahatu, Nd, nd, Nhat, w,
                              domain.ord, grid2, Ds, sh, etahat, lwr, upp, vsc, acov, amat_0)
        etahat <- ans_im$muhat
        lwr <- ans_im$lwr
        upp <- ans_im$upp
        vsc <- ans_im$vsc
        acov_cp <- ans_im$acov_cp
        domain.ord <- ans_im$domain.ord
      }

      sig1 <- vsc
      sig2 <- NULL
      if (Nd >= 1 && nsm >= 1) {
        new_obs_cells <- sort(unique(c(empty_cells, obs_cells)))
        if (Nd == 1) {
          ans_im2 <- .impute_em2(small_cells, new_obs_cells,
                                 M = (M + zeros), etahatu,
                                 Nd, nd, Nhat, w, domain.ord, grid2, Ds, sh, etahat,
                                 lwr, upp, vsc, amat_0
          )
          sig2 <- ans_im2$vsc
        }
        if (Nd >= 2) {
          ans_im2 <- .impute_em3(small_cells, M = (M + zeros), etahatu, Nd, nd, Nhat, w, domain.ord, grid2,
                                 Ds, sh, etahat, lwr, upp, vsc, new_obs_cells, amat_0)
          sig2 <- ans_im2$vsc
        }
      }

      vsc_mix <- sig1
      vsc_mix_imp <- NULL
      if (Nd >= 1 && nsm >= 1) {
        imp_ps <- sort(c(zeros_ps, small_cells))
        nd_imp_ps <- nd[imp_ps]
        l_imp_ps <- length(nd_imp_ps)

        sig1_imp <- sig1[imp_ps]
        sig2_imp <- sig2[imp_ps]

        vsc_mix_imp <- 1:l_imp_ps * 0
        # tmp:
        if (Nd > 1) {
          nbors_lst_0 <- .fn_nbors_2(small_cells, amat_0, etahat)
        }

        k_sm <- 1
        for (k in 1:l_imp_ps) {
          ndek <- nd_imp_ps[k]
          s1k <- sig1_imp[k]
          s2k <- sig2_imp[k]
          # tmp
          if (Nd > 1 && ndek > 0) {
            max_ps <- nbors_lst_0$nb_lst_maxs[[k_sm]]
            min_ps <- nbors_lst_0$nb_lst_mins[[k_sm]]
            # max_ps = min_ps = NULL
            if (length(max_ps) > 0) {
              mu_max <- etahat[max_ps]
            }
            if (length(min_ps) > 0) {
              mu_min <- etahat[min_ps]
            }
            k_sm <- k_sm + 1
          } else if (Nd == 1 || ndek == 0) {
            min_ps <- max_ps <- 1
          }

          # when min_ps or max_ps = NULL: kth cell is in a 2D case, small cell, at the edge of the
          # grid and there is only one constraint in amat
          if (length(min_ps) == 0 || length(max_ps) == 0) {
            varek <- s1k
          } else {
            if (ndek == 0) {
              # varek = s1k
              varek <- s2k
            } else if (ndek == 1) {
              varek <- s1k * .2 + s2k * .8
            } else if (ndek %in% c(2, 3)) {
              varek <- s1k * .3 + s2k * .7
            } else if (ndek %in% c(4, 5)) {
              varek <- s1k * .6 + s2k * .4
            } else if (ndek %in% c(6, 7)) {
              varek <- s1k * .6 + s2k * .4
            } else if (ndek %in% c(8, 9, 10)) {
              varek <- s1k * .8 + s2k * .2
            } # else if(ndek == 10){
            #  varek = s1k*.8 + s2k*.2
            # }
          }
          vsc_mix_imp[k] <- varek
        }
        vsc_mix[imp_ps] <- vsc_mix_imp
      }
    }
    # eta = etahat[ps]; se = vsc[ps]
    # learnt from predict.svyglm:
    z.mult <- qnorm((1 - level) / 2, lower.tail = FALSE)
    fit <- etahat[ps]
    se.fit <- sqrt(vsc_mix[ps])
    lwr <- etahat[ps] - z.mult * se.fit
    upp <- etahat[ps] + z.mult * se.fit

    fit <- switch(type,
                  link = fit,
                  response = object$family$linkinv(fit)
    )
    lwr <- switch(type,
                  link = lwr,
                  response = object$family$linkinv(lwr)
    )
    upp <- switch(type,
                  link = upp,
                  response = object$family$linkinv(upp)
    )
    se.fit <- switch(type,
                     link = se.fit,
                     response = object$family$linkinv(se.fit)
    )

    # no unconstrained:...
    if (type == "link") {
      ans <- list(fit = fit, lwr = lwr, upp = upp)#, se.fit = se.fit)
    }
    if (type == "response") {
      ans <- list(fit = fit, lwr = lwr, upp = upp)
    }
  }
  return(ans)
}

#############################################
# print method
#############################################
print.csvy <- function(x,...){
  print(x$survey.design, varnames = FALSE, design.summaries = FALSE,...)
  NextMethod()
}


###################################################################################
# extract mixture variance-covariance matrix
# basically the same as vcov.svyglm
# vcov.svyby is more complicated but seems unnecessary
# common for csvby and csvyglm?
###################################################################################
vcov.csvy <- function(object,...){
  v <- object$acov_cp
  dimnames(v) <- list(names(coef(object)), names(coef(object)))
  return(v)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#routines modified from methods of svyby or svyglm, svyrepglm
#---------------------------------------------------------------------------------------------------------------------------------------
# dotchart.csvy <- function(x,...,pch = 19){
#   xx <- x$ans.unc_cp
#   dotchart(xx,...,pch = pch)
# }

#---------------------------------------------------------------------------------------------------------------------------------------
deff.csvy <- function(object,...) {
  x <- object$ans.unc_cp
  deff(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
barplot.csvy <- function(height, beside = TRUE,...){
  x <- height$ans.unc_cp
  #survey:::barplot.svyby(xx, beside = TRUE,...)
  barplot(x, beside = TRUE,...)
}

#plot.csvy <- barplot.csvy
#---------------------------------------------------------------------------------------------------------------------------------------
SE.csvy <- function(object,...){
  x <- object$ans.unc_cp
  SE(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
coef.csvy <- function(object, ...) {
  x <- object$ans.unc_cp
  coef(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
svycontrast.csvy <- function(stat, contrasts,...) {
  x <- stat$ans.unc_cp
  svycontrast(x, contrasts,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#fix
ftable.csvy <- function(x,...) {
  xx <- x$ans.unc_cp
  ftable(xx,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
# confint.csvyby <- function(object, parm, level = 0.95, df = Inf,...){
#   x <- object$ans.unc_cp
#   confint(x, parm, level, df,...)
# }

#---------------------------------------------------------------------------------------------------------------------------------------
# deff.csvyby <- function(object,...){
#   x <- object$ans.unc_cp
#   if(missing(x) | is.null(x)){
#     stop("only works when the fit is a class of svyby!")
#   }
#   deff(x,...)
# }


#---------------------------------------------------------------------------------------------------------------------------------------
#subroutines not exported
#---------------------------------------------------------------------------------------------------------------------------------------
#########################################
# subroutines for confidence interval
#########################################
.makebin <- function(x) {
  k <- length(x)
  r <- 0
  for (i in 1:k) {
    r <- r + x[k - i + 1] * 2^(i - 1)
  }
  r
}

#########################
# empty cell imputation
##########################
.impute_em2 <- function(empty_cells, obs_cells, M, yvec, Nd, nd, Nhat, w,
                        domain.ord, grid2, Ds = NULL, sh = NULL, muhat = NULL,
                        lwr = NULL, upp = NULL, vsc = NULL, acov = NULL, amat_0 = NULL) {
  ne <- length(empty_cells)
  # observed sample size in cells with 0 or 1 n
  nd_ne <- nd[empty_cells]
  zeros <- ones <- 0
  zeros_ps <- NULL
  if (ne >= 1) {
    # ignore nd=1 now
    if (any(nd == 1)) {
      ones_ps <- which(nd == 1)
      ones <- length(ones_ps)
      #   yvec_ones = yvec[which(obs_cells %in% ones_ps)]
      #   Nhat_ones = Nhat[which(obs_cells %in% ones_ps)]
    }
    if (any(nd == 0)) {
      zeros_ps <- which(nd == 0)
      zeros <- length(zeros_ps)
    }
    zeros <- ne - ones
  }

  lc <- rc <- 1:ne * 0
  acov_cp <- NULL
  if(!is.null(acov)) {
    acov_cp <- matrix(0, M, M)
    acov_cp[obs_cells, obs_cells] <- acov
  }
  #new: imputate covariance for empty cells: average of l cell and r cell
  if (Nd == 1) {
    vsc_all <- lwr_all <- upp_all <- muhat_all <- domain.ord_all <- 1:M * 0
    nslst <- pslst <- vector("list", ne)
    vsce <- lwre <- uppe <- muhate <- domain.ord_em <- NULL
    cve <- vector("list", length = ne)
    for (k in 1:ne) {
      e <- empty_cells[k]
      ndek <- nd_ne[k]
      # one_k = 1
      # new: consider n = 1
      # at_max = all(obs_cells < e | (ndek == 1 & e == max(obs_cells))
      # at_min = all(obs_cells > e | (ndek == 1 & e == min(obs_cells))
      at_max <- all(obs_cells < e | e == max(obs_cells))
      at_min <- all(obs_cells > e | e == min(obs_cells))
      if (at_max) {
        lc <- max(obs_cells[obs_cells < e])
        rc <- lc
        dmem <- obs_cells[which(obs_cells %in% lc)]
      }
      if (at_min) {
        rc <- min(obs_cells[obs_cells > e])
        lc <- rc
        dmem <- obs_cells[which(obs_cells %in% rc)]
      }
      if (!any(at_max) & !any(at_min)) {
        lc <- max(obs_cells[obs_cells < e])
        rc <- min(obs_cells[obs_cells > e])
        dmem <- obs_cells[which(obs_cells %in% lc)]
      }
      n_lc <- nd[lc]
      n_rc <- nd[rc]

      # yl = yvec[which(obs_cells%in%lc)]
      # yr = yvec[which(obs_cells%in%rc)]
      #
      # Nl = Nhat[lc]
      # Nr = Nhat[rc]
      # new: switch lc and rc if sh=2; add 6 and 8
      if (sh == 2 | sh == 6 | sh == 8) {
        lc_0 <- rc
        rc_0 <- lc
        lc <- lc_0
        rc <- rc_0
      }
      if (!is.null(muhat)) {
        mul <- muhat[which(obs_cells %in% lc)]
        mur <- muhat[which(obs_cells %in% rc)]
      }

      # variance
      if (!is.null(vsc)) {
        vl <- vsc[which(obs_cells %in% lc)]
        vr <- vsc[which(obs_cells %in% rc)]
      }

      #new: covariance
      if(!is.null(acov)) {
        cvl <- acov[which(obs_cells %in% lc), ]
        cvr <- acov[which(obs_cells %in% rc), ]
        cvek <- (cvl + cvr) / 2
        cve[[k]] <- cvek
      }

      # muhate = c(muhate, (mul + mur)/2)
      # lwre = c(lwre, (mul + mur)/2 - 2 * sqrt(vl))
      # uppe = c(uppe, (mul + mur)/2 + 2 * sqrt(vr))

      if (!is.null(muhat)) {
        muhatek <- (mul - 2 * sqrt(vl) + mur + 2 * sqrt(vr)) / 2
        # check if muhatek follows the shape constraint
        lr <- c(lc, rc)
        # sh=10 is for block ave
        # sh=9 is for block ord
        #new: add incr.conc 7, incr.conv 5
        if (sh %in% c(1, 5, 7, 9, 10)) {
          if (length(lr) > 1) {
            bool <- mul <= muhatek & mur >= muhatek
          } else if (length(lr) == 1) {
            if (at_min) {
              bool <- mur >= muhatek
            } else if (at_max) {
              bool <- mul <= muhatek
            }
          }
        }
        #new: add decr.conc 8 and decr.conv 6
        if (sh %in% c(2, 6, 8)) {
          if (length(lr) > 1) {
            bool <- mul >= muhatek & mur <= muhatek
          } else if (length(lr) == 1) {
            if (at_min) {
              bool <- mur <= muhatek
            } else if (at_max) {
              bool <- mul >= muhatek
            }
          }
        }
        if (!bool) {
          if (!any(at_max) & !any(at_min)) {
            muhatek <- (mul + mur) / 2
          } else if (at_min) {
            muhatek <- mul
          } else if (at_max) {
            muhatek <- mur
          }
        }
        muhate <- c(muhate, muhatek)
      }

      if (!is.null(vsc)) {
        # varek = (mur + 2*sqrt(vr) - mul + 2*sqrt(vl))^2 / 16
        rpt <- mur + 2 * sqrt(vr)
        lpt <- mul - 2 * sqrt(vl)
        varek <- (rpt - lpt)^2 / 16
        vsce <- c(vsce, varek)
        if (at_min & (ndek == 0)) {
          lwrek <- -1e+5
          uppek <- muhatek + 2 * sqrt(varek)
        } else if (at_max & (ndek == 0)) {
          uppek <- 1e+5
          lwrek <- muhatek - 2 * sqrt(varek)
        } else if (at_min & (ndek >= 1) & (ndek <= 10)) {
          varek_lwr <- max(vsc)
          lwrek <- muhatek - 2 * sqrt(varek_lwr)
          uppek <- muhatek + 2 * sqrt(varek)
        } else if (at_max & (ndek >= 1) & (ndek <= 10)) {
          varek_upp <- max(vsc)
          lwrek <- muhatek - 2 * sqrt(varek)
          uppek <- muhatek + 2 * sqrt(varek_upp)
        } else {
          lwrek <- muhatek - 2 * sqrt(varek)
          uppek <- muhatek + 2 * sqrt(varek)
        }
        lwre <- c(lwre, lwrek)
        uppe <- c(uppe, uppek)
      }

      ps <- unique(c(lc, rc))
      ns <- unique(c(n_lc, n_rc))
      pslst[[k]] <- ps
      nslst[[k]] <- ns
      domain.ord_em <- c(domain.ord_em, dmem)
    }

    domain.ord_all[obs_cells] <- domain.ord
    # domain.ord_all[which(nd > 0)] = domain.ord
    domain.ord_all[empty_cells] <- domain.ord_em
    #
    if (!is.null(muhat)) {
      muhat_all[obs_cells] <- muhat
      muhat_all[empty_cells] <- muhate
      muhat <- muhat_all
    }

    if (!is.null(vsc)) {
      lwr_all[obs_cells] <- lwr
      lwr_all[empty_cells] <- lwre

      upp_all[obs_cells] <- upp
      upp_all[empty_cells] <- uppe

      upp <- upp_all
      lwr <- lwr_all

      vsc_all[obs_cells] <- vsc
      vsc_all[empty_cells] <- vsce
      vsc <- vsc_all
    }

    #new:
    if(!is.null(acov)) {
      acov_cp[obs_cells, obs_cells] <- acov
      diag(acov_cp) <- vsc
      for (k in 1:ne) {
        e <- empty_cells[k]
        acov_cp[!(1:M) %in% empty_cells, e] <- cve[[k]]
        acov_cp[e, !(1:M) %in% empty_cells] <- cve[[k]]
        #acov_cp[!(1:M) %in% e, e] <- cve[[k]]
        #acov_cp[e, !(1:M) %in% e] <- cve[[k]]
      }
    }
    # yvec = yvec_all
    # w = w_all
    domain.ord <- domain.ord_all
  } else if (Nd >= 2) {
    vsc_all <- lwr_all <- upp_all <- muhat_all <- 1:M * 0
    yvec_all <- w_all <- domain.ord_all <- 1:M * 0
    empty_grids <- grid2[empty_cells, ]
    maxs <- apply(grid2, 2, max)
    mins <- apply(grid2, 2, min)
    at_max <- (empty_grids == maxs)
    at_min <- (empty_grids == mins)

    vsce <- lwre <- uppe <- muhate <- domain.ord_em <- NULL
    #new:
    cve <- vector("list", length = ne)

    ye <- we <- NULL
    nslst <- pslst <- vector("list", 2)
    nbors <- NULL
    nbors_lst_0 <- .fn_nbors_3(empty_cells, amat_0, muhat)
    for (k in 1:length(empty_cells)) {
      # new: if is new because of rm_id2
      nbors_k_maxs <- nbors_lst_0$nb_lst_maxs[[k]]
      nbors_k_mins <- nbors_lst_0$nb_lst_mins[[k]]
      em_k <- empty_cells[k]

      # revised: nbors_k_mins and nbors_k_maxs must be in obs_cells
      # if(length(nbors_k_mins) > 0 | length(nbors_k_maxs) > 0){
      if (length(nbors_k_mins) > 0 & any(nbors_k_mins %in% obs_cells) | length(nbors_k_maxs) > 0 & any(nbors_k_maxs %in% obs_cells)) {
        max_ps <- min_ps <- NULL
        if (length(nbors_k_maxs) > 0) {
          mu_max <- max(muhat[obs_cells %in% nbors_k_maxs])
          max_ps <- (nbors_k_maxs[which(muhat[obs_cells %in% nbors_k_maxs] == mu_max)])[1]
        } else {
          # mu_max = muhat[em_k]
          # max_ps = NULL
          mu_min <- min(muhat[obs_cells %in% nbors_k_mins])
          min_ps <- (nbors_k_mins[which(muhat[obs_cells %in% nbors_k_mins] == mu_min)])[1]
          mu_max <- mu_min
          max_ps <- min_ps
          dmem <- obs_cells[which(obs_cells %in% min_ps)]
        }

        if (length(nbors_k_mins) > 0) {
          mu_min <- min(muhat[obs_cells %in% nbors_k_mins])
          min_ps <- (nbors_k_mins[which(muhat[obs_cells %in% nbors_k_mins] == mu_min)])[1]
        } else {
          # mu_min = muhat[em_k]
          # min_ps = NULL
          mu_max <- max(muhat[obs_cells %in% nbors_k_maxs])
          max_ps <- (nbors_k_maxs[which(muhat[obs_cells %in% nbors_k_maxs] == mu_max)])[1]
          mu_min <- mu_max
          min_ps <- max_ps
          dmem <- obs_cells[which(obs_cells %in% max_ps)]
        }

        if (!is.null(max_ps)) {
          v_max <- vsc[which(obs_cells %in% max_ps)]
        }
        if (!is.null(min_ps)) {
          v_min <- vsc[which(obs_cells %in% min_ps)]
        }

        if (!is.null(max_ps) & !is.null(min_ps)) {
          varek_new <- (mu_max + 2 * sqrt(v_max) - mu_min + 2 * sqrt(v_min))^2 / 16
          muhatek <- (mu_min - 2 * sqrt(v_min) + mu_max + 2 * sqrt(v_max)) / 2
          # ?mu_min
          dmem <- obs_cells[which(obs_cells %in% min_ps)]
        }

        varek <- varek_new
        vsce <- c(vsce, varek)
        muhate <- c(muhate, muhatek)
        domain.ord_em <- c(domain.ord_em, dmem)

        if (length(nbors_k_mins) == 0) {
          lwrek <- -1e+5
          uppek <- muhatek + 2 * sqrt(varek)
        } else if (length(nbors_k_maxs) == 0) {
          uppek <- 1e+5
          lwrek <- muhatek - 2 * sqrt(varek)
        } else if (length(nbors_k_mins) > 0 & length(nbors_k_maxs) > 0) {
          uppek <- muhatek + 2 * sqrt(varek)
          lwrek <- muhatek - 2 * sqrt(varek)
        }

        uppe <- c(uppe, uppek)
        lwre <- c(lwre, lwrek)

        #new: covariance
        if(!is.null(acov)) {
          cv_max <- acov[which(obs_cells %in% max_ps), ]
          cv_min <- acov[which(obs_cells %in% min_ps), ]
          cvek <- (cv_max + cv_min) / 2
          cve[[k]] <- cvek
        }
      } else {
        # ye = c(ye, rev(ye)[1])
        # we = c(we, rev(we)[1])
        # need to test more
        if (!is.null(rev(domain.ord_em)[1])) {
          domain.ord_em <- c(domain.ord_em, rev(domain.ord_em)[1])
        }
        if (!is.null(muhat)) {
          muhate <- c(muhate, rev(muhate)[1])
        }
        if (!is.null(vsc)) {
          vsce <- c(vsce, rev(vsce)[1])
          lwre <- c(lwre, rev(lwre)[1])
          uppe <- c(uppe, rev(uppe)[1])
        }
        #new: covariance
        #need to test more
        if(!is.null(acov)){
          cve[[k]] <- cvek[[k-1]]
        }
      }
    }

    domain.ord_all[obs_cells] <- domain.ord
    # temp:
    # domain.ord_em = obs_cells[which(obs_cells%in%ps[1:length(empty_cells)])]
    domain.ord_all[empty_cells] <- domain.ord_em
    domain.ord <- domain.ord_all

    if (!is.null(muhat)) {
      muhat_all[obs_cells] <- muhat
      muhat_all[empty_cells] <- muhate
      muhat <- muhat_all
    }

    if (!is.null(vsc)) {
      lwr_all[obs_cells] <- lwr
      lwr_all[empty_cells] <- lwre
      upp_all[obs_cells] <- upp
      upp_all[empty_cells] <- uppe
      vsc_all[obs_cells] <- vsc
      vsc_all[empty_cells] <- vsce

      upp <- upp_all
      lwr <- lwr_all
      vsc <- vsc_all
    }

    #new:
    if(!is.null(acov)) {
      acov_cp[obs_cells, obs_cells] <- acov
      diag(acov_cp) <- vsc
      for (k in 1:ne) {
        em_k <- empty_cells[k]
        acov_cp[!(1:M) %in% empty_cells, em_k] <- cve[[k]]
        acov_cp[em_k, !(1:M) %in% empty_cells] <- cve[[k]]
      }
    }

  }

  rslt <- list(muhat = muhat, lwr = lwr, upp = upp, domain.ord = domain.ord, vsc = vsc, acov_cp = acov_cp)
  # rslt = list(ps = ps, ns = ns, muhat = muhat, lwr = lwr, upp = upp, domain.ord = domain.ord, yvec = yvec, w = w)
  return(rslt)
}

###############################
# subrountine to find nbors (>=2D)
###############################
.fn_nbors <- function(empty_grids, grid2, mins, maxs) {
  nr <- nrow(empty_grids)
  nb_lst <- vector("list", length = nr)
  to_merge <- list()
  iter <- 1
  for (k in 1:nr) {
    # ndek = nd_ne[k]
    # axes = apply(empty_grids[k,], 2, function(e) e + c(-1,1))
    nbors_k <- NULL
    ptk <- empty_grids[k, ]
    ptk <- as.numeric(ptk)
    nc <- ncol(grid2)
    for (i in 1:nc) {
      pt1 <- ptk[i]
      pt2 <- ptk[-i]
      axes <- pt1 + c(-1, 1)
      pair_i <- matrix(0, nrow = 2, ncol = nc)
      pair_i[1, i] <- axes[1]
      pair_i[1, -i] <- pt2
      pair_i[2, i] <- axes[2]
      pair_i[2, -i] <- pt2
      nbors_k <- rbind(nbors_k, pair_i)
    }
    rm_id <- unique(c(which(apply(nbors_k, 1, function(e) any(e < mins))), which(apply(nbors_k, 1, function(e) any(e > maxs)))))
    if (length(rm_id) > 0) {
      nbors_k <- nbors_k[-rm_id, , drop = FALSE]
    }
    # new: label the cells whose neighbors should be merged later
    # rm_id2 = NULL
    if (nrow(nbors_k) > 0) {
      for (h in 1:nrow(nbors_k)) {
        nbh <- nbors_k[h, ]
        bool <- any(apply(empty_grids, 1, function(e) all(e == nbh)))
        if (bool) {
          # rm_id2 = c(rm_id2, h)
          ps <- as.numeric(which(apply(empty_grids, 1, function(e) all(e == nbh))))
          to_merge[[iter]] <- sort(c(k, ps))
          iter <- iter + 1
        }
      }
    }
    nb_lst[[k]] <- nbors_k
  }

  # need to be more efficient
  to_merge <- unique(to_merge)
  nm <- length(to_merge)
  if (nm > 1) {
    for (i in 1:(nm - 1)) {
      to_merge_i <- to_merge[[i]]
      vec <- to_merge_i
      vec_ps <- i
      for (j in (i + 1):nm) {
        to_merge_j <- to_merge[[j]]
        if (any(to_merge_j %in% to_merge_i)) {
          vec <- sort(unique(c(vec, to_merge_j)))
          vec_ps <- c(vec_ps, j)
        }
      }
      if (length(vec_ps) > 1) {
        to_merge[vec_ps] <- lapply(to_merge[vec_ps], function(x) x <- vec)
      }
    }
    to_merge <- unique(to_merge)
  }

  if (nm > 0) {
    for (i in 1:length(to_merge)) {
      ps <- to_merge[[i]]
      nbors_ps <- unique(do.call(rbind, nb_lst[ps]))
      rm_id2 <- NULL
      for (h in 1:nrow(nbors_ps)) {
        nbh <- nbors_ps[h, ]
        bool <- any(apply(empty_grids, 1, function(e) all(e == nbh)))
        if (bool) {
          rm_id2 <- c(rm_id2, h)
        }
      }
      if (length(rm_id2) > 0) {
        nbors_ps <- nbors_ps[-rm_id2, , drop = FALSE]
      }
      nb_lst[ps] <- lapply(nb_lst[ps], function(x) x <- nbors_ps)
    }
  }
  return(nb_lst)
}

####################################################################################
# compute the variance for n=2...10 cells
# keep muhat, only compute weighted variance
####################################################################################
.impute_em3 <- function(small_cells, M, yvec, Nd, nd, Nhat, w, domain.ord, grid2,
                        Ds = NULL, sh = NULL, muhat = NULL, lwr = NULL, upp = NULL,
                        vsc = NULL, new_obs_cells = NULL, amat_0 = NULL) {
  ne <- length(small_cells)
  # observed sample size for n=2...10 cells
  # the `small_cells` are not empty
  nd_ne <- nd[small_cells]
  # vsc_ne = vsc[which(new_obs_cells%in%small_cells)]
  vsc_ne <- vsc[small_cells]

  lc <- rc <- 1:ne * 0
  if (Nd >= 2) {
    vsc_all <- lwr_all <- upp_all <- muhat_all <- 1:M * 0
    yvec_all <- w_all <- domain.ord_all <- 1:M * 0
    empty_grids <- grid2[small_cells, ]
    maxs <- apply(grid2, 2, max)
    mins <- apply(grid2, 2, min)
    at_max <- (empty_grids == maxs)
    at_min <- (empty_grids == mins)

    vsce <- lwre <- uppe <- muhate <- domain.ord_em <- NULL
    ye <- we <- NULL
    nslst <- pslst <- vector("list", 2)
    nbors <- NULL
    nbors_lst_0 <- .fn_nbors_2(small_cells, amat_0, muhat)
    for (k in 1:length(small_cells)) {
      # new: if ... is new because of rm_id2
      # nbors_k_maxs = nbors_lst_0$nb_lst_maxs[[k]]
      # nbors_k_mins = nbors_lst_0$nb_lst_mins[[k]]
      #max_ps <- (nbors_lst_0$nb_lst_maxs[[k]])[1]
      #min_ps <- (nbors_lst_0$nb_lst_mins[[k]])[1]
      #new: avoid the case when there is only integer(0) and [1] will give NA
      if (length(nbors_lst_0$nb_lst_maxs[[k]]) >= 1) {
        max_ps <- nbors_lst_0$nb_lst_maxs[[k]][1]
      } else if (length(nbors_lst_0$nb_lst_maxs[[k]]) == 0) {
        max_ps <- nbors_lst_0$nb_lst_maxs[[k]]
      }
      if (length(nbors_lst_0$nb_lst_mins[[k]]) >= 1) {
        min_ps <- nbors_lst_0$nb_lst_mins[[k]][1]
      } else if (length(nbors_lst_0$nb_lst_mins[[k]]) == 0) {
        min_ps <- nbors_lst_0$nb_lst_mins[[k]]
      }

      sm_k <- small_cells[k]
      if (length(max_ps) > 0) {
        mu_max <- muhat[max_ps]
        v_max <- vsc[max_ps]
      }
      if (length(min_ps) > 0) {
        mu_min <- muhat[min_ps]
        v_min <- vsc[min_ps]
      }
      # new:
      vsc_nek <- vsc_ne[k]
      if (length(min_ps) > 0 & length(max_ps) > 0) {
        varek_new <- (mu_max + 2 * sqrt(v_max) - mu_min + 2 * sqrt(v_min))^2 / 16
      } else {
        varek_new <- vsc_nek
      }
      # the variance by survey
      # the if else can be replaced by vectorized computation?
      vsce <- c(vsce, varek_new)
      # print(c(k, length(vsce)))
    }

    if (!is.null(vsc)) {
      vsc_all <- vsc
      vsc_all[small_cells] <- vsce
      vsc <- vsc_all
    }
  }
  # rslt = list(lwr = lwr, upp = upp, vsc = vsc)
  rslt <- list(vsc = vsc)
  return(rslt)
}

###############################
# subrountine to find nbors (>=2D)
###############################
.fn_nbors_2 <- function(small_cells, amat, muhat) {
  nr <- length(small_cells)
  nb_lst_mins <- nb_lst_maxs <- vector("list", length = nr)
  # to_merge = list()
  # iter = 1
  for (k in 1:nr) {
    ptk <- small_cells[k]
    col_ptk <- amat[, ptk]
    rows_ps <- which(col_ptk != 0)
    amat_sub <- amat[rows_ps, , drop = FALSE]
    # nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0)) %>% as.vector() %>% unique()
    nbors_k <- apply(amat_sub, 1, function(.x) which(.x != 0))

    #for a constraint: -1,2,-1, need to remove the 2nd -1
    if(is.list(nbors_k)){
      nbors_k <- NULL
      sgn_ptk <- as.vector(amat[rows_ps, ptk, drop = FALSE])
      for(ir in 1:nrow(amat_sub)){
        amat_sub_ir <- amat_sub[ir, ]
        sng_ptk_ir <- sgn_ptk[ir]
        for(ic in 1:ncol(amat_sub)){
          if (ic == ptk) {
            nbors_k <- c(nbors_k, ptk)
          } else if (amat_sub_ir[ic] != sng_ptk_ir & amat_sub_ir[ic] != 0) {
            nbors_k <- c(nbors_k, ic)
          }
        }
      }
    }

    nbors_k <- unique(as.vector(nbors_k))
    nbors_k <- nbors_k[-which(nbors_k == ptk)]
    # 1 and -1 show the position of the cell being imputed
    # 1: to search for min candidates; -1: to search for max candidates
    sgn_nbors_k <- col_ptk[which(col_ptk != 0)]
    #nbors_k_mins <- nbors_k[sgn_nbors_k == 1]
    #nbors_k_maxs <- nbors_k[sgn_nbors_k == -1]
    #new: include conc or incr.conc; if conv or decr, then let amat = -amat
    nbors_k_mins <- nbors_k[sgn_nbors_k == 1 | sgn_nbors_k == 2]
    nbors_k_maxs <- nbors_k[sgn_nbors_k == -1]
    nbors_k_min <- nbors_k_max <- integer(0L)
    if (length(nbors_k_maxs) > 0) {
      max_mu <- max(muhat[nbors_k_maxs])
      nbors_k_max <- nbors_k_maxs[which(muhat[nbors_k_maxs] == max_mu)]
    }
    if (length(nbors_k_mins) > 0) {
      min_mu <- min(muhat[nbors_k_mins])
      nbors_k_min <- nbors_k_mins[which(muhat[nbors_k_mins] == min_mu)]
    }

    nb_lst_mins[[k]] <- nbors_k_min
    nb_lst_maxs[[k]] <- nbors_k_max
  }

  rslt <- list(nb_lst_maxs = nb_lst_maxs, nb_lst_mins = nb_lst_mins)
  return(rslt)
}

##############################################################
# subrountine to find nbors (>=2D)
# new neighbor routine for >= 2D empty cells
##############################################################
.fn_nbors_3 <- function(empty_cells, amat_0, muhat) {
  nr <- length(empty_cells)
  nb_lst_mins <- nb_lst_maxs <- vector("list", length = nr)
  to_merge_mins <- to_merge_maxs <- list()
  iter <- 1
  for (k in 1:nr) {
    nbors_k <- NULL
    ptk <- empty_cells[k]
    col_ptk <- amat_0[, ptk]
    rows_ps <- which(col_ptk != 0)
    amat_sub <- amat_0[rows_ps, , drop = FALSE]
    # nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0)) %>% as.vector() %>% unique()
    nbors_k <- apply(amat_sub, 1, function(.x) which(.x != 0))
    #new: for conc: -1,2,-1 will give three cells and -1,1 will give two cells and this creates a list
    #if(is.list(nbors_k)){
    #  nbors_k <- unlist(nbors_k)
    #}
    #for a constraint: -1,2,-1, need to remove the 2nd -1
    if(is.list(nbors_k)){
      nbors_k <- NULL
      sgn_ptk <- as.vector(amat_0[rows_ps, ptk, drop = FALSE])
      for(ir in 1:nrow(amat_sub)){
        amat_sub_ir <- amat_sub[ir, ]
        sng_ptk_ir <- sgn_ptk[ir]
        for(ic in 1:ncol(amat_sub)){
          if (ic == ptk) {
            nbors_k <- c(nbors_k, ptk)
          } else if (amat_sub_ir[ic] != sng_ptk_ir & amat_sub_ir[ic] != 0) {
            nbors_k <- c(nbors_k, ic)
          }
        }
      }
    }

    nbors_k <- unique(as.vector(nbors_k))
    nbors_k <- nbors_k[-which(nbors_k == ptk)]
    # 1 and -1 show the position of the cell being imputed
    # 1: to search for min candidates; -1: to search for max candidates
    sgn_nbors_k <- col_ptk[which(col_ptk != 0)]
    #new: include conc or incr.conc; if conv or decr, then let amat = -amat
    nbors_k_mins <- nbors_k[sgn_nbors_k == 1 | sgn_nbors_k == 2]
    #nbors_k_mins <- nbors_k[sgn_nbors_k == 1]
    nbors_k_maxs <- nbors_k[sgn_nbors_k == -1]

    # new: need to remove empty_cells from nbors_k_maxs and nbors_k_mins
    nbors_k_mins <- setdiff(nbors_k_mins, empty_cells)
    nbors_k_maxs <- setdiff(nbors_k_maxs, empty_cells)

    if (length(nbors_k_maxs) > 0) {
      nbors_k_maxs <- nbors_k_maxs[which(muhat[nbors_k_maxs] == max(muhat[nbors_k_maxs]))]
    }
    if (length(nbors_k_mins) > 0) {
      nbors_k_mins <- nbors_k_mins[which(muhat[nbors_k_mins] == min(muhat[nbors_k_mins]))]
    }

    # new: label the cells whose neighbors should be merged later
    rm_id2 <- NULL
    if (length(nbors_k_maxs) > 0) {
      for (h in 1:length(nbors_k_maxs)) {
        nbh <- nbors_k_maxs[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          ps <- which(empty_cells == nbh)
          to_merge_maxs[[iter]] <- sort(c(k, ps))
          iter <- iter + 1
        }
      }
    }

    if (length(nbors_k_mins) > 0) {
      for (h in 1:length(nbors_k_mins)) {
        nbh <- nbors_k_mins[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          ps <- which(empty_cells == nbh)
          to_merge_mins[[iter]] <- sort(c(k, ps))
          iter <- iter + 1
        }
      }
    }

    nb_lst_mins[[k]] <- nbors_k_mins
    nb_lst_maxs[[k]] <- nbors_k_maxs
  }

  # need to be more efficient
  to_merge_maxs <- unique(to_merge_maxs)
  to_merge_mins <- unique(to_merge_mins)
  nm_maxs <- length(to_merge_maxs)
  nm_mins <- length(to_merge_mins)
  if (nm_maxs > 1) {
    for (i in 1:(nm_maxs - 1)) {
      to_merge_i <- to_merge_maxs[[i]]
      vec <- to_merge_i
      vec_ps <- i
      for (j in (i + 1):nm_maxs) {
        to_merge_j <- to_merge_maxs[[j]]
        if (any(to_merge_j %in% to_merge_i)) {
          vec <- sort(unique(c(vec, to_merge_j)))
          vec_ps <- c(vec_ps, j)
        }
      }
      if (length(vec_ps) > 1) {
        to_merge_maxs[vec_ps] <- lapply(to_merge_maxs[vec_ps], function(x) x <- vec)
      }
    }
    to_merge_maxs <- unique(to_merge_maxs)
  }

  if (nm_maxs > 0) {
    for (i in 1:length(to_merge_maxs)) {
      ps <- to_merge_maxs[[i]]
      # nbors_ps_maxs = unlist(nb_lst_maxs[ps]) %>% unique()
      nbors_ps_maxs <- unique(unlist(nb_lst_maxs[ps]))
      rm_id2_maxs <- NULL
      for (h in 1:length(nbors_ps_maxs)) {
        nbh <- nbors_ps_maxs[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_maxs <- c(rm_id2_maxs, h)
        }
      }
      if (length(rm_id2_maxs) > 0) {
        nbors_ps_maxs <- nbors_ps_maxs[-rm_id2_maxs]
      }
      nb_lst_maxs[ps] <- lapply(nb_lst_maxs[ps], function(x) x <- nbors_ps_maxs)

      ###
      for (h in 1:length(nbors_ps_maxs)) {
        nbh <- nbors_ps_maxs[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_maxs <- c(rm_id2_maxs, h)
        }
      }
      if (length(rm_id2_maxs) > 0) {
        nbors_ps_maxs <- nbors_ps_maxs[-rm_id2_maxs]
      }
      nb_lst_maxs[ps] <- lapply(nb_lst_maxs[ps], function(x) x <- nbors_ps_maxs)
    }
  }

  if (nm_mins > 0) {
    for (i in 1:length(to_merge_mins)) {
      ps <- to_merge_mins[[i]]
      # nbors_ps_mins = unlist(nb_lst_mins[ps]) %>% unique()
      nbors_ps_mins <- unique(unlist(nb_lst_mins[ps]))
      rm_id2_mins <- NULL
      for (h in 1:length(nbors_ps_mins)) {
        nbh <- nbors_ps_mins[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_mins <- c(rm_id2_mins, h)
        }
      }
      if (length(rm_id2_mins) > 0) {
        nbors_ps_mins <- nbors_ps_mins[-rm_id2_mins]
      }
      nb_lst_mins[ps] <- lapply(nb_lst_mins[ps], function(x) x <- nbors_ps_mins)

      ###
      for (h in 1:length(nbors_ps_mins)) {
        nbh <- nbors_ps_mins[h]
        bool <- any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_mins <- c(rm_id2_mins, h)
        }
      }
      if (length(rm_id2_mins) > 0) {
        nbors_ps_mins <- nbors_ps_mins[-rm_id2_mins]
      }
      nb_lst_mins[ps] <- lapply(nb_lst_mins[ps], function(x) x <- nbors_ps_mins)
    }
  }
  rslt <- list(nb_lst_maxs = nb_lst_maxs, nb_lst_mins = nb_lst_mins)
  return(rslt)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#subroutines not exported
#---------------------------------------------------------------------------------------------------------------------------------------
#combine all makeamat functions into one
.MakeAmat <- function(sh, Nd, block.ave.lst, block.ord.lst, zeros_ps,
                      M, M_0, noord_cells, xvs2, grid2, Ds) {
  #hkeep <- NULL
  zeros <- length(zeros_ps)
  amat_0 <- NULL # to be used in Nd > 1
  #if (is.null(amat)) {
  if (Nd == 1) {
    # sh = 0 means user must define amat
    if (sh > 0) {
      # Ds is not correct when >= 1 empty cell
      # new: try to find the constrained est. first
      block.ave <- NULL
      block.ord <- NULL
      if (length(block.ave.lst) > 0) {
        # print ('a')
        block.ave <- block.ave.lst[[1]]
      }
      if (length(block.ord.lst) > 0) {
        # print ('a')
        block.ord <- block.ord.lst[[1]]
      }
      if (length(zeros_ps) > 0) {
        M2 <- seq(1, M + zeros)[-zeros_ps]
      } else {
        M2 <- 1:M
      }
      amat <- .makeamat(
        x = M2, sh = sh, Ds = M, interp = TRUE, relaxes = FALSE,
        block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
        zeros_ps = zeros_ps, noord_cells = noord_cells, xvs2 = xvs2, grid2 = grid2
      )
    } #else if (sh == 0 && is.null(amat)) {
    # print (sh)
    #stop("User must define a constraint matrix!")
    #}
  } else if (Nd >= 2) {
    # temp:
    if (any(sh > 0)) {
      amat <- .makeamat_2D(
        x = NULL, sh = sh, grid2 = grid2, xvs2 = xvs2,
        zeros_ps = zeros_ps, noord_cells = noord_cells,
        Ds = Ds, interp = TRUE,
        block.ave.lst = block.ave.lst,
        block.ord.lst = block.ord.lst, M_0 = M_0
      )
      # to be used in impute_2 or impute_3
      amat_0 <- .makeamat_2D(
        x = NULL, sh = sh, grid2 = grid2, xvs2 = xvs2,
        zeros_ps = NULL, noord_cells = noord_cells,
        Ds = Ds, interp = TRUE,
        block.ave.lst = block.ave.lst,
        block.ord.lst = block.ord.lst, M_0 = M_0
      )
    } #else if (all(sh == 0) && is.null(amat)) {
    #stop("User must define a constraint matrix!")
    #}
  }
  #}
  ans <- list(amat = amat, amat_0 = amat_0)
  return(ans)
}


######
# 1D
######
.makeamat <- function(x, sh, Ds = NULL, suppre = FALSE, interp = FALSE,
                      relaxes = FALSE, h = NULL, block.ave = NULL,
                      block.ord = NULL, zeros_ps = NULL, noord_cells = NULL,
                      D_index = 1, xvs2 = NULL, grid2 = NULL) {
  n <- length(x)
  xu <- sort(unique(x))
  n1 <- length(xu)
  sm <- 1e-7
  ms <- NULL
  obs <- 1:n
  Nd <- length(Ds)
  M <- rev(x)[1]
  if (relaxes) {
    # create a list of amat's first
    d <- l <- x
    amat <- matrix(0, nrow = (M - 1), ncol = M)
    for (k in 1:(M - 1)) {
      amat[k, ] <- exp(-abs(l - (k + 1)) / h) / sum(exp(-abs(d - (k + 1)) / h)) - exp(-abs(l - k) / h) / sum(exp(-abs(d - k) / h))
    }
  } else {
    # block.ave and block.ord = -1 if sh is between 1 and 9
    if (all(block.ave == -1) & all(block.ord == -1)) {
      # print ('a')
      if (sh < 3) {
        # if sh=0, then keep the zero matrix
        amat <- matrix(0, nrow = n1 - 1, ncol = n)
        if (sh > 0) {
          for (i in 1:(n1 - 1)) {
            c1 <- min(obs[abs(x - xu[i]) < sm])
            c2 <- min(obs[abs(x - xu[i + 1]) < sm])
            amat[i, c1] <- -1
            amat[i, c2] <- 1
          }
          if (sh == 2) {
            amat <- -amat
          }
        }
        # return (amat)
      } else if (sh == 3 | sh == 4) {
        #  convex or concave
        amat <- matrix(0, nrow = n1 - 2, ncol = n)
        for (i in 1:(n1 - 2)) {
          c1 <- min(obs[x == xu[i]])
          c2 <- min(obs[x == xu[i + 1]])
          c3 <- min(obs[x == xu[i + 2]])
          amat[i, c1] <- xu[i + 2] - xu[i + 1]
          amat[i, c2] <- xu[i] - xu[i + 2]
          amat[i, c3] <- xu[i + 1] - xu[i]
        }
        if (sh == 4) {
          amat <- -amat
        }
        # return (amat)
      } else if (sh > 4 & sh < 9) {
        amat <- matrix(0, nrow = n1 - 1, ncol = n)
        for (i in 1:(n1 - 2)) {
          c1 <- min(obs[x == xu[i]])
          c2 <- min(obs[x == xu[i + 1]])
          c3 <- min(obs[x == xu[i + 2]])
          amat[i, c1] <- xu[i + 2] - xu[i + 1]
          amat[i, c2] <- xu[i] - xu[i + 2]
          amat[i, c3] <- xu[i + 1] - xu[i]
        }
        if (sh == 5) { ### increasing convex
          c1 <- min(obs[x == xu[1]])
          c2 <- min(obs[x == xu[2]])
          amat[n1 - 1, c1] <- -1
          amat[n1 - 1, c2] <- 1
          # return (amat)
        } else if (sh == 6) { ## decreasing convex
          c1 <- min(obs[x == xu[n1]])
          c2 <- min(obs[x == xu[n1 - 1]])
          amat[n1 - 1, c1] <- -1
          amat[n1 - 1, c2] <- 1
          # return (amat)
        } else if (sh == 7) { ## increasing concave
          amat <- -amat
          c1 <- min(obs[x == xu[n1]])
          c2 <- min(obs[x == xu[n1 - 1]])
          amat[n1 - 1, c1] <- 1
          amat[n1 - 1, c2] <- -1
          # return (amat)
        } else if (sh == 8) { ## decreasing concave
          amat <- -amat
          c1 <- min(obs[x == xu[1]])
          c2 <- min(obs[x == xu[2]])
          amat[n1 - 1, c1] <- 1
          amat[n1 - 1, c2] <- -1
          # return (amat)
        }
      }
      return(amat)
    } else if (all(block.ave == -1) & (length(block.ord) > 1)) {
      # nbcks is the total number of blocks
      block.ord_nz <- block.ord
      noord_ps <- which(block.ord == 0)

      # if(ncol(grid2) == 1){
      rm_id <- unique(sort(c(zeros_ps, noord_ps)))
      # } else {
      #  rm_id = unique(sort(noord_ps))
      # }
      if (length(rm_id) > 0) {
        block.ord_nz <- block.ord[-rm_id]
      }

      ubck <- unique(block.ord_nz)
      # use values in block.ord as an ordered integer vector
      ubck <- sort(ubck)

      nbcks <- length(table(block.ord_nz))
      szs <- as.vector(table(block.ord_nz))

      # new:
      if (is.list(xvs2)) {
        ps <- sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
      } else if (is.matrix(xvs2)) {
        # D_index=1
        # ps = lapply(xvs2, function(.x) which(as.vector(grid2) %in% .x)[1])
        # ps = as.matrix(grid2)
        # test! wrong to use x
        ps <- 1:M
      }

      # if(length(noord_ps) > 0){
      #  ps = ps[-noord_ps]
      # }
      if (length(rm_id) > 0) {
        ps <- ps[-rm_id]
      }
      # amat dimension: l1*l2...*lk by M
      amat <- NULL
      for (i in 1:(nbcks - 1)) {
        # bck_1 = bck_lst[[i]]; bck_2 = bck_lst[[i+1]]
        # l1 = length(bck_1); l2 = length(bck_2)
        ubck1 <- ubck[i]
        ubck2 <- ubck[i + 1]
        bck_1 <- which(block.ord_nz == ubck1)
        bck_2 <- which(block.ord_nz == ubck2)
        l1 <- length(bck_1)
        l2 <- length(bck_2)

        amat_bcki <- NULL
        rm_l <- length(zeros_ps) + length(noord_cells)
        M.ord <- length(block.ord)

        # amatj = matrix(0, nrow = l2, ncol = M.ord)
        amatj <- matrix(0, nrow = l2, ncol = M.ord - rm_l)
        bck_1k <- bck_1[1]
        # bck_1k = ps_bck1[1]
        amatj[, bck_1k] <- -1
        row_pointer <- 1
        for (j in 1:l2) {
          bck_2j <- bck_2[j]
          # bck_2j = ps_bck2[j]
          amatj[row_pointer, bck_2j] <- 1
          row_pointer <- row_pointer + 1
        }
        amatj0 <- amatj
        amat_bcki <- rbind(amat_bcki, amatj)

        if (l1 >= 2) {
          for (k in 2:l1) {
            bck_1k <- bck_1[k]
            bck_1k0 <- bck_1[k - 1]
            # bck_1k = ps_bck1[k]; bck_1k0 = ps_bck1[k-1]
            amatj <- amatj0
            # set the value for the kth element in block 1 to be -1
            # set the value for the (k-1)th element in block 1 back to 0
            # keep the values for block 2
            amatj[, bck_1k] <- -1
            amatj[, bck_1k0] <- 0
            amatj0 <- amatj
            amat_bcki <- rbind(amat_bcki, amatj)
          }
        }
        amat <- rbind(amat, amat_bcki)
      }

      if (length(noord_cells) > 0) {
        nr <- nrow(amat)
        nc <- ncol(amat)
        amat_tmp <- matrix(0, nrow = nr, ncol = (M.ord - rm_l + length(noord_cells)))
        if (length(zeros_ps) > 0) {
          ps <- which(block.ord[-zeros_ps] != 0)
        } else {
          ps <- which(block.ord != 0)
        }
        amat_tmp[, ps] <- amat
        amat <- amat_tmp
      }
      return(amat)
    } else if ((length(block.ave) > 1) & all(block.ord == -1)) {
      block.ave_nz <- block.ave
      # rm_id = unique(sort(c(zeros_ps, noord_cells)))

      noord_ps <- which(block.ave == 0)
      rm_id <- unique(sort(c(zeros_ps, noord_ps)))

      if (length(rm_id) > 0) {
        block.ave_nz <- block.ave[-rm_id]
      }

      ubck <- unique(block.ave_nz)
      # use values in block.ord as an ordered integer vector
      ubck <- sort(ubck)
      nbcks <- length(table(block.ave_nz))

      szs <- as.numeric(table(block.ave_nz))

      rm_l <- length(zeros_ps) + length(noord_ps)
      M.ave <- length(block.ave)
      amat <- matrix(0, nrow = (nbcks - 1), ncol = (M.ave - rm_l))
      for (i in 1:(nbcks - 1)) {
        ubck1 <- ubck[i]
        ubck2 <- ubck[i + 1]
        ps1 <- which(block.ave_nz == ubck1)
        ps2 <- which(block.ave_nz == ubck2)
        sz1 <- length(ps1)
        sz2 <- length(ps2)
        amat[i, ps1] <- -1 / sz1
        amat[i, ps2] <- 1 / sz2
      }

      if (length(noord_ps) > 0) {
        amat_tmp <- matrix(0, nrow = (nbcks - 1), ncol = (M.ave - rm_l + length(noord_ps)))
        if (length(zeros_ps) > 0) {
          ps <- which(block.ave[-zeros_ps] != 0)
        } else {
          ps <- which(block.ave != 0)
        }
        amat_tmp[, ps] <- amat
        amat <- amat_tmp
      }
      return(amat)
    } else if ((length(block.ave) > 1) & (length(block.ord) > 1)) {
      stop("only one of block average and block ordering can be used! \n")
    }
  }
  # return (amat)
}


######
#>=2D
######
.makeamat_2D <- function(x = NULL, sh, grid2 = NULL, xvs2 = NULL, zeros_ps = NULL, noord_cells = NULL,
                         Ds = NULL, suppre = FALSE, interp = FALSE,
                         relaxes = FALSE, block.ave.lst = NULL,
                         block.ord.lst = NULL, M_0 = NULL) {
  # n = length(x)
  # xu = sort(unique(x))
  # n1 = length(xu)
  sm <- 1e-7
  ms <- NULL
  Nd <- length(Ds)
  ne <- length(zeros_ps)
  if (Nd == 2) {
    if (any(sh > 0)) {
      D1 <- Ds[1]
      D2 <- Ds[2]
      obs <- 1:(D1 * D2)
      # new: group sh = 9 | 10 with other shapes
      if (sh[1] <= 10) {
        # new:
        if (sh[1] == 0) {
          amat1 <- NULL
        } else {
          if (length(zeros_ps) > 0) {
            amat0_lst <- list()
            grid3 <- grid2[-zeros_ps, , drop = FALSE]
            row_id <- as.numeric(rownames(grid3))
            cts <- row_id[which(row_id %% D1 == 0)]
            # temp: check if zeros_ps contains some point %% D1 = 0
            if (any((zeros_ps %% D1) == 0)) {
              zeros_ps_at_cts <- zeros_ps[which((zeros_ps %% D1) == 0)]
              cts_add <- row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
              cts <- unique(sort(c(cts, cts_add)))
            }
            obs <- 1:nrow(grid3)
            for (i in 1:length(cts)) {
              if (i == 1) {
                st <- 1
                # ed = cts[1]
                ed <- obs[row_id == cts[1]]
              } else {
                st <- obs[row_id == cts[i - 1]] + 1
                ed <- obs[row_id == cts[i]]
              }
              gridi <- grid3[st:ed, , drop = FALSE]
              na_ps <- which(apply(gridi, 1, function(x) all(is.na(x))))
              if (length(na_ps) > 0) {
                gridi <- gridi[-na_ps, , drop = FALSE]
              }
              obsi <- as.numeric(rownames(gridi))
              vec <- (D1 * (i - 1) + 1):(D1 * i)
              # vec = 1:D1
              zerosi <- vec[-which(vec %in% obsi)]
              non_zerosi <- vec[which(vec %in% obsi)]
              D1i <- nrow(gridi)
              # if ed == st, it could be: 1,2,3,4 with 0,0,0,92 observations
              # in this case, amat0i = NULL, and it will be removed before calling bdiag
              if (length(zerosi) >= 1 & (ed > st)) {
                if (sh[1] < 9) {
                  amat0i <- .makeamat(1:D1i, sh[1])
                } else {
                  # new:
                  zeros_ps_i <- zerosi %% D1
                  if (any(zeros_ps_i %in% 0)) {
                    zeros_ps_i[which(zeros_ps_i %in% 0)] <- D1
                  } else {
                    zeros_ps_i <- zerosi %% D1
                  }
                  # not sure....
                  amat0i <- .makeamat(1:D1i, sh[1],
                                      block.ave = block.ave.lst[[1]],
                                      block.ord = block.ord.lst[[1]],
                                      zeros_ps = zeros_ps_i,
                                      noord_cells = noord_cells[[1]],
                                      D_index = 1, xvs2 = xvs2, grid2 = grid2
                  )
                }
                amat0_lst[[i]] <- amat0i
              } else if (length(zerosi) >= 1 & (ed == st)) {
                # next
                amat0i <- 0
                amat0_lst[[i]] <- amat0i
              } else if (length(zerosi) == 0) {
                if (sh[1] < 9) {
                  amat0i <- .makeamat(1:D1, sh[1])
                } else {
                  # new:
                  amat0i <- .makeamat(1:D1i, sh[1],
                                      block.ave = block.ave.lst[[1]],
                                      block.ord = block.ord.lst[[1]],
                                      zeros_ps = NULL,
                                      noord_cells = noord_cells[[1]],
                                      D_index = 1,
                                      xvs2 = xvs2, grid2 = grid2
                  )
                }
                amat0_lst[[i]] <- amat0i
              }
              # amat0_lst[[i]] = amat0i
            }
          } else {
            # if(sh[1] < 9){
            amat0 <- .makeamat(1:D1,
                               sh = sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
                               zeros_ps = NULL, noord_cells = noord_cells[[1]], D_index = 1, xvs2 = xvs2, grid2 = grid2
            )
            # } else {
            #  amat0 = makeamat_2d_block(x=1:D1, sh=sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
            #        zeros_ps=NULL, noord_cells = noord_cells[[1]], D_index = 1, xvs2 = xvs2, grid2 = grid2, M_0 = M_0)
            # }
            amat0_lst <- rep(list(amat0), D2)
          }
          amat0_lst <- Filter(Negate(is.null), amat0_lst)
          amat1 <- bdiag(amat0_lst)
          amat1 <- as.matrix(amat1)
        }
      }
      # 2D
      amat2 <- NULL
      if (sh[2] < 3) {
        nr2 <- D1 * (D2 - 1)
        nc2 <- D1 * D2
        if (length(zeros_ps) == 0) {
          amat2 <- matrix(0, nrow = nr2, ncol = nc2)
          # new: if sh[2] == 0, keep the zero matrix
          if (sh[2] > 0) {
            for (i in 1:nr2) {
              amat2[i, i] <- -1
              amat2[i, (i + D1)] <- 1
            }
          }
        } else {
          # temp
          # amat2 = matrix(0, nrow=1, ncol=nc2)
          amat2_lst <- vector("list", length = nr2)
          rm_id <- NULL
          for (i in 1:nr2) {
            if (i %in% zeros_ps) {
              if (i <= D1) {
                amat2_lst[[i]] <- 1:nc2 * 0
                rm_id <- c(rm_id, i)
              } else if (i <= (nc2 - D1) & i > D1) {
                # ch_ps = c(ch_ps, i-D1)
                amat2_lst[[i]] <- 1:nc2 * 0
                rm_id <- c(rm_id, i)
                rowi <- 1:nc2 * 0
                rowi[i - D1] <- -1
                rowi[i + D1] <- 1
                amat2_lst[[i - D1]] <- rowi
              }
            } else {
              rowi <- 1:nc2 * 0
              rowi[i] <- -1
              rowi[i + D1] <- 1
              amat2_lst[[i]] <- rowi
              # amat2 = rbind(amat2, amat2i)
            }
          }
          for (i in (nr2 + 1):nc2) {
            if (i %in% zeros_ps) {
              amat2_lst[[i - D1]] <- 1:nc2 * 0
              rm_id <- c(rm_id, i - D1)
            }
          }
          if (!is.null(rm_id)) {
            amat2_lst <- amat2_lst[-rm_id]
          }
          amat2 <- NULL
          for (i in 1:length(amat2_lst)) {
            amat2 <- rbind(amat2, amat2_lst[[i]])
          }
          amat2 <- amat2[, -zeros_ps, drop = FALSE]
        }
        # new: if sh[2] == 0, keep the zero matrix
        if (sh[2] == 0) {
          # nr2 = nrow(amat2)
          # nc2 = ncol(amat2)
          # amat2 = matrix(0, nrow=nr2, ncol=nc2)
          amat2 <- NULL
        }
        if (sh[2] == 2) {
          amat2 <- -amat2
        }
      } else if (sh[2] == 3 | sh[2] == 4) {
        nr2 <- D1 * (D2 - 2)
        nc2 <- D1 * D2
        amat2 <- matrix(0, nrow = nr2, ncol = nc2)
        xu <- 1:D2
        for (i in 1:nr2) {
          # c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
          amat2[i, i] <- 1
          amat2[i, (i + D1)] <- -2
          amat2[i, (i + 2 * D1)] <- 1
        }
        if (sh[2] == 4) {
          amat2 <- -amat2
        }
      } else if (sh[2] > 4 & sh[2] < 11) {
        nr2 <- D1 * (D2 - 2)
        nc2 <- D1 * D2
        amat2 <- matrix(0, nrow = nr2, ncol = nc2)
        xu <- 1:D2
        for (i in 1:nr2) {
          # c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
          # amat2[i, c1] = 1; amat2[i, (c2+D1)] = -2; amat2[i, (c3+2*D1)] = 1
          amat2[i, i] <- 1
          amat2[i, (i + D1)] <- -2
          amat2[i, (i + 2 * D1)] <- 1
        }
        if (sh[2] == 5) { ### increasing convex
          amat2_add <- matrix(0, nrow = D1, ncol = nc2)
          # c1 = min(obs[xu == xu[1]])
          # c2 = min(obs[xu == xu[2]])
          # amat2[nr2, c1] = -1; amat2[nr2, (c1+D1)] = 1
          for (i in 1:D1) {
            amat2_add[i, i] <- -1
            amat2_add[i, i + D1] <- 1
          }
          amat2 <- rbind(amat2, amat2_add)
        } else if (sh[2] == 6) { ## decreasing convex
          # c1 = min(obs[xu == xu[D2]]); c2 = min(obs[xu == xu[D2 - 1]])
          amat2_add <- matrix(0, nrow = D1, ncol = nc2)
          for (i in 1:D1) {
            amat2_add[i, i + 2 * D1] <- -1
            amat2_add[i, i + D1] <- 1
          }
          amat2 <- rbind(amat2, amat2_add)
        } else if (sh[2] == 7) { ## increasing concave
          amat2 <- -amat2
          # c1 = min(obs[xu == xu[D2]]); c2 = min(obs[xu == xu[D2 - 1]])
          # amat2[nr2, c1+D1] = 1; amat2[nr2, c2] = -1
          amat2_add <- matrix(0, nrow = D1, ncol = nc2)
          for (i in 1:D1) {
            amat2_add[i, i + 2 * D1] <- 1
            amat2_add[i, i + D1] <- -1
          }
          amat2 <- rbind(amat2, amat2_add)
        } else if (sh[2] == 8) { ## decreasing concave
          amat2 <- -amat2
          # c1 = min(obs[xu == xu[1]]); c2 = min(obs[xu == xu[2]])
          # amat2[nr2, c1] = 1; amat2[nr2, (c2+D1)] = -1
          amat2_add <- matrix(0, nrow = D1, ncol = nc2)
          for (i in 1:D1) {
            amat2_add[i, i] <- 1
            amat2_add[i, i + D1] <- -1
          }
          amat2 <- rbind(amat2, amat2_add)
        } else if (sh[2] == 9 | sh[2] == 10) { ## block ordering or block average
          amat2 <- .makeamat_2d_block(
            x = NULL, sh[2], block.ave = block.ave.lst[[2]],
            block.ord = block.ord.lst[[2]],
            zeros_ps = zeros_ps,
            noord_cells = noord_cells[[2]],
            D_index = 2, xvs2 = xvs2, grid2 = grid2, M_0 = M_0
          )
        }
      }
      amat <- rbind(amat1, amat2)
    }
  } else if (Nd >= 3) {
    M <- length(Ds)
    D1 <- Ds[1]
    cump <- cumprod(Ds)
    nc <- cump[M]
    amat1 <- amat2 <- amatm <- NULL

    # 1D: group sh=9 or sh=10 with other shapes in 1D
    if (sh[1] == 0) {
      amat1 <- NULL
    } else {
      if (length(zeros_ps) > 0) {
        amat0_lst <- list()
        grid3 <- grid2[-zeros_ps, , drop = FALSE]
        row_id <- as.numeric(rownames(grid3))
        cts <- row_id[which(row_id %% D1 == 0)]
        if (any((zeros_ps %% D1) == 0)) {
          zeros_ps_at_cts <- zeros_ps[which((zeros_ps %% D1) == 0)]
          cts_add <- row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
          # cts = sort(c(cts, cts_add))
          cts <- unique(sort(c(cts, cts_add)))
        }
        obs <- 1:nrow(grid3)

        for (i in 1:length(cts)) {
          if (i == 1) {
            st <- 1
            ed <- obs[row_id == cts[1]]
          } else {
            st <- obs[row_id == cts[i - 1]] + 1
            ed <- obs[row_id == cts[i]]
          }
          gridi <- grid3[st:ed, , drop = FALSE]
          na_ps <- which(apply(gridi, 1, function(x) all(is.na(x))))
          if (length(na_ps) > 0) {
            gridi <- gridi[-na_ps, , drop = FALSE]
          }
          obsi <- as.numeric(rownames(gridi))
          vec <- (D1 * (i - 1) + 1):(D1 * i)
          zerosi <- vec[-which(vec %in% obsi)]
          non_zerosi <- vec[which(vec %in% obsi)]
          D1i <- nrow(gridi)
          if (length(zerosi) >= 1 & (ed > st)) {
            # amat0i = makeamat(1:D1i, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            # new:
            zeros_ps_i <- zerosi %% D1
            if (any(zeros_ps_i %in% 0)) {
              zeros_ps_i[which(zeros_ps_i %in% 0)] <- D1
            } # else {
            # zeros_ps_i = zerosi%%D1
            # }
            amat0i <- .makeamat(1:D1i, sh[1],
                                block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
                                zeros_ps = zeros_ps_i, noord_cells = noord_cells[[1]], xvs2 = xvs2, grid2 = grid2
            )
            amat0_lst[[i]] <- amat0i
          } else if (length(zerosi) >= 1 & (ed == st)) {
            amat0i <- 0
            amat0_lst[[i]] <- amat0i
          } else if (length(zerosi) == 0) {
            # amat0i = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            amat0i <- .makeamat(1:D1, sh[1],
                                block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
                                zeros_ps = zeros_ps, noord_cells = noord_cells[[1]], xvs2 = xvs2, grid2 = grid2
            )
            amat0_lst[[i]] <- amat0i
          }
        }
        amat0_lst <- Filter(Negate(is.null), amat0_lst)
        amat1 <- bdiag(amat0_lst)
        amat1 <- as.matrix(amat1)
      } else {
        n1 <- nc / D1
        amat1_0 <- .makeamat(1:D1, sh[1],
                             block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
                             zeros_ps = zeros_ps, noord_cells = noord_cells[[1]], xvs2 = xvs2, grid2 = grid2
        )
        amat1_0_lst <- rep(list(amat1_0), n1)
        amat1 <- bdiag(amat1_0_lst)
        amat1 <- as.matrix(amat1)
      }
    }
    amat <- amat1

    # 2D
    for (i in 2:(M - 1)) {
      Di <- Ds[i]
      cump <- cumprod(Ds[1:i])
      nci <- cump[i]
      gap <- cump[i - 1]
      ni <- nc / nci

      if (sh[i] %in% c(9, 10)) {
        amati <- .makeamat_2d_block(
          x = NULL, sh[i], block.ave = block.ave.lst[[i]],
          block.ord = block.ord.lst[[i]],
          zeros_ps = zeros_ps,
          noord_cells = noord_cells[[i]],
          D_index = i, xvs2 = xvs2, grid2 = grid2, M_0 = M_0
        )
      } else {
        # new: sh[i] == 0
        if (sh[i] == 0) {
          amati <- NULL
        }
        if (sh[i] > 0) {
          if ((sh[i] %in% c(1, 2))) {
            nr <- gap * (Di - 1)
          }
          if (sh[i] > 2) {
            if (Di < 3) {
              stop("For monotonicity + convexity constraints, number of domains must be >= 3! \n")
            } else {
              nr <- gap * (Di - 2)
            }
          }
          amati_0 <- .makeamat_2d(sh = sh[i], nr2 = nr, nc2 = nci, gap = gap, D2 = Di, Ds = Ds)

          if (length(zeros_ps) == 0) {
            amati_0_lst <- rep(list(amati_0), ni)
          } else {
            amati_0_lst <- vector("list", length = ni)
            cts <- seq(nci, nc, length = ni)
            bool_lst <- map(zeros_ps, .f = function(.x) as.numeric(cts >= .x))
            # ps_lst keeps track of the cut points such that the point giving -1 is the upp bd
            # of the interval containining one zero cell
            ps <- map_dbl(bool_lst, .f = function(.x) min(which(.x == 1)))
            cts_check <- cts[ps]
            for (k in 1:ni) {
              if (!(cts[k] %in% cts_check)) {
                # no zero cell in this interval
                amati_0_lst[[k]] <- amati_0
              } else {
                # find out the empty cells within an interval
                if (k == 1) {
                  zeros_ctsk <- zeros_ps[zeros_ps <= cts[k]]
                } else {
                  zeros_ctsk <- zeros_ps[zeros_ps <= cts[k] & zeros_ps > cts[k - 1]]
                }
                nek <- length(zeros_ctsk)
                rm_id <- NULL
                amati_0_cp <- amati_0
                nci <- ncol(amati_0)
                col_ps_vec <- NULL
                for (k2 in 1:nek) {
                  zero_k2 <- zeros_ctsk[k2]
                  # check: zero_k2 %% nr = 0
                  if (zero_k2 != cts[k]) {
                    col_ps <- zero_k2 %% nci
                  } else {
                    col_ps <- nci
                  }
                  col_ps_vec <- c(col_ps_vec, col_ps)

                  row_ps <- col_ps - gap
                  # if((zero_k2 %% nr + gap) > nr){
                  if ((col_ps + gap) > nci) {
                    rm_id <- c(rm_id, row_ps)
                  } else {
                    rm_id <- c(rm_id, col_ps)
                    # test:
                    if (row_ps > 0) {
                      amati_0_cp[row_ps, (col_ps + gap)] <- 1
                      amati_0_cp[row_ps, col_ps] <- 0
                    }
                  }
                }
                amati_0_cp <- amati_0_cp[-rm_id, ]
                # amati_0_cp = amati_0_cp[,-(zeros_ctsk%%nr)]
                amati_0_cp <- amati_0_cp[, -col_ps_vec]
                amati_0_lst[[k]] <- amati_0_cp
              }
            }
          }
          amati <- bdiag(amati_0_lst)
          amati <- as.matrix(amati)
        }
      }
      amat <- rbind(amat, amati)
    }

    # 3D
    cump <- cumprod(Ds[1:(M - 1)])
    gap <- cump[M - 1]
    Dm <- Ds[M]
    if ((sh[M] %in% c(9, 10))) {
      amatm <- .makeamat_2d_block(
        x = NULL, sh[M], block.ave = block.ave.lst[[M]],
        block.ord = block.ord.lst[[M]],
        zeros_ps = zeros_ps,
        noord_cells = noord_cells[[M]],
        D_index = M, xvs2 = xvs2, grid2 = grid2, M_0 = M_0
      )
    } else {
      if ((sh[M] %in% c(1, 2))) {
        nr <- gap * (Dm - 1)
      }
      if (sh[M] >= 3) {
        nr <- gap * (Dm - 2)
      }
      # new: sh[M] == 0
      if (sh[M] == 0) {
        # nrm = nrow(amatm)
        # ncm = ncol(amatm)
        # amatm = matrix(0, nrow=nrm, ncol=ncm)
        amatm <- NULL
      } else {
        amatm <- .makeamat_2d(sh = sh[M], nr2 = nr, nc2 = nc, gap = gap, D2 = Dm, Ds = Ds)
        ncm <- ncol(amatm)
        if (length(zeros_ps) > 0) {
          amatm_cp <- amatm
          rm_id <- NULL
          for (k in 1:ne) {
            zero_k <- zeros_ps[k]
            col_ps <- zero_k
            row_ps <- col_ps - gap
            if ((col_ps + gap) > ncm) {
              rm_id <- c(rm_id, row_ps)
            } else {
              rm_id <- c(rm_id, col_ps)
              # test:
              if (row_ps > 0) {
                amatm_cp[row_ps, (col_ps + gap)] <- 1
                amatm_cp[row_ps, col_ps] <- 0
              }
            }
          }
          amatm_cp <- amatm_cp[-rm_id, ]
          amatm_cp <- amatm_cp[, -zeros_ps]
          amatm <- amatm_cp
        }
      }
    }
    amat <- rbind(amat, amatm)
  }
  return(amat)
}

############################################
# local function to be called in makeamat_2D
############################################
.makeamat_2d <- function(x = NULL, sh, nr2, nc2, gap, D1 = NULL, D2, Ds = NULL, suppre = FALSE, interp = FALSE) {
  sm <- 1e-7
  ms <- NULL
  Nd <- length(Ds)

  amat2 <- NULL
  if (sh < 3) {
    amat2 <- matrix(0, nrow = nr2, ncol = nc2)
    for (i in 1:nr2) {
      amat2[i, i] <- -1
      amat2[i, (i + gap)] <- 1
    }
    if (sh == 2) {
      amat2 <- -amat2
    }
    return(amat2)
  } else if (sh == 3 | sh == 4) {
    amat2 <- matrix(0, nrow = nr2, ncol = nc2)
    # xu = 1:D2
    for (i in 1:nr2) {
      # c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
      # amat2[i, c1] = 1; amat2[i, (c2+gap)] = -2; amat2[i, (c3+2*gap)] = 1
      amat2[i, i] <- 1
      amat2[i, (i + gap)] <- -2
      amat2[i, (i + 2 * gap)] <- 1
    }
    if (sh == 4) {
      amat2 <- -amat2
    }
    return(amat2)
  } else if (sh > 4 & sh < 9) {
    amat2 <- matrix(0, nrow = nr2, ncol = nc2)
    # xu = 1:D2
    for (i in 1:nr2) {
      # c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
      # amat2[i, c1] = 1; amat2[i, (c2+gap)] = -2; amat2[i, (c3+2*gap)] = 1
      amat2[i, i] <- 1
      amat2[i, (i + gap)] <- -2
      amat2[i, (i + 2 * gap)] <- 1
    }
    if (sh == 5) { ### increasing convex
      nr2_add <- gap
      amat2_add <- matrix(0, nrow = nr2_add, ncol = nc2)
      for (i in 1:nr2_add) {
        amat2_add[i, i] <- -1
        amat2_add[i, i + gap] <- 1
      }
      amat2 <- rbind(amat2, amat2_add)
      return(amat2)
    } else if (sh == 6) { ## decreasing convex
      nr2_add <- gap
      amat2_add <- matrix(0, nrow = nr2_add, ncol = nc2)
      for (i in 1:nr2_add) {
        amat2_add[i, i + 2 * gap] <- -1
        amat2_add[i, i + gap] <- 1
      }
      amat2 <- rbind(amat2, amat2_add)
      return(amat2)
    } else if (sh == 7) { ## increasing concave
      nr2_add <- gap
      amat2 <- -amat2
      amat2_add <- matrix(0, nrow = nr2_add, ncol = nc2)
      for (i in 1:nr2_add) {
        amat2_add[i, i + 2 * gap] <- 1
        amat2_add[i, i + gap] <- -1
      }
      amat2 <- rbind(amat2, amat2_add)
      return(amat2)
    } else if (sh == 8) { ## decreasing concave
      nr2_add <- gap
      amat2 <- -amat2
      amat2_add <- matrix(0, nrow = nr2_add, ncol = nc2)
      for (i in 1:nr2_add) {
        amat2_add[i, i] <- 1
        amat2_add[i, i + gap] <- -1
      }
      amat2 <- rbind(amat2, amat2_add)
      return(amat2)
    }
  }
  # return (amat2)
  # rslt = list(amat = amat,  ms = ms)
  # rslt
}


################################################################
# local function to be called in makeamat_2D for block ordering
################################################################
.makeamat_2d_block <- function(x = NULL, sh, block.ave = NULL, block.ord = NULL, zeros_ps = NULL,
                               noord_cells = NULL, D_index = 1, xvs2 = NULL,
                               grid2 = NULL, M_0 = NULL) {
  # block average is not added yet
  if (all(block.ave == -1) & (length(block.ord) > 1)) {
    if (length(zeros_ps) > 0) {
      iter <- 0
      ps <- sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
      ps0 <- ps
      ps1_ed <- ps[2] - 1

      amat <- NULL
      # obs = 1:length(block.ord)
      while (iter < ps1_ed) {
        ps <- ps0 + iter

        block.ord_nz <- block.ord
        noord_ps <- which(block.ord == 0)
        rm_id <- noord_ps

        # test:
        z_rm_id <- NULL
        if (length(zeros_ps) > 0) {
          for (e in zeros_ps) {
            z_rm_id <- c(z_rm_id, which(ps == e))
          }
        }

        if (length(z_rm_id) > 1) {
          rm_id <- sort(c(z_rm_id, noord_ps))
        }

        if (length(rm_id) > 0) {
          block.ord_nz <- block.ord[-rm_id]
        }

        ubck <- unique(block.ord_nz)
        # use values in block.ord as an ordered integer vector
        ubck <- sort(ubck)

        nbcks <- length(table(block.ord_nz))
        szs <- as.vector(table(block.ord_nz))

        if (length(rm_id) > 0) {
          ps <- ps[-rm_id]
        }

        # amati dimension: l1*l2...*lk by M
        amati <- NULL
        for (i in 1:(nbcks - 1)) {
          ubck1 <- ubck[i]
          ubck2 <- ubck[i + 1]
          bck_1 <- which(block.ord_nz == ubck1)
          bck_2 <- which(block.ord_nz == ubck2)
          l1 <- length(bck_1)
          l2 <- length(bck_2)

          ps_bck1 <- ps[which(block.ord_nz == ubck1)]
          ps_bck2 <- ps[which(block.ord_nz == ubck2)]
          # for each element in block1, the value for each element in block 2 will be the same
          # for each element in block1, the number of comparisons is l2
          amat_bcki <- NULL

          # tmp:
          if (D_index == 1) {
            M.ord <- length(block.ord)
          } else {
            M.ord <- M_0
            # M.ord = nrow(grid2)
          }

          amatj <- matrix(0, nrow = l2, ncol = M.ord)
          # amatj = matrix(0, nrow = l2, ncol = M.ord-rm_l)
          # bck_1k = bck_1[1]
          bck_1k <- ps_bck1[1]
          amatj[, bck_1k] <- -1
          row_pointer <- 1
          for (j in 1:l2) {
            # bck_2j = bck_2[j]
            bck_2j <- ps_bck2[j]
            amatj[row_pointer, bck_2j] <- 1
            row_pointer <- row_pointer + 1
          }
          amatj0 <- amatj
          amat_bcki <- rbind(amat_bcki, amatj)

          if (l1 >= 2) {
            for (k in 2:l1) {
              # bck_1k = bck_1[k]; bck_1k0 = bck_1[k-1]
              bck_1k <- ps_bck1[k]
              bck_1k0 <- ps_bck1[k - 1]
              amatj <- amatj0
              # set the value for the kth element in block 1 to be -1
              # set the value for the (k-1)th element in block 1 back to 0
              # keep the values for block 2
              amatj[, bck_1k] <- -1
              amatj[, bck_1k0] <- 0
              amatj0 <- amatj
              amat_bcki <- rbind(amat_bcki, amatj)
            }
          }
          amati <- rbind(amati, amat_bcki)
        }

        amat <- rbind(amat, amati)
        iter <- iter + 1
      }
      amat <- amat[, -zeros_ps]
      return(amat)
    } else {
      ps <- sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
      # nbcks is the total number of blocks
      block.ord_nz <- block.ord
      noord_ps <- which(block.ord == 0)
      rm_id <- unique(sort(noord_ps))
      # test:
      # z_rm_id = NULL
      # if(length(zeros_ps) > 0){
      #  for(e in zeros_ps){
      #    z_rm_id = c(z_rm_id, max(which(ps <= e)))
      #  }
      # }
      # if(ncol(grid2) == 1){
      # rm_id = unique(sort(c(z_rm_id, noord_ps)))
      # } else {
      #  rm_id = unique(sort(noord_ps))
      # }
      if (length(rm_id) > 0) {
        block.ord_nz <- block.ord[-rm_id]
      }

      ubck <- unique(block.ord_nz)
      # use values in block.ord as an ordered integer vector
      ubck <- sort(ubck)

      nbcks <- length(table(block.ord_nz))
      szs <- as.vector(table(block.ord_nz))

      if (length(rm_id) > 0) {
        ps <- ps[-rm_id]
      }
      # amat dimension: l1*l2...*lk by M
      amat <- NULL
      for (i in 1:(nbcks - 1)) {
        ubck1 <- ubck[i]
        ubck2 <- ubck[i + 1]
        bck_1 <- which(block.ord_nz == ubck1)
        bck_2 <- which(block.ord_nz == ubck2)
        l1 <- length(bck_1)
        l2 <- length(bck_2)

        ps_bck1 <- ps[which(block.ord_nz == ubck1)]
        ps_bck2 <- ps[which(block.ord_nz == ubck2)]
        # for each element in block1, the value for each element in block 2 will be the same
        # for each element in block1, the number of comparisons is l2
        amat_bcki <- NULL

        # tmp:
        if (D_index == 1) {
          M.ord <- length(block.ord)
        } else {
          M.ord <- M_0
          # M.ord = nrow(grid2)
        }

        amatj <- matrix(0, nrow = l2, ncol = M.ord)
        # amatj = matrix(0, nrow = l2, ncol = M.ord-rm_l)
        # bck_1k = bck_1[1]
        bck_1k <- ps_bck1[1]
        amatj[, bck_1k] <- -1
        row_pointer <- 1
        for (j in 1:l2) {
          # bck_2j = bck_2[j]
          bck_2j <- ps_bck2[j]
          amatj[row_pointer, bck_2j] <- 1
          row_pointer <- row_pointer + 1
        }
        amatj0 <- amatj
        amat_bcki <- rbind(amat_bcki, amatj)

        if (l1 >= 2) {
          for (k in 2:l1) {
            # bck_1k = bck_1[k]; bck_1k0 = bck_1[k-1]
            bck_1k <- ps_bck1[k]
            bck_1k0 <- ps_bck1[k - 1]
            amatj <- amatj0
            # set the value for the kth element in block 1 to be -1
            # set the value for the (k-1)th element in block 1 back to 0
            # keep the values for block 2
            amatj[, bck_1k] <- -1
            amatj[, bck_1k0] <- 0
            amatj0 <- amatj
            amat_bcki <- rbind(amat_bcki, amatj)
          }
        }
        amat <- rbind(amat, amat_bcki)
      }

      if (D_index > 1) {
        iter <- 1
        amat0 <- amat
        nr <- nrow(amat0)
        nc <- ncol(amat0)
        ps <- sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
        ps1_ed <- ps[2] - 1

        if (length(noord_cells) > 0) {
          ps <- ps[-noord_cells]
        }

        while (iter < ps1_ed) {
          ps1 <- ps + iter
          amati <- matrix(0, nrow = nr, ncol = nc)
          amati[, ps1] <- amat0[, ps]
          amat <- rbind(amat, amati)
          iter <- iter + 1
        }
      }
    }
    return(amat)
  }
}
