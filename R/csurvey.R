csvy <- function(formula, data, design, nD=NULL, family=gaussian, amat=NULL, 
                level=0.95, n.mix=100L, test=TRUE) {
  cl <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("'family' not recognized!")
  labels <- NULL
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  shapes_add <- NULL
  xmat_add0 <- NULL; xmat_add <- NULL; xnms_add <- NULL
  nums_add <- NULL; xid_add <- 1
  relaxes <- NULL
  block.ave.lst <- vector("list", length = (ncol(mf)-1))
  block.ord.lst <- vector("list", length = (ncol(mf)-1))
  cnms <- colnames(mf)
  for (i in 2:ncol(mf)) {
    #if (is.numeric(attributes(mf[,i])$shape) && (attributes(mf[,i])$categ == "additive")) {
    if (is.numeric(attributes(mf[,i])$shape)) {
      labels <- c(labels, "additive")
      shapes_add <- c(shapes_add, attributes(mf[,i])$shape)
      reli <- attributes(mf[,i])$relax
      if (is.null(reli)) {
        reli <- FALSE
      } else {
        reli <- TRUE
      }
      relaxes <- c(relaxes, reli)
      if (is.factor(mf[,i]) && is.character(levels(mf[,i]))) {
        #some factor var will have more levels than observed groups
        xmat_add <- cbind(xmat_add, as.numeric(mf[,i]))
      } else {
        xmat_add <- cbind(xmat_add, mf[,i])
      }
      xmat_add0 <- cbind(xmat_add0, mf[,i])
      xnms_add <- c(xnms_add, attributes(mf[,i])$nm)
      xid_add <- xid_add + 1
      
      block.ord.lst[[i-1]] <- -1
      block.ave.lst[[i-1]] <- -1
      #print(i-1)
      if ((attributes(mf[,i])$categ == "block.ord")) {
        block.ord.lst[[i-1]] <- attributes(mf[,i])$order
        #xid_ord <- xid_ord + 1
      }
      if ((attributes(mf[,i])$categ == "block.ave")) {
        block.ave.lst[[i-1]] <- attributes(mf[,i])$order
        #xid_ave <- xid_ave + 1
      }
    } else if (is.null(attributes(mf[,i])$shape)) {
      #new: if shape is NULL, then code it as 0
      block.ord.lst[[i-1]] <- -1
      block.ave.lst[[i-1]] <- -1
      reli <- attributes(mf[,i])$relax
      if (is.null(reli)) {
        reli <- FALSE
      } else {
        reli <- TRUE
      }
      relaxes <- c(relaxes, reli)
      shapes_add <- c(shapes_add, 0)
      labels <- c(labels, "additive")
      xmat_add <- cbind(xmat_add, mf[,i])
      xmat_add0 <- cbind(xmat_add0, mf[,i])
      #xnms_add <- c(xnms_add, attributes(mf[,i])$nm)
      xnms_add <- c(xnms_add, cnms[i])
      xid_add <- xid_add + 1
    }
  }
  fit <- eval(call("csvy.fit", design = design, M = nD, relaxes = relaxes, xm = xmat_add, 
                   sh = shapes_add, ynm = ynm, amat = amat, 
                   block.ave.lst = block.ave.lst, block.ord.lst = block.ord.lst, 
                   level = level, n.mix = n.mix, cl = cl, test = test))
  rslt <- list(design = design, data = data, muhat = fit$muhat, muhat.s = fit$muhat.s, muhat.un = fit$muhat.un, lwr = fit$lwr, upp = fit$upp, 
               lwr.s = fit$lwr.s, upp.s = fit$upp.s, lwru = fit$lwru, uppu = fit$uppu, domain = fit$domain, domain.ord = fit$domain.ord, amat = fit$amat, 
               grid = fit$grid, grid_ps = fit$grid_ps, nd = fit$nd, Ds = fit$Ds, cov.c = fit$cov.c, cov.un = fit$cov.un, ec = fit$ec, xmat_add = xmat_add, xmat_add0 = xmat_add0, xnms_add = xnms_add, 
               shapes_add = shapes_add, ynm = ynm, zeros_ps = fit$zeros_ps, pval = fit$pval, CIC = fit$CIC)
  rslt$call <- cl
  class(rslt) <- "csvy"
  return (rslt)
}

############
## csvy.fit#
############
csvy.fit = function(design, M=NULL, relaxes=NULL, xm=NULL, sh=1, ynm=NULL, nd=NULL, 
                    ID=NULL, block.ave.lst=NULL, block.ord.lst=NULL, amat=NULL, 
                    additive=TRUE, level=0.95, n.mix=100L, cl=NULL, h_grid=NULL, test=TRUE,...){
  #make sure each column in xm is numeric:
  #xm = map_dbl(xm, .f = function(.x) as.numeric(.x))
  bool = apply(xm, 2, function(x) is.numeric(x))
  if(any(!bool)){
    xm = apply(xm, 2, function(x) as.numeric(x))
  }
  Ds = apply(xm, 2, function(x) length(unique(x)))
  xm.red = unique(xm)
  #xvs2 returns labeled domain names from 1 to max x
  #xvs2 = apply(xm.red, 2, function(x) 1:length(unique(x)))
  xvs2 = apply(xm.red, 2, function(x) sort(unique(x)))
  if(is.matrix(xvs2)){
    grid2 = expand.grid(as.data.frame(xvs2))
  }else if(is.list(xvs2)){
    grid2 = expand.grid(xvs2)
  }
  df = design$variables
  ds = design
  #sample sizes for each observed domain
  if(is.null(nd)){
    xm.df = as.data.frame(xm)
    nd = aggregate(list(nd = rep(1, nrow(xm.df))), xm.df, length)$nd
  }
  if(is.null(ID)){
    #ID = apply(xm, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
    #new: for empty cell
    if (length(Ds) == 1) {
      ID = apply(xm, 1, function(elem) grid2[which(apply(grid2, 1, function(gi) all(gi == elem))),])
    }
    if (length(Ds) > 1) {
      ID = apply(xm, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
    }
  }
  ID = as.ordered(ID)
  
  uds = unique(ID)
  domain.ord = order(uds)
  uds = sort(uds)
  #new M:
  obs_cells = as.numeric(levels(uds))[uds]
  #let user to define M?
  if (is.null(M)) {
    #some x won't start from 1:
    #M = length(levels(uds))#wrong: cannot include empty cells
    M = max(obs_cells)
  }
  if (any(class(ds) == "survey.design")) {
    weights = 1/(ds$prob)
  }
  if (any(class(ds) == "svyrep.design")) {
    weights = ds$pweights
  }
  #Nhat and stratawt will follow the order: 1 to M
  y = df[[ynm]]
  n = length(y)
  ys = stratawt = vector("list", length=M)
  Nhat = nd = 1:M*0
  #for (i in 1:M) {
  for (udi in uds){
    #for(udi in obs_cells){
    udi = as.numeric(udi)
    ps = which(ID %in% udi)
    ysi = y[ps]
    wi = weights[ps]
    Nhi = sum(wi)
    
    Nhat[udi] = Nhi
    nd[udi] = length(wi)
    ys[[udi]] = ysi
    stratawt[[udi]] = wi
  }
  if (any(class(design) == "survey.design")) {
    rds = as.svrepdesign(ds, type="JKn")
  }
  if (any(class(ds) == "svyrep.design")) {
    rds = ds
  }
  fo = as.formula(paste0("~", ynm))
  ans.unc = svyby(formula=fo, by=~ID, design=rds, FUN=svymean, covmat=TRUE, multicore=FALSE, drop.empty.groups=FALSE)
  v1 = vcov(ans.unc)
  yvecu = yvec = ans.unc[[ynm]]
  vsu = diag(v1)
  
  #new: replace n=1 variance with the average
  #not to change vsu, which will be used for 
  ps = which(round(diag(v1), 7) == 0L)
  diag(v1)[ps] = mean(diag(v1)[-ps])
  
  w = Nhat/sum(Nhat)
  #need to add incr + decr; decr + decr; 2 decr
  #Nd = length(sh)
  #amat and w will follow the order
  Nd = length(Ds)
  
  empty_cells = which(nd == 0)
  #temp
  #block.ave is not considered now
  #noord_cells will be a list, each elememt may have different length
  noord_cells = NULL
  #if(length(block.ave.lst) > 0 & Nd == 1){
  bool_ave = map_dbl(block.ave.lst, .f = function(.x) all(.x == -1))
  bool_ord = map_dbl(block.ord.lst, .f = function(.x) all(.x == -1))
  if(any(bool_ave == 0)){
    noord_cells = map(block.ave.lst, .f = function(.x) which(.x == 0))
    if(length(Ds) == 1){
      noord_cells = noord_cells[[1]]
    }
  }
  #if(length(block.ord.lst) > 0 & Nd == 1){
  if(any(bool_ord == 0)){
    noord_cells = map(block.ord.lst, .f = function(.x) which(.x == 0))
    if(length(Ds) == 1){
      noord_cells = noord_cells[[1]]
    }
  }
  ne = length(empty_cells)
  nd_ne = nd[empty_cells]
  #if (length(block.ave.lst) == 0 & length(block.ord.lst) == 0) {
  #new: impute for n=1 cells 
  #new: re-compute the variance for n=2...10 cells
  small_cells = which(nd >= 1 & nd <= 10)
  nsm = length(small_cells)
  #new: leave nd = 1 alone....
  #empty_cells = which(nd == 0)
  #observed sample size in cells with 0 or 1 n
  #} else {ne = 0; nsm = 0}
  
  #zeros_ps = empty_cells
  #zeros = ne
  zeros = ones = 0
  zeros_ps = NULL
  
  if (ne >= 1) {
    if (any(nd == 0)) {
      zeros_ps = which(nd == 0)
      zeros = length(zeros_ps)
    }
  }
  
  pskeep = ps = NULL
  M_0 = M
  w_0 = w
  
  M = M - zeros
  imat = diag(M)
  if (!is.null(zeros_ps)) {
    w = w[-zeros_ps]
  }
  #define weight mat and its inverse here to avoid repeating it
  #need to consider zeros
  wt = diag(w)
  wtinv = diag(1/w)
  hkeep = NULL
  if (is.null(amat)) {
    if (Nd == 1) {
      #sh = 0 means user must define amat 
      if (sh > 0) {
        #Ds is not correct when >= 1 empty cell
        #new: try to find the constrained est. first
        block.ave = NULL
        block.ord = NULL
        if (length(block.ave.lst) > 0) {
          #print ('a')
          block.ave = block.ave.lst[[1]]
        }
        if (length(block.ord.lst) > 0) {
          #print ('a')
          block.ord = block.ord.lst[[1]]
        }
           
        #relaxed monotone is not fully debugged
        if (relaxes) {
          #amats = makeamat(1:M, sh=sh, Ds=M, interp=TRUE, relaxes=relaxes)
          if (is.null(h_grid)) {
            n_h = 60 
            h_grid = seq(0.01, 6, length = n_h)
          } else {
            n_h = length(h_grid)
          }
          amats = vector("list", length = n_h)
          for (i in seq_len(n_h)) {
            h = h_grid[i]
            amat = makeamat(1:M, sh=sh, Ds=M, interp=TRUE, relaxes=relaxes, h=h)
            amats[[i]] = amat
          }
          #CIC to choose h
          #find the h_hat besed on CIC criterion for each simulated y
          B1 = 10 #??
          #ysim = MASS::mvrnorm(B1, mu = yvec, Sigma = v1)
          ysim = MASS::mvrnorm(B1, mu = muhat, Sigma = v1)
          index = rep(0, B1)
          #wtinv = diag(1/w)
          for (j in 1:B1) {
            CIC_h = rep(0, n_h)
            #rme = matrix(0, M, n_h)
            for (i in seq_len(n_h)) {
              amati = amats[[i]]
              ansi = coneA(ysim[j, ], amati, w, msg=FALSE)
              #rme[,i] = ansi$thetahat
              facei = ansi$face
              
              if (length(facei) == 0){
                #pmat_is = diag(rep(0, M))
                mat = wt %*% v1
              } else {
                dd = amati[facei, ,drop=FALSE]
                #check solve
                wtinvdd = wtinv %*% t(dd)
                pmat_is = wtinvdd %*% solve(dd %*% wtinvdd) %*% dd
                mat = wt %*% (imat - pmat_is) %*% v1
              }
              evec = ysim[j,] - ansi$thetahat
              CIC_h[i] = t(evec) %*% wt %*% evec + 2*(sum(diag(mat)))
            }
            CIC_h = round(CIC_h, 6)
            index[j] = which.min(CIC_h)[1]
          }
          #ps = max(index)
          #try the 4th largest
          ps = (sort(index, decreasing = TRUE))[4]
          #choose amat_h 
          amat = amats[[ps]]
          hkeep = h_grid[ps]
        } else {
          #amat = makeamat(1:M, sh=sh, Ds=M, interp=TRUE, relaxes=FALSE)
          if(length(zeros_ps) > 0){
            M2 = seq(1, M+zeros)[-zeros_ps]
          }else{
            M2 = 1:M
          }
          amat = makeamat(x=M2, sh=sh, Ds=M, interp=TRUE, relaxes=FALSE, 
                          block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]], 
                          zeros_ps=zeros_ps, noord_cells=noord_cells)
          #amat = makeamat(1:M, sh=sh, Ds=M, interp=TRUE, relaxes=FALSE, block.ave=block.ave, block.ord=block.ord)
        }  
      } else if (sh == 0 & is.null(amat)) {
        stop ("User must define a constraint matrix!")
      }
    } else if (Nd >= 2) {
      #temp:
      if (any(sh > 0)) {
        #if (!relaxes) {
        amat = makeamat_2D(x=NULL, sh=sh, grid2=grid2, 
                           zeros_ps=zeros_ps, noord_cells=noord_cells, 
                           Ds=Ds, interp=TRUE, 
                           block.ave.lst=block.ave.lst,
                           block.ord.lst=block.ord.lst)
      } else if (all(sh == 0) & is.null(amat)) {
        stop ("User must define a constraint matrix!")
      }
    }
  }
  
  #if (Nd >= 2 || (Nd == 1 & sh == 0 & !is.null(amat))) {
  if (!is.null(amat)) {
    ans.polar = coneA(yvec, amat, w=w, msg=FALSE)
    muhat = round(ans.polar$thetahat, 10)
    face = ans.polar$face
    if (length(face) == 0){
      mat = wt %*% v1
    } else {
      dd = amat[face, ,drop=FALSE]
      wtinvdd = wtinv %*% t(dd)
      pmat_is = wtinvdd %*% solve(dd %*% wtinvdd) %*% dd
      mat = wt %*% (imat - pmat_is) %*% v1
    }
    evec = yvec - muhat
    CIC = t(evec) %*% wt %*% evec + 2*(sum(diag(mat)))
    CIC = as.numeric(CIC)
  }
  
  lwr = upp = NULL
  acov = v1
  dp = -t(amat)
  dp = apply(dp, 2, function(e) e / (sum(e^2))^(.5))
  # M is revised if there's empty cell
  if(n.mix >= 1){
    #M = M-zeros
    m_acc = M 
    sector = NULL
    times = NULL
    df.face = NULL
    iter = 1
    #vms = list()
    obs = 1:M
    
    #not finished: v1_tmp
    ysims = MASS::mvrnorm(n.mix, mu=muhat, Sigma=v1)
    for (iloop in 1:n.mix) {
      #ysim is the group mean
      #new:
      ysim = ysims[iloop, ]
      ansi = coneA(ysim, amat, w=w, msg=FALSE)
      muhati = round(ansi$thetahat, 10)
      
      facei = ansi$face
      if (length(facei) == 0) {
        next
      } else {
        sec = 1:nrow(amat)*0
        sec[facei] = 1
        
        r = makebin(sec) + 1
        if (iter == 1) {
          df.face = rbind(df.face, c(r, 1))
          sector = rbind(sector, sec)
        } else {
          if (r %in% df.face[,1]) {
            ps = which(df.face[,1] %in% r)
            df.face[ps,2] = df.face[ps,2] + 1
          } else {
            df.face = rbind(df.face, c(r, 1))
            sector = rbind(sector, sec)
          }
        }
        iter = iter+1
      }
    }
    
    if (!is.null(df.face)) {
      sm_id = which((df.face[,2]/n.mix) < 1e-3)
      if (any(sm_id)) {
        df.face = df.face[-sm_id, ,drop=FALSE]
        sector = sector[-sm_id, ,drop=FALSE]
      }
      nsec = nrow(df.face)
      bsec = df.face
      bsec[,2] = bsec[,2] / sum(bsec[,2])
      
      acov = matrix(0, nrow=M, ncol=M)
      
      for(is in 1:nsec) {
        jvec = sector[is, ]
        smat = dp[,which(jvec==1),drop=FALSE]
        wtinvs = wtinv %*% smat
        pmat_is_p = wtinv %*% smat %*% solve(t(smat) %*% wtinv %*% smat) %*% t(smat)
        pmat_is = (imat-pmat_is_p)
        acov = acov + bsec[is,2]*pmat_is%*%v1%*%t(pmat_is)
      }
    } else {
      acov = v1
    }
    z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
    #z.mult = 2
    vsc = diag(acov)
    hl = z.mult*sqrt(vsc)
    lwr = muhat - hl
    upp = muhat + hl
  }
  #new: if n.mix = 0, then use the faces for the final projection of the real y
  if (n.mix == 0L) {
    vmat = qr.Q(qr(t(amat)), complete = TRUE)[, -(1:(qr(t(amat))$rank)), drop = FALSE]
    pmat_is = vmat %*% solve(crossprod(vmat), t(vmat))
    acov = wtinv^(.5)%*%pmat_is%*%wt^(.5)%*%v1%*%wt^(.5)%*%pmat_is%*%wtinv^(.5)
    z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
    #z.mult = 2
    vsc = diag(acov)
    hl = z.mult*sqrt(vsc)
    lwr = muhat - hl
    upp = muhat + hl
  }
  muhatu = yvecu
  z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
  #z.mult = 2
  hl = z.mult*sqrt(vsu)
  lwru = muhatu - hl
  uppu = muhatu + hl
  
  #new: if there's empty cells, just augment muhatu to include NA's
  muhatu_all = lwru_all = uppu_all = 1:(M+zeros)*0
  muhatu_all[zeros_ps] = NA
  muhatu_all[obs_cells] = muhatu
  muhatu = muhatu_all
  
  lwru_all[zeros_ps] = NA
  lwru_all[obs_cells] = lwru
  lwru = lwru_all
  
  uppu_all[zeros_ps] = NA
  uppu_all[obs_cells] = uppu
  uppu = uppu_all
  
  #for monotone only: to find the c.i. for empty/1 cells by using smallest L and largest U of neighbors
  #new: create grid2 again because the routines to do imputation need to have each x start from 1
  xvs2 = apply(xm.red, 2, function(x) 1:length(unique(x)))
  if(is.matrix(xvs2)){
    grid2 = expand.grid(as.data.frame(xvs2))
  }else if(is.list(xvs2)){
    grid2 = expand.grid(xvs2)
  }
  if (ne >= 1) {
    ans_im = impute_em2(empty_cells, obs_cells, M=(M+zeros), yvec, Nd, nd, Nhat, w, 
                        domain.ord, grid2, Ds, sh, muhat, lwr, upp, vsc)
    muhat = ans_im$muhat
    lwr = ans_im$lwr
    upp = ans_im$upp
    vsc = ans_im$vsc
    domain.ord = ans_im$domain.ord
  }

  sig1 = vsc
  sig2 = NULL
  if (Nd >= 1 & nsm >= 1) {
    new_obs_cells = sort(unique(c(empty_cells, obs_cells)))
    if(Nd == 1){
      ans_im2 = impute_em2(small_cells, new_obs_cells, M=(M+zeros), yvec, 
                           Nd, nd, Nhat, w, domain.ord, grid2, Ds, sh, muhat, lwr, upp, vsc)
    }
    if(Nd >= 2){
      ans_im2 = impute_em3(small_cells, M=(M+zeros), yvec, Nd, nd, Nhat, w, domain.ord,
                           grid2, Ds, sh, muhat, lwr, upp, vsc, new_obs_cells)
    }
    lwr = ans_im2$lwr
    upp = ans_im2$upp
    sig2 = ans_im2$vsc
  }
  
  vsc_mix = sig1
  if (Nd >= 1 & nsm >= 1) {
    imp_ps = sort(c(zeros_ps, small_cells))
    nd_imp_ps = nd[imp_ps]
    l_imp_ps = length(nd_imp_ps)
    
    sig1_imp = sig1[imp_ps]
    sig2_imp = sig2[imp_ps]
    
    vsc_mix_imp = 1:l_imp_ps*0
    for(k in 1:l_imp_ps){
      ndek = nd_imp_ps[k]
      s1k = sig1_imp[k]
      s2k = sig2_imp[k]
      if(ndek == 0){
        varek = s1k
      } else if (ndek == 1) {
        varek = s1k*.2 + s2k*.8
      } else if(ndek == 2 | ndek == 3){
        varek = s1k*.3 + s2k*.7
      }else if(ndek == 4 | ndek == 5){
        varek = s1k*.6 + s2k*.4
      }else if(ndek == 6 | ndek == 7){
        varek = s1k*.6 + s2k*.4
      }else if(ndek == 8 | ndek == 9){
        varek = s1k*1 + s2k*0
      }else if(ndek == 10){
        varek = s1k*.8 + s2k*.2
      }
      vsc_mix_imp[k] = varek 
    }
    
    vsc_mix[imp_ps] = vsc_mix_imp
  }
  hl = z.mult*sqrt(vsc_mix)
  lwr = muhat - hl
  upp = muhat + hl

  #new: test of H_0: theta in V and H_1: theta in C
  #do lsq fit with Atheta = 0 to get fit under H0
  #use the 1st option to compute T1 and T2
  if(test){
    #m = nrow(amat)
    m = qr(amat)$rank
    umat = chol(v1)
    uinv = ginv(umat)
    vmat = qr.Q(qr(t(amat)), complete = TRUE)[, -(1:(qr(t(amat))$rank)), drop = FALSE]
    #if(ncol(vmat>0)){
    #  pvmat = vmat %*% solve(crossprod(vmat), t(vmat))
    #}
    #imat = diag(M)
    
    #1st option
    vec = amat %*% yvecu
    aumat = amat %*% t(umat)
    auinv = ginv(aumat)
    vec2 = auinv %*% vec 
    T1_hat = t(vec2) %*% vec2
    #T1_hat = t(yvecu) %*% t(amat) %*% ginv(amat %*% v1 %*% t(amat)) %*% amat %*% yvecu
    #2nd option (works ok...)
    #if(ncol(vmat>0)){
    #  T1_hat =  t(yvecu) %*% t(imat - pvmat) %*% uinv %*% t(uinv) %*% (imat - pvmat) %*% yvecu
    #}else{
    #  T1_hat =  t(yvecu) %*% uinv %*% t(uinv) %*% yvecu
    #}
    
    eans = eigen(v1, symmetric=TRUE)
    evecs = eans$vectors
    evals = eans$values
    Z_s = uinv %*% yvecu
    L = evecs %*% diag(sqrt(evals))
    #Linv = solve(L)
    Linv = ginv(L)
    atil = amat %*% L
    Z_s = Linv %*% yvecu
    
    theta_hat = coneA(Z_s, atil, msg=FALSE)$thetahat
    T2_hat = t(Z_s-theta_hat) %*% (Z_s-theta_hat)
    
    #atil = amat %*% umat 
    #2nd option (works ok...)
    #T2_hat = t(yvecu - muhat) %*% uinv %*% t(uinv) %*% (yvecu - muhat) 
    bval = (T1_hat-T2_hat)/T1_hat
    sm = 1e-8
    
    if (bval > sm) {
      nloop = 100
      zerovec = rep(0, M)
      zsims = mvrnorm(nloop, mu = zerovec, Sigma = imat)
      
      J_length = rep(0, nloop)
      for (i in 1:nloop) {
        J_length[i] = length(coneA(zsims[i, ], atil, msg=FALSE)$face)
      }
      
      mdist = 0:m*0
      for (i in 1:(m+1)){
        mdist[i] = length(which(J_length == (i-1)))
      }
      mdist = mdist/nloop
      
      d0 = qr(vmat)$rank
      d = m-d0
      pval = 1 - sum(c(0, pbeta(bval, d:1/2, 1:d/2), 1)*mdist[1:(d+2)])
   } else {pval = 1}
  }
  ans = list(muhat = muhat[domain.ord], muhat.s = muhat, muhat.un = muhatu, 
             lwr = lwr[domain.ord], upp = upp[domain.ord], lwr.s = lwr, 
             upp.s = upp, lwru = lwru, uppu = uppu, domain = ID, 
             grid_ps = pskeep, amat = amat, Ds = Ds, domain.ord = domain.ord, 
             nd = nd, grid = grid, cov.un = vcov(ans.unc), cov.c = vsc_mix,
             ec = empty_cells, hkeep = hkeep, h_grid = h_grid, zeros_ps=zeros_ps, pval=pval, CIC=CIC)
  class(ans) = "csvy"
  invisible(ans)
  #return (ans)
}


####################################
#apply plotpersp on a csvy.fit object
#####################################
plotpersp <- function(object, x1 = NULL, x2 = NULL,...) {
  x1nm <- deparse(substitute(x1))
  x2nm <- deparse(substitute(x2))
  if (inherits(object, "cgamp")) {
    xnms_add <- object$object$xnms_add
  } else { #cgam and csvy
    xnms_add <- object$xnms_add
  }
  if (inherits(object, "wpsp")) {
    xnms_wp <- object$object$xnms_wp
  } else {
    xnms_wp <- object$xnms_wp
  }
  if (inherits(object, "trisplp")) {
    xnms_tri <- object$object$xnms_tri
  } else {
    xnms_tri <- object$xnms_tri
  }
  is_add <- all(c(any(grepl(x1nm, xnms_add, fixed = TRUE)), any(grepl(x2nm, xnms_add, fixed = TRUE))))
  is_wps <- all(c(any(grepl(x1nm, xnms_wp, fixed = TRUE)), any(grepl(x2nm, xnms_wp, fixed = TRUE))))
  is_tri <- all(c(any(grepl(x1nm, xnms_tri, fixed = TRUE)), any(grepl(x2nm, xnms_tri, fixed = TRUE))))
  if (missing(x1) | missing(x2)) {
    UseMethod("plotpersp")
  } else {
    cs = class(object)
    if (length(cs) == 1 & is.null(x1nm) & is.null(x2nm)) {
      UseMethod("plotpersp")
    } else {
      plotpersp.csvy(object, x1, x2, x1nm, x2nm,...)
    }
  }
}

####################################
plotpersp.csvy = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL, data=NULL, ci="none", transpose=FALSE, main=NULL, categ=NULL, col="white", cex.main=.8, xlab=NULL, 
                          ylab=NULL, zlab=NULL, zlim=NULL, box=TRUE, axes=TRUE, th=NULL, ltheta=NULL, ticktype="detailed", NCOL=NULL,...) {
  if (!inherits(object, "csvy")) {
    warning("calling plotpersp(<fake-csvy-object>) ...")
  }
  t_col = function(color, percent = 80, name = NULL) {
    rgb.val = col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                maxColorValue = 255,
                alpha = (100-percent)*255/100,names = name)
    ## Save the color
    invisible(t.col)
  }
  col.upp = t_col("green", perc = 90, name = "lt.green")
  col.lwr = t_col("pink", perc = 80, name = "lt.pink")
  Ds = object$Ds
  muhat = object$muhat.s
  lwr = object$lwr.s
  upp = object$upp.s
  xnms = object$xnms_add
  xmat = object$xmat_add
  bool = apply(xmat, 2, function(x) is.numeric(x))
  if(any(!bool)){
    xmat = apply(xmat, 2, function(x) as.numeric(x))
  }
  xmat0 = object$xmat_add0
  bool = apply(xmat0, 2, function(x) is.numeric(x))
  if(any(!bool)){
    xmat0 = apply(xmat0, 2, function(x) as.numeric(x))
  }
  shp = object$shapes_add
  ynm = object$ynm
  dm = object$domain
  #new: to avoid empty cell
  dm = as.numeric(levels(dm))[dm]
  nd = object$nd
  ec = object$ec
  data = object$data
  
  knms = length(xnms)
  obs = 1:knms
  if (is.null(x1nm) && is.null(x2nm)) {
    if (length(xnms) >= 2) {
      if (is.null(categ)) {
        x1nm = xnms[1]
        x2nm = xnms[2]
        x1id = 1
        x2id = 2
      } else {
        if (!is.character(categ)) {
          stop("categ must be a character argument!")
        } else if (any(grepl(categ, xnms))) {
          id = which(grepl(categ, xnms))
          x12nms = xnms[-id]
          x12ids = obs[-id]
          x1nm = x12nms[1]
          x2nm = x12nms[2]
          x1id = x12ids[1]
          x2id = x12ids[2]
        }  
      }
    } else {stop ("Number of non-parametric predictors must >= 2!")}
  } else {
    #new: use the data as the environment for x1 and x2
    x1 = data[[x1nm]]; x2 = data[[x2nm]]
    #print (length(x1))
    #print (length(x2))
    if (all(xnms != x1nm)) {
      if (length(x1) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x1!")
      }
      bool = apply(xmat, 2, function(x) all(x1 == x))
      if (any(bool)) {
        x1id = obs[bool]
        #change x1nm to be the one in formula
        x1nm = xnms[bool]
      } else {
        stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the model!"))
      }
    } else {
      x1id = obs[xnms == x1nm]
    }
    
    if (all(xnms != x2nm)) {
      if (length(x2) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x2!")
      }
      bool = apply(xmat, 2, function(x) all(x2 == x))
      if (any(bool)) {
        x2id = obs[bool]
        x2nm = xnms[bool]
      } else {
        stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the model!"))
      }
    } else {
      x2id = obs[xnms == x2nm]
    }
  }
  #xmat keeps the order of x's in the original data frame, but x1 and x2 may change depending on x1id and x2id
  x1 = xmat[, x1id]
  x2 = xmat[, x2id]
  
  if (is.null(xlab)) {
    #xlab = deparse(x1nm)
    xlab = x1nm
  }
  if (is.null(ylab)) {
    #ylab = deparse(x2nm)
    ylab = x2nm
  }
  if (is.null(zlab)) {
    zlab = paste("Est mean of", ynm)
  }
  
  D1 = Ds[x1id]
  D2 = Ds[x2id]
  
  grid2 = cbind(x1, x2, dm)
  ord = order(dm)
  #order x1 and x2 in ascending order as how we estimate muhat
  grid2 = grid2[ord, ,drop=FALSE]
  dm = sort(dm)
  
  ps = NULL
  ps_id = 1:nrow(grid2)
  knms = length(xnms)
  obs = 1:knms
  if (knms >= 3) {
    if (is.null(categ)) {
      x3id = obs[-c(x1id, x2id)]
      kx3 = length(x3id)
      for (i in 1:kx3) {
        x3i = xmat[, x3id[i]]
        x3i = x3i[ord]
        x3i0 = xmat0[, x3id[i]]
        x3i0 = x3i0[ord]
        x3i_use = max(x3i[x3i <= median(x3i)])
        ps = c(ps, x3i_use)
        ps_idi = which(x3i == x3i_use)
        ps_id = intersect(ps_id, ps_idi)
      }
      domain_id = unique(dm[ps_id])
      muhat_use = muhat[domain_id]
      upp_use = upp[domain_id]
      lwr_use = lwr[domain_id]
    } else {
      x3nms = xnms[-c(x1id, x2id)]
      x3mat = xmat[,-c(x1id, x2id),drop=FALSE]
      x3mat0 = xmat0[,-c(x1id, x2id),drop=FALSE]
      #print (head(x3mat))
      if (!is.character(categ)) {
        stop("categ must be a character argument!")
      } else if (any(grepl(categ, x3nms))) {
        id = which(grepl(categ, x3nms))
        x3nmi = x3nms[id]
        if (x3nmi %in% c(x1nm, x2nm)) {
          stop("categ must be different than x1 and x2!")
        }
        x3i = x3mat[,id]
        x3i = x3i[ord]
        x3i0 = x3mat0[,id]
        x3i0 = x3i0[ord]
        
        ps_id4 = NULL
        ps_id = 1:nrow(grid2)
        if (ncol(x3mat) > 1) {
          x4s = x3mat[,-id,drop=FALSE]
          x4s_use = apply(x4s, 2, min)
          for (i in 1:ncol(x4s)) {
            #print (i)
            #print (x4s)
            x4i = x4s[,i,drop=FALSE]
            x4i = x4i[ord]
            x4i_use = x4s_use[i]
            ps_id4 = which(x4i == x4i_use)
            ps_id = intersect(ps_id, ps_id4)
          }
        }
        
        surfs = list()
        dms = list()
        ux3i = unique(x3i)
        ux3i0 = unique(x3i0)
        kz_add = length(ux3i)
        mins = maxs = NULL
        
        dm = object$domain
        dm = sort(dm)
        for (iz in 1:kz_add) {
          x3i_use = ux3i[iz]
          ps_idi = which(x3i == x3i_use)
          ps_id_iz = intersect(ps_id, ps_idi)
          domain_id = unique(dm[ps_id_iz])
          dms[[iz]] = domain_id
          muhat_use = muhat[domain_id]
          surfs[[iz]] = muhat_use
          mins = c(mins, min(muhat_use))
          maxs = c(maxs, max(muhat_use))
        }
      } else {print(paste(categ, "is not an exact character name defined in the csurvey fit!"))}
    }
  } else {
    domain_id = unique(dm)
    muhat_use = muhat
    upp_use = upp
    lwr_use = lwr
  }
  
  if (is.null(categ)) {
    #if (!is.null(upp)&&!is.null(lwr)) {
    if (ci != 'none') {
      upp_use = upp
      lwr_use = lwr
      rg = max(upp_use, na.rm = TRUE) - min(lwr_use, na.rm = TRUE)
      z.lwr = min(muhat_use, na.rm = TRUE) - rg/5
      z.upp = max(muhat_use, na.rm = TRUE) + rg/5
    } else if (ci == 'none') {
      rg = max(muhat_use, na.rm = TRUE) - min(muhat_use, na.rm = TRUE)
      z.lwr = min(muhat_use, na.rm = TRUE) - rg/5
      z.upp = max(muhat_use, na.rm = TRUE) + rg/5
    }
  } else {
    z.lwr = min(mins, na.rm = TRUE)
    z.upp = max(maxs, na.rm = TRUE)
  }
  if (is.null(zlim)) {
    zlim0 = c(z.lwr, z.upp)
  } else {zlim0 = zlim}
  
  ang = NULL
  if (is.null(th)) {
    if (all(shp == 1)) {
      ang = -40
    } else if (all(shp == 2)){
      ang = 40
    } else if (shp[1] == 1 & shp[2] == 2) {
      ang = -50
    } else if (shp[1] == 2 & shp[2] == 1) {
      ang = -230
    } else {
      ang = -40
    }
  } else {
    ang = th
  }
  if (is.null(categ)){
    x1p = 1:D1
    x2p = 1:D2
    surf1 = surf2 = surf3 = matrix(0, nrow = D1, ncol = D2)
    #new: some x won't start from 1
    ux1 = sort(unique(x1))
    ux2 = sort(unique(x2))
    if (knms == 1 | knms == 2) {
      for(di in domain_id){
        rid = (grid2[grid2[,3] == di, ,drop=FALSE])[1, ]
        #new: some x won't start from 1
        ps1 = which(ux1 %in% rid[1])
        ps2 = which(ux2 %in% rid[2])
        surf1[ps1, ps2] = muhat_use[di]
        #surf1[rid[1], rid[2]] = muhat_use[di]
        if(!is.null(upp)&&!is.null(lwr)){
          surf2[ps1, ps2] = upp_use[di]
          surf3[ps1, ps2] = lwr_use[di]
          #surf2[rid[1], rid[2]] = upp_use[di]
          #surf3[rid[1], rid[2]] = lwr_use[di]
        }
      }
    } else if (knms >= 3) {
      for(di in domain_id){
        rid = (grid2[grid2[,3] == di, ,drop=FALSE])[1, ]
        #new: some x won't start from 1
        ps1 = which(ux1 %in% rid[1])
        ps2 = which(ux2 %in% rid[2])
        surf1[ps1, ps2] = muhat_use[which(domain_id %in% di)]
        #surf1[rid[1], rid[2]] = muhat_use[which(domain_id %in% di)]
        if(!is.null(upp)&&!is.null(lwr)){
          surf2[ps1, ps2] = upp_use[which(domain_id %in% di)]
          surf3[ps1, ps2] = lwr_use[which(domain_id %in% di)]
          #surf2[rid[1], rid[2]] = upp_use[which(domain_id %in% di)]
          #surf3[rid[1], rid[2]] = lwr_use[which(domain_id %in% di)]
        }
      }
    }
    
    #empty cell: test 3D
    if (length(ec) >= 1) {
      if (knms == 1 | knms == 2) {
        surf1[ec] = muhat_use[ec]
        if(!is.null(upp) && !is.null(lwr)){
          surf2[ec] = upp_use[ec]
          surf3[ec] = lwr_use[ec]
        }
      } else if (knms >= 3) {
        surf1[which(domain_id %in% ec)] = muhat_use[which(domain_id %in% ec)]
        if(!is.null(upp) && !is.null(lwr)){
          surf2[which(domain_id %in% ec)] = upp_use[which(domain_id %in% ec)]
          surf3[which(domain_id %in% ec)] = lwr_use[which(domain_id %in% ec)]
        }
      }
    }
    #suggested by cran: not to change users' setting
    oldpar <- par(no.readonly = TRUE)  
    on.exit(par(oldpar))    
    if (ci == "none") {
      persp(x1p, x2p, surf1, zlim = zlim0, col=col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes)
      par(new=FALSE)
    } else if (ci == 'up') {
      persp(x1p, x2p, surf1, zlim = zlim0, col=col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes)
      par(new = TRUE)
      persp(x1p, x2p, surf2, zlim = zlim0, col=col.upp, theta = ang, box=FALSE, axes=FALSE)
      par(new=FALSE)
    } else if (ci == 'lwr') {
      persp(x1p, x2p, surf1, zlim = zlim0, col=col, xlab = xlab, ylab = ylab, zlab = zlab, main = main, theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes)
      par(new = TRUE)
      persp(x1p, x2p, surf3, zlim = zlim0, col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
      par(new=FALSE)
    }
  } else {
    #new: some x won't start from 1
    ux1 = sort(unique(x1))
    ux2 = sort(unique(x2))
    palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
    ks = length(surfs)
    if (ks > 1 && ks < 11) {
      col = palette[1:ks]
    } else {
      #new: use rainbow
      col = topo.colors(ks)
    }
    width = ks^.5
    wd = round(width)
    if (width > wd) {
      wd = wd + 1
    }  
    if ((wd^2 - ks) >= wd ) {
      fm = c(wd, wd - 1)
    } else {
      fm = rep(wd, 2)
    }
    if (wd > 3) {
      cex = .7
      cex.main = .8
    }
    if (transpose) {
      fm = rev(fm)
    }
    #suggested by cran: not to change users' setting
    oldpar <- par(no.readonly = TRUE)  
    on.exit(par(oldpar))   
    if(is.null(NCOL)){
      par(mfrow = c(fm[1], fm[2]))
    } else {
      nr = fm[1]*fm[2]/NCOL
      par(mfrow = c(nr, NCOL))
    }
    
    par(mar = c(4, 1, 1, 1))
    par(cex.main = cex.main)
    
    #print (ks)
    x1p = 1:D1
    x2p = 1:D2
    for(i in 1:ks){
      surfi = surfs[[i]]
      dmi = dms[[i]]
      #print (surfi)
      surf1 = matrix(0, nrow = D1, ncol = D2)
      for(j in 1:(D1*D2)){
        di = dmi[j]
        rid = (grid2[grid2[,3] == di, ,drop=FALSE])[1, ]
        #new: some x won't start from 1
        ps1 = which(ux1 %in% rid[1])
        ps2 = which(ux2 %in% rid[2])
        surf1[ps1, ps2] = surfi[j]
        #surf1[rid[1], rid[2]] = surfi[j]
      }
      #new: avoid thick labs
      #box = TRUE
      #axes = TRUE
      #if (i > 1) {
      #  xlab = ylab = zlab = " "
      #  box = FALSE
      #  axes = FALSE
      #}
      persp(x1p, x2p, surf1, zlim = zlim0, col=col[i], xlab = xlab, ylab = ylab, zlab = zlab, main = paste0(x3nmi, " = ", ux3i0[i]), theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes)
      #par(new = TRUE)
    }
    par(new=FALSE)
  }
  #par(new=FALSE)
}


#########################################
#subroutines for confidence interval
#########################################
makebin = function(x) {
  k = length(x)
  r = 0
  for (i in 1:k) {
    r = r + x[k-i+1]*2^(i-1)
  }
  r
}

getvars = function(num) {
  i = num
  digits = 0
  power = 0
  while (digits == 0) {
    if (i < 2^power) {
      digits = power
    }
    power = power+1
  }
  binry = 1:digits*0
  if (num > 0) {
    binry[1] = 1
  }
  i = i - 2^(digits - 1)
  power = digits - 2
  for (p in power:0) {
    if (i >= 2^p) {
      i = i - 2^p
      binry[digits-p] = 1
    }
  }
  binry
}

getbin = function(num, capl) {
  br = getvars(num-1)
  digits = length(br)
  binrep = 1:capl*0
  binrep[(capl-digits+1):capl] = br
  binrep
}

##############################################
#eight shape functions
#add constr for the user-defined amat case
##############################################
constr <- function(x, numknots = 0, knots = 0, space = "E")
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 0
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

incr <- function(x, numknots = 0, knots = 0, space = "E")
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 1
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

#relaxed increasing
# relax.incr <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 1
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   attr(x, "relax") <- TRUE
#   #class(x) <- "additive"
#   x
# }


decr <- function(x, numknots = 0, knots = 0, space = "E")
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 2
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
  x
}

#relaxed decreasing
# relax.decr <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 2
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   attr(x, "relax") <- TRUE
#   #class(x) <- "additive"
#   x
# }

# conv <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 3
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

# conc <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 4
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

# incr.conv <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 5
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

# decr.conv <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 6
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

# incr.conc <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 7
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

# decr.conc <- function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl <- match.call()
#   pars <- match.call()[-1]
#   attr(x, "nm") <- deparse(pars$x)
#   attr(x, "shape") <- 8
#   attr(x, "numknots") <- numknots
#   attr(x, "knots") <- knots
#   attr(x, "space") <- space
#   attr(x, "categ") <- "additive"
#   #class(x) <- "additive"
#   x
# }

block.Ord <- function(x, order = NULL, numknots = 0, knots = 0, space = "E")
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 9
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "block.ord"
  attr(x, "order") <- order
  #class(x) <- "additive"
  x
}

block.Ave <- function(x, order = NULL, numknots = 0, knots = 0, space = "E")
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- 10
  attr(x, "numknots") <- numknots
  attr(x, "knots") <- knots
  attr(x, "space") <- space
  attr(x, "categ") <- "block.ave"
  attr(x, "order") <- order
  #class(x) <- "additive"
  x
}

fitted.csvy <- function(object,...) {
  ans <- object$muhat
  return(ans)
}

######
#1D
######
makeamat = function(x, sh, Ds = NULL, suppre = FALSE, interp = FALSE, 
                    relaxes = FALSE, h = NULL, block.ave = NULL, 
                    block.ord = NULL, zeros_ps = NULL, noord_cells = NULL) {
  n = length(x)
  xu = sort(unique(x))
  n1 = length(xu)
  sm = 1e-7
  ms = NULL
  obs = 1:n
  Nd = length(Ds)
  M = rev(x)[1]
  if (relaxes) {
    #create a list of amat's first
    d = l = x
    amat = matrix(0, nrow=(M-1), ncol=M)
    for(k in 1:(M-1)){
      amat[k, ] = exp(-abs(l-(k+1))/h)/sum(exp(-abs(d-(k+1))/h))-exp(-abs(l-k)/h)/sum(exp(-abs(d-k)/h))
    }
    # return (amat)
  } else {
    #block.ave and block.ord = -1 if sh is between 1 and 9
    if (all(block.ave == -1) & all(block.ord == -1)) {
      #print ('a')
      if (sh < 3) {
        #if sh=0, then keep the zero matrix
        amat = matrix(0, nrow = n1 - 1, ncol = n)
        if (sh > 0) {
          for (i in 1:(n1-1)) {
            c1 = min(obs[abs(x - xu[i]) < sm]); c2 <- min(obs[abs(x - xu[i + 1]) < sm])
            amat[i, c1] = -1; amat[i, c2] = 1
          }
          if (sh == 2) {amat = -amat}
        }
        #return (amat)
      } else if (sh == 3 | sh == 4) {
        #  convex or concave
        amat = matrix(0, nrow = n1 - 2 ,ncol = n)
        for (i in 1: (n1 - 2)) {
          c1 = min(obs[x == xu[i]]); c2 = min(obs[x == xu[i+1]]); c3 = min(obs[x == xu[i+2]])
          amat[i, c1] = xu[i+2] - xu[i+1]; amat[i, c2] = xu[i] - xu[i+2]; amat[i, c3] = xu[i+1] - xu[i]
        }
        if (sh == 4) {amat = -amat}
        #return (amat)
      } else if (sh > 4 & sh < 9) {
        amat = matrix(0, nrow = n1 - 1, ncol = n)
        for (i in 1:(n1 - 2)) {
          c1 = min(obs[x == xu[i]]); c2 = min(obs[x == xu[i + 1]]); c3 = min(obs[x == xu[i + 2]])
          amat[i, c1] = xu[i + 2] - xu[i + 1]; amat[i, c2] = xu[i] - xu[i + 2]; amat[i, c3] = xu[i + 1] - xu[i]
        }
        if (sh == 5) { ### increasing convex
          c1 = min(obs[x == xu[1]]); c2 = min(obs[x == xu[2]])
          amat[n1 - 1, c1] = -1; amat[n1 - 1, c2] = 1
          #return (amat)
        } else if (sh == 6) {  ## decreasing convex
          c1 = min(obs[x == xu[n1]]); c2 = min(obs[x == xu[n1 - 1]])
          amat[n1 - 1, c1] = -1; amat[n1 - 1, c2] = 1
          #return (amat)
        } else if (sh == 7) { ## increasing concave
          amat = -amat
          c1 = min(obs[x == xu[n1]]); c2 = min(obs[x == xu[n1 - 1]])
          amat[n1 - 1, c1] = 1; amat[n1 - 1, c2] = -1
          #return (amat)
        } else if (sh == 8) {## decreasing concave
          amat = -amat
          c1 = min(obs[x == xu[1]]); c2 = min(obs[x == xu[2]])
          amat[n1 - 1, c1] = 1; amat[n1 - 1, c2] = -1
          #return (amat)
        }
      } 
      return (amat)
    } else if (all(block.ave == -1) & (length(block.ord) > 1)) {
      #nbcks is the total number of blocks
      block.ord_nz = block.ord
      
      noord_ps = which(block.ord == 0)
      rm_id = unique(sort(c(zeros_ps, noord_ps)))
      #rm_id = unique(sort(c(zeros_ps, noord_cells)))
      if(length(rm_id)>0){
        block.ord_nz = block.ord[-rm_id]
      }
      
      ubck = unique(block.ord_nz)
      #use values in block.ord as an ordered integer vector
      ubck = sort(ubck)
      
      nbcks = length(table(block.ord_nz))
      szs = as.vector(table(block.ord_nz))
      #bck_lst = f_ecl(x, nbcks, szs)
      ubck = unique(block.ord_nz)
      #use values in block.ord as an ordered integer vector
      ubck = sort(ubck)
      
      #amat dimension: l1*l2...*lk by M
      amat = NULL
      for(i in 1:(nbcks-1)) {
        #bck_1 = bck_lst[[i]]; bck_2 = bck_lst[[i+1]]
        #l1 = length(bck_1); l2 = length(bck_2)
        ubck1 = ubck[i]
        ubck2 = ubck[i+1]
        bck_1 = which(block.ord_nz == ubck1)
        bck_2 = which(block.ord_nz == ubck2)
        l1 = length(bck_1)
        l2 = length(bck_2)
        
        #nr = l1*l2; nc = M
        #for each element in block1, the value for each element in block 2 will be the same
        #for each element in block1, the number of comparisons is l2
        amat_bcki = NULL
        #rm_l = length(zeros_ps) + length(noord_cells)
        rm_l = length(zeros_ps) + length(noord_ps)
        M.ord = length(block.ord)
        amatj = matrix(0, nrow = l2, ncol = M.ord-rm_l)
        bck_1k = bck_1[1]
        amatj[, bck_1k] = -1
        row_pointer = 1
        for(j in 1:l2){
          bck_2j = bck_2[j]
          amatj[row_pointer, bck_2j] = 1
          row_pointer = row_pointer + 1
        }
        amatj0 = amatj
        amat_bcki = rbind(amat_bcki, amatj)
        
        for(k in 2:l1) {
          bck_1k = bck_1[k]; bck_1k0 = bck_1[k-1]
          amatj = amatj0
          #set the value for the kth element in block 1 to be -1
          #set the value for the (k-1)th element in block 1 back to 0
          #keep the values for block 2
          amatj[, bck_1k] = -1; amatj[, bck_1k0] = 0
          amatj0 = amatj
          amat_bcki = rbind(amat_bcki, amatj)
        }
        amat = rbind(amat, amat_bcki)
      }
      
      if(length(noord_ps)>0){
        nr = nrow(amat); nc = ncol(amat)
        amat_tmp = matrix(0, nrow = nr, ncol = (M.ord-rm_l+length(noord_ps)))
        if(length(zeros_ps) > 0){
          ps = which(block.ord[-zeros_ps] != 0)
        }else{
          ps = which(block.ord != 0)
        }
        amat_tmp[, ps] = amat
        amat = amat_tmp
      }
      #print (amat[,1])
      return (amat)
    } else if ((length(block.ave) > 1) & all(block.ord == -1)) {
      block.ave_nz = block.ave
      #rm_id = unique(sort(c(zeros_ps, noord_cells)))
      
      noord_ps = which(block.ave == 0)
      rm_id = unique(sort(c(zeros_ps, noord_ps)))
      
      if(length(rm_id)>0){
        block.ave_nz = block.ave[-rm_id]
      }

      ubck = unique(block.ave_nz)
      #use values in block.ord as an ordered integer vector
      ubck = sort(ubck)
      nbcks = length(table(block.ave_nz))
      
      szs = as.numeric(table(block.ave_nz))
      
      rm_l = length(zeros_ps) + length(noord_ps)
      M.ave = length(block.ave)
      amat = matrix(0, nrow = (nbcks-1), ncol = (M.ave-rm_l))
      for(i in 1:(nbcks-1)) {
        ubck1 = ubck[i]
        ubck2 = ubck[i+1]
        ps1 = which(block.ave_nz == ubck1)
        ps2 = which(block.ave_nz == ubck2)
        sz1 = length(ps1)
        sz2 = length(ps2)
        amat[i, ps1] = -1/sz1
        amat[i, ps2] = 1/sz2
      }
      
      if(length(noord_ps)>0){
        amat_tmp = matrix(0, nrow = (nbcks-1), ncol = (M.ave-rm_l+length(noord_ps)))
        if(length(zeros_ps) > 0){
          ps = which(block.ave[-zeros_ps] != 0)
        }else{
          ps = which(block.ave != 0)
        }
        amat_tmp[, ps] = amat
        amat = amat_tmp
      }
      return (amat)
    } else if ((length(block.ave) > 1) & (length(block.ord) > 1)) {
      stop ("only one of block average and block ordering can be used!")
    }
  }
  #return (amat)
}

######
#>=2D
######
makeamat_2D = function(x=NULL, sh, grid2=NULL, zeros_ps=NULL, noord_cells=NULL,
                       Ds = NULL, suppre = FALSE, interp = FALSE, 
                       relaxes = FALSE, block.ave.lst = NULL, block.ord.lst = NULL) {
  sm = 1e-7
  ms = NULL
  Nd = length(Ds)
  ne = length(zeros_ps)
  if (Nd == 2) {
    if (any(sh > 0)) {
      D1 = Ds[1]
      D2 = Ds[2]
      obs = 1:(D1*D2)
      if (sh[1] < 9) {
        if (length(zeros_ps) > 0) {
          amat0_lst = list()
          grid3 = grid2[-zeros_ps, ,drop=FALSE]
          row_id = as.numeric(rownames(grid3))
          cts = row_id[which(row_id%%D1 == 0)]
          #temp: check if zeros_ps contains some point %% D1 = 0
          if (any((zeros_ps %% D1) == 0)) {
            zeros_ps_at_cts = zeros_ps[which((zeros_ps %% D1) == 0)]
            cts_add = row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
            cts = sort(c(cts, cts_add))
          }
          obs = 1:nrow(grid3)
          for (i in 1:length(cts)) {
            if (i == 1) {
              st = 1
              #ed = cts[1]
              ed = obs[row_id == cts[1]]
            } else {
              st = obs[row_id == cts[i-1]]+1
              ed = obs[row_id == cts[i]]
            }
            gridi = grid3[st:ed,,drop=FALSE]
            na_ps = which(apply(gridi, 1, function(x) all(is.na(x))))
            if (length(na_ps) > 0) {
              gridi = gridi[-na_ps,,drop=FALSE]
            }
            obsi = as.numeric(rownames(gridi))
            vec = (D1*(i-1) + 1):(D1*i)
            #vec = 1:D1
            zerosi = vec[-which(vec %in% obsi)]
            non_zerosi = vec[which(vec %in% obsi)]
            D1i = nrow(gridi)
            #if ed == st, it could be: 1,2,3,4 with 0,0,0,92 observations
            #in this case, amat0i = NULL, and it will be removed before calling bdiag
            if (length(zerosi) >= 1 & (ed > st)) {
              amat0i = makeamat(1:D1i, sh[1])
              amat0_lst[[i]] = amat0i
            } else if (length(zerosi) >= 1 & (ed == st)) {
              #next
              amat0i = 0
              amat0_lst[[i]] = amat0i
            } else if (length(zerosi) == 0) {
              amat0i = makeamat(1:D1, sh[1])
              amat0_lst[[i]] = amat0i
            }
            #amat0_lst[[i]] = amat0i
          }
        } else {
          amat0 = makeamat(1:D1, sh=sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]], 
                           zeros_ps = zeros_ps, noord_cells = noord_cells[[1]])
          amat0_lst = rep(list(amat0), D2)
        }
        amat0_lst = Filter(Negate(is.null), amat0_lst)
        amat1 = bdiag(amat0_lst)
        amat1 = as.matrix(amat1)
      }
      #print (dim(amat1))
      #new
      if (sh[1] == 9 | sh[1] == 10) {
        #revise order vector when D >= 2, temp
        block.ave1 = block.ave.lst[[1]]
        block.ord1 = block.ord.lst[[1]]
        if (length(block.ave1) > 1) {
          block.ave1 = rep(block.ave1, D2)
        }
        if (length(block.ord1) > 1) {
          block.ord1 = rep(block.ord1, D2)
        }
        amat1 = makeamat(1:D1, sh[1], block.ave = block.ave1, block.ord = block.ord1, 
                         zeros_ps = zeros_ps, noord_cells = noord_cells[[1]])
        #amat1 = makeamat(1:D1, sh[1], block.ave = block.ave.lst[[1]], 
        #                 block.ord = block.ord.lst[[1]], zeros_ps = zeros_ps, noord_cells = noord_cells[[1]])
      }
      amat2 = NULL
      if (sh[2] < 3) {
        nr2 = D1*(D2-1)
        nc2 = D1*D2
        if (length(zeros_ps) == 0) {
          amat2 = matrix(0, nrow=nr2, ncol=nc2)
          #new: if sh[2] == 0, keep the zero matrix
          if (sh[2]>0){
            for(i in 1:nr2){
              amat2[i,i] = -1; amat2[i,(i+D1)] = 1
            }
          }
        } else {
          #temp
          #amat2 = matrix(0, nrow=1, ncol=nc2)
          amat2_lst = vector('list', length = nr2)
          rm_id = NULL
          for(i in 1:nr2) {
            if (i %in% zeros_ps) {
              if (i <= D1) {
                amat2_lst[[i]] = 1:nc2*0
                rm_id = c(rm_id, i)
              } else if (i <= (nc2-D1) & i > D1) {
                #ch_ps = c(ch_ps, i-D1)
                amat2_lst[[i]] = 1:nc2*0
                rm_id = c(rm_id, i)
                rowi = 1:nc2*0
                rowi[i-D1] = -1; rowi[i+D1] = 1
                amat2_lst[[i-D1]] = rowi
              }
            } else {
              rowi = 1:nc2*0
              rowi[i] = -1; rowi[i+D1] = 1
              amat2_lst[[i]] = rowi
              #amat2 = rbind(amat2, amat2i)
            }
          }
          for(i in (nr2+1):nc2) {
            if (i %in% zeros_ps)  {
              amat2_lst[[i-D1]] = 1:nc2*0
              rm_id = c(rm_id, i-D1)
            }
          }
          if (!is.null(rm_id)){
            amat2_lst = amat2_lst[-rm_id]
          }
          amat2 = NULL
          for(i in 1:length(amat2_lst)){
            amat2 = rbind(amat2, amat2_lst[[i]])
          }
          amat2 = amat2[,-zeros_ps,drop=FALSE]
        }
        #new: if sh[2] == 0, keep the zero matrix
        if (sh[2] == 0) {
          #nr2 = nrow(amat2)
          #nc2 = ncol(amat2)
          #amat2 = matrix(0, nrow=nr2, ncol=nc2)
          amat2 = NULL
        }
        if (sh[2] == 2) {amat2 = -amat2}
      } else if (sh[2] == 3 | sh[2] == 4) {
        nr2 = D1*(D2-2)
        nc2 = D1*D2
        amat2 = matrix(0, nrow=nr2, ncol=nc2)
        xu = 1:D2
        for (i in 1:nr2) {
          #c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
          amat2[i, i] = 1; amat2[i, (i+D1)] = -2; amat2[i, (i+2*D1)] = 1
        }
        if (sh[2] == 4) {amat2 = -amat2}
      } else if (sh[2] > 4 & sh[2] < 11) {
        nr2 = D1*(D2-2)
        nc2 = D1*D2
        amat2 = matrix(0, nrow=nr2, ncol=nc2)
        xu = 1:D2
        for (i in 1:nr2) {
          #c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
          #amat2[i, c1] = 1; amat2[i, (c2+D1)] = -2; amat2[i, (c3+2*D1)] = 1
          amat2[i, i] = 1; amat2[i, (i+D1)] = -2; amat2[i, (i+2*D1)] = 1
        }
        if (sh[2] == 5) { ### increasing convex
          amat2_add = matrix(0, nrow=D1, ncol=nc2)
          #c1 = min(obs[xu == xu[1]])
          #c2 = min(obs[xu == xu[2]])
          #amat2[nr2, c1] = -1; amat2[nr2, (c1+D1)] = 1
          for(i in 1:D1) {
            amat2_add[i,i] = -1; amat2_add[i,i+D1] = 1
          }
          amat2 = rbind(amat2, amat2_add)
        } else if (sh[2] == 6) {  ## decreasing convex
          #c1 = min(obs[xu == xu[D2]]); c2 = min(obs[xu == xu[D2 - 1]])
          amat2_add = matrix(0, nrow=D1, ncol=nc2)
          for(i in 1:D1) {
            amat2_add[i, i+2*D1] = -1; amat2_add[i, i+D1] = 1
          }
          amat2 = rbind(amat2, amat2_add)
        } else if (sh[2] == 7) { ## increasing concave
          amat2 = -amat2
          #c1 = min(obs[xu == xu[D2]]); c2 = min(obs[xu == xu[D2 - 1]])
          #amat2[nr2, c1+D1] = 1; amat2[nr2, c2] = -1
          amat2_add = matrix(0, nrow=D1, ncol=nc2)
          for(i in 1:D1) {
            amat2_add[i, i+2*D1] = 1; amat2_add[i, i+D1] = -1
          }
          amat2 = rbind(amat2, amat2_add)
        } else if (sh[2] == 8) {## decreasing concave
          amat2 = -amat2
          #c1 = min(obs[xu == xu[1]]); c2 = min(obs[xu == xu[2]])
          #amat2[nr2, c1] = 1; amat2[nr2, (c2+D1)] = -1
          amat2_add = matrix(0, nrow=D1, ncol=nc2)
          for(i in 1:D1) {
            amat2_add[i,i] = 1; amat2_add[i,i+D1] = -1
          }
          amat2 = rbind(amat2, amat2_add)
        } else if (sh[2] == 9 | sh[2] == 10) {## block ordering or block average
          #print (block.ave.lst[[2]])
          #print (block.ord.lst[[2]])
          #print (D2)
          #temp: revise the order vector
          block.ave2 = block.ave.lst[[2]]
          block.ord2 = block.ord.lst[[2]]
          block.ave2 = rep(block.ave2, each=D1)
          block.ord2 = rep(block.ord2, each=D1)
          amat2 = makeamat(1:D2, sh[2], block.ave = block.ave2, 
                           block.ord = block.ord2, 
                           zeros_ps = zeros_ps, 
                           noord_cells = noord_cells[[2]])
        }
      }
      amat = rbind(amat1, amat2)
    } 
  } else if (Nd >= 3) {
    M = length(Ds)
    D1 = Ds[1]
    cump = cumprod(Ds)
    nc = cump[M]
    amat1 = amat2 = amatm = NULL
    
    # 1D
    if(sh[1] == 9 | sh[2] == 10){
      #temp: revise the order vector
      block.ave1 = block.ave.lst[[1]]
      block.ord1 = block.ord.lst[[1]]
      rpts = rev(cumprod(Ds[-1]))[1]
      if (length(block.ave1) > 1) {
        block.ave1 = rep(block.ave1, rpts)
      }
      if (length(block.ord1) > 1) {
        block.ord1 = rep(block.ord1, rpts)
      }
      block.ave.lst[[1]] = block.ave1
      block.ord.lst[[1]] = block.ord1
      amat1 = makeamat(1:D1, sh[1], block.ave = block.ave1, block.ord = block.ord1, 
                       zeros_ps = zeros_ps, noord_cells = noord_cells[[1]])
    } else {
      if(length(zeros_ps)>0){
        amat0_lst = list()
        grid3 = grid2[-zeros_ps, ,drop=FALSE]
        row_id = as.numeric(rownames(grid3))
        cts = row_id[which(row_id%%D1 == 0)]
        if(any((zeros_ps %% D1) == 0)){
          zeros_ps_at_cts = zeros_ps[which((zeros_ps %% D1) == 0)]
          cts_add = row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
          cts = sort(c(cts, cts_add))
        }
        obs = 1:nrow(grid3)
        
        for (i in 1:length(cts)) {
          if (i == 1) {
            st = 1
            ed = obs[row_id == cts[1]]
          } else {
            st = obs[row_id == cts[i-1]]+1
            ed = obs[row_id == cts[i]]
          }
          gridi = grid3[st:ed,,drop=FALSE]
          na_ps = which(apply(gridi, 1, function(x) all(is.na(x))))
          if (length(na_ps) > 0) {
            gridi = gridi[-na_ps, ,drop=FALSE]
          }
          obsi = as.numeric(rownames(gridi))
          vec = (D1*(i-1) + 1):(D1*i)
          zerosi = vec[-which(vec %in% obsi)]
          non_zerosi = vec[which(vec %in% obsi)]
          D1i = nrow(gridi)
          if (length(zerosi) >= 1 & (ed > st)) {
            amat0i = makeamat(1:D1i, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            amat0_lst[[i]] = amat0i
          } else if (length(zerosi) >= 1 & (ed == st)) {
            amat0i = 0
            amat0_lst[[i]] = amat0i
          } else if (length(zerosi) == 0) {
            amat0i = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            amat0_lst[[i]] = amat0i
          }
        }
        amat0_lst = Filter(Negate(is.null), amat0_lst)
        amat1 = bdiag(amat0_lst)
        amat1 = as.matrix(amat1)
      }else{
        n1 = nc/D1
        amat1_0 = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
        amat1_0_lst = rep(list(amat1_0), n1)
        amat1 = bdiag(amat1_0_lst)
        amat1 = as.matrix(amat1)
      }
    }
    amat = amat1
    
    #2D
    for (i in 2:(M-1)) {
      Di = Ds[i]
      cump = cumprod(Ds[1:i])
      nci = cump[i]
      gap = cump[i-1]
      ni = nc/nci
      
      #temp: revise the order vector
      block.avei = block.ave.lst[[i]]
      block.ordi = block.ord.lst[[i]]
      if (length(block.avei) > 1) {
        tmp.avei = rep(block.avei, each = gap)
        block.avei = rep(tmp.avei, ni)
      }
      if (length(block.ordi) > 1) {
        tmp.ordi = rep(block.ordi, each = gap)
        block.ordi = rep(tmp.ordi, ni)
      }
      block.ave.lst[[i]] = block.avei
      block.ord.lst[[i]] = block.ordi
      
      if (sh[i] %in% c(9,10)) {
        amati = makeamat(1:Di, sh[i], block.ave = block.ave.lst[[i]], 
                         block.ord = block.ord.lst[[i]], zeros_ps = zeros_ps, noord_cells = noord_cells[[i]])
      } else {
        if ((sh[i] %in% c(1,2))) {
          nr = gap*(Di - 1)
        } else {
          if (Di < 3) {
            stop ("For monotonicity + convexity constraints, number of domains must be >= 3!")
          } else {
            nr = gap*(Di - 2)
          }
        }
        amati_0 = makeamat_2d(sh=sh[i], nr2=nr, nc2=nci, gap=gap, D2=Di, Ds=Ds)
        
        if(length(zeros_ps) == 0) {
          amati_0_lst = rep(list(amati_0), ni)
        } else {
          amati_0_lst = vector("list", length = ni)
          cts = seq(nci, nc, length=ni)
          bool_lst = map(zeros_ps, .f = function(.x) as.numeric(cts >= .x))
          #ps_lst keeps track of the cut points such that the point giving -1 is the upp bd 
          #of the interval containining one zero cell
          ps = map_dbl(bool_lst, .f = function(.x) min(which(.x == 1)))
          cts_check = cts[ps]
          for(k in 1:ni){
            if(!(cts[k] %in% cts_check)){
              #no zero cell in this interval
              amati_0_lst[[k]] = amati_0
            } else {
              #find out the empty cells within an interval
              if(k == 1) {
                zeros_ctsk = zeros_ps[zeros_ps <= cts[k]]
              }else{
                zeros_ctsk = zeros_ps[zeros_ps <= cts[k] & zeros_ps > cts[k-1]]
              }
              nek = length(zeros_ctsk)
              rm_id = NULL
              amati_0_cp = amati_0; nci = ncol(amati_0)
              col_ps_vec = NULL
              for(k2 in 1:nek){
                zero_k2 = zeros_ctsk[k2]
                #check: zero_k2 %% nr = 0
                if(zero_k2 != cts[k]){
                  col_ps = zero_k2 %% nci
                }else{
                  col_ps = nci
                }
                col_ps_vec = c(col_ps_vec, col_ps)
                
                row_ps = col_ps - gap
                #if((zero_k2 %% nr + gap) > nr){
                if((col_ps + gap) > nci){
                  rm_id = c(rm_id, row_ps)
                }else{
                  rm_id = c(rm_id, col_ps)
                  amati_0_cp[row_ps, (col_ps+gap)] = 1
                  amati_0_cp[row_ps, col_ps] = 0
                }
              }
              
              amati_0_cp = amati_0_cp[-rm_id, ]
              #amati_0_cp = amati_0_cp[,-(zeros_ctsk%%nr)]
              amati_0_cp = amati_0_cp[,-col_ps_vec]
              amati_0_lst[[k]] = amati_0_cp
            }
          }
        }
        
        amati = bdiag(amati_0_lst)
        amati = as.matrix(amati)
        #new: sh[i] == 0
        if(sh[i] == 0){
          #nri = nrow(amati)
          #nci = ncol(amati)
          #amati = matrix(0, nrow=nri, ncol=nci)
          amati = NULL
        }
      }
      amat = rbind(amat, amati)
    }
    
    #3D
    cump = cumprod(Ds[1:(M-1)])
    gap = cump[M-1]
    Dm = Ds[M]
    if((sh[M] %in% c(9,10))){
      #temp: revise the order vector
      block.avei = block.ave.lst[[M]]
      block.ordi = block.ord.lst[[M]]
      #ni=1
      if (length(block.avei) > 1) {
        block.avei = rep(block.avei, each = gap)
        #block.avei = rep(tmp.avei, ni)
      }
      if (length(block.ordi) > 1) {
        block.ordi = rep(block.ordi, each = gap)
        #block.ordi = rep(tmp.ordi, ni)
      }
      block.ave.lst[[M]] = block.avei
      block.ord.lst[[M]] = block.ordi
      
      amatm = makeamat(1:Dm, sh[M], block.ave = block.ave.lst[[M]],
                       block.ord = block.ord.lst[[M]], zeros_ps = zeros_ps, 
                       noord_cells = noord_cells[[M]])
    } else {
      if ((sh[M] %in% c(1,2))) {
        nr = gap*(Dm-1)
      } 
      if (sh[M] >= 3) {
        nr = gap*(Dm-2)
      }
      #new: sh[M] == 0
      if(sh[M] == 0){
        #nrm = nrow(amatm)
        #ncm = ncol(amatm)
        #amatm = matrix(0, nrow=nrm, ncol=ncm)
        amatm = NULL
      } else {
        amatm = makeamat_2d(sh=sh[M], nr2=nr, nc2=nc, gap=gap, D2=Dm, Ds=Ds)
        ncm = ncol(amatm)
        if(length(zeros_ps) > 0){
          amatm_cp = amatm
          rm_id = NULL
          for(k in 1:ne){
            zero_k = zeros_ps[k]
            col_ps = zero_k 
            row_ps = col_ps - gap
            if((col_ps + gap) > ncm){
              rm_id = c(rm_id, row_ps)
            }else{
              rm_id = c(rm_id, col_ps)
              amatm_cp[row_ps, (col_ps+gap)] = 1
              amatm_cp[row_ps, col_ps] = 0
            }
          }
          amatm_cp = amatm_cp[-rm_id, ]
          amatm_cp = amatm_cp[,-zeros_ps]
          amatm = amatm_cp
        }
      }
    }
    amat = rbind(amat, amatm)
  }
  return (amat)
}

####
#local function to be called in makeamat_2D
####
makeamat_2d = function(x=NULL, sh, nr2, nc2, gap, D1=NULL, D2, Ds = NULL, suppre = FALSE, interp = FALSE) {
  sm = 1e-7
  ms = NULL
  Nd = length(Ds)
  
  amat2 = NULL
  if (sh < 3) {
    amat2 = matrix(0, nrow=nr2, ncol=nc2)
    for(i in 1:nr2){
      amat2[i,i] = -1; amat2[i,(i+gap)] = 1
    }
    if (sh == 2) {amat2 = -amat2}
    return (amat2)
  } else if (sh == 3 | sh == 4) {
    amat2 = matrix(0, nrow=nr2, ncol=nc2)
    #xu = 1:D2
    for (i in 1:nr2) {
      #c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
      #amat2[i, c1] = 1; amat2[i, (c2+gap)] = -2; amat2[i, (c3+2*gap)] = 1
      amat2[i, i] = 1; amat2[i, (i+gap)] = -2; amat2[i, (i+2*gap)] = 1
    }
    if (sh == 4) {amat2 = -amat2}
    return (amat2)
  } else if (sh > 4 & sh < 9) {
    amat2 = matrix(0, nrow=nr2, ncol=nc2)
    #xu = 1:D2
    for (i in 1:nr2) {
      #c1 = min(obs[xu == xu[i]]); c2 = min(obs[xu == xu[i+1]]); c3 = min(obs[xu == xu[i+2]])
      #amat2[i, c1] = 1; amat2[i, (c2+gap)] = -2; amat2[i, (c3+2*gap)] = 1
      amat2[i, i] = 1; amat2[i, (i+gap)] = -2; amat2[i, (i+2*gap)] = 1
    }
    if (sh == 5) { ### increasing convex
      nr2_add = gap
      amat2_add = matrix(0, nrow=nr2_add, ncol=nc2)
      for(i in 1:nr2_add) {
        amat2_add[i,i] = -1; amat2_add[i,i+gap] = 1
      }
      amat2 = rbind(amat2, amat2_add)
      return (amat2)
    } else if (sh == 6) {  ## decreasing convex
      nr2_add = gap
      amat2_add = matrix(0, nrow=nr2_add, ncol=nc2)
      for(i in 1:nr2_add) {
        amat2_add[i, i+2*gap] = -1; amat2_add[i, i+gap] = 1
      }
      amat2 = rbind(amat2, amat2_add)
      return (amat2)
    } else if (sh == 7) { ## increasing concave
      nr2_add = gap
      amat2 = -amat2
      amat2_add = matrix(0, nrow=nr2_add, ncol=nc2)
      for(i in 1:nr2_add) {
        amat2_add[i, i+2*gap] = 1; amat2_add[i, i+gap] = -1
      }
      amat2 = rbind(amat2, amat2_add)
      return (amat2)
    } else if (sh == 8) {## decreasing concave
      nr2_add = gap
      amat2 = -amat2
      amat2_add = matrix(0, nrow=nr2_add, ncol=nc2)
      for(i in 1:nr2_add) {
        amat2_add[i,i] = 1; amat2_add[i,i+gap] = -1
      }
      amat2 = rbind(amat2, amat2_add)
      return (amat2)
    }
  }
  #return (amat2)
  #rslt = list(amat = amat,  ms = ms)
  #rslt
}


###############
#print method
###############
print.csvy = function(x, ...) 
{
  print(x$family)
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$design)
  invisible(x)
}

#########################
#empty cell imputation 
##########################
impute_em2 = function(empty_cells, obs_cells, M, yvec, Nd, nd, Nhat, w, domain.ord, grid2, Ds=NULL, sh=NULL, muhat=NULL, lwr=NULL, upp=NULL, vsc=NULL)
{
  ne = length(empty_cells)
  #observed sample size in cells with 0 or 1 n
  nd_ne = nd[empty_cells]
  zeros = ones = 0
  zeros_ps = NULL
  if (ne >= 1) {
    #ignore nd=1 now
    if (any(nd == 1)) {
       ones_ps = which(nd == 1)
       ones = length(ones_ps)
    #   yvec_ones = yvec[which(obs_cells %in% ones_ps)]
    #   Nhat_ones = Nhat[which(obs_cells %in% ones_ps)]
    } 
    if (any(nd == 0)) {
      zeros_ps = which(nd == 0)
      zeros = length(zeros_ps)
    }
    zeros = ne - ones
  }
  
  lc = rc = 1:ne*0
  if (Nd == 1) {
    vsc_all = lwr_all = upp_all = muhat_all = domain.ord_all = 1:M*0
    nslst = pslst = vector('list', ne)
    vsce = lwre = uppe = muhate = domain.ord_em = NULL
    for(k in 1:ne){
      e = empty_cells[k]
      ndek = nd_ne[k]
      #one_k = 1
      #new: consider n = 1 
      #at_max = all(obs_cells < e | (ndek == 1 & e == max(obs_cells))
      #at_min = all(obs_cells > e | (ndek == 1 & e == min(obs_cells))
      at_max = all(obs_cells < e | e == max(obs_cells))
      at_min = all(obs_cells > e | e == min(obs_cells))
      if (at_max) {
        lc = max(obs_cells[obs_cells < e])
        rc = lc       
        dmem = obs_cells[which(obs_cells%in%lc)]
      }
      if (at_min) {
        rc = min(obs_cells[obs_cells > e])
        lc = rc
        dmem = obs_cells[which(obs_cells%in%rc)]
      }
      if (!any(at_max) & !any(at_min)) {
        lc = max(obs_cells[obs_cells < e])
        rc = min(obs_cells[obs_cells > e])
        dmem = obs_cells[which(obs_cells%in%lc)]
      }
      n_lc = nd[lc]
      n_rc = nd[rc]
      
      # yl = yvec[which(obs_cells%in%lc)]
      # yr = yvec[which(obs_cells%in%rc)]
      # 
      # Nl = Nhat[lc]
      # Nr = Nhat[rc]
      #new: switch lc and rc if sh=2
      if (sh == 2) {
        lc_0 = rc; rc_0 = lc
        lc = lc_0; rc = rc_0
      }
      if (!is.null(muhat)) {
        mul = muhat[which(obs_cells%in%lc)]
        mur = muhat[which(obs_cells%in%rc)]
      }
      
      #variance
      if (!is.null(vsc)) {
        vl = vsc[which(obs_cells%in%lc)]
        vr = vsc[which(obs_cells%in%rc)]
      }
      
      #muhate = c(muhate, (mul + mur)/2)
      #lwre = c(lwre, (mul + mur)/2 - 2 * sqrt(vl))
      #uppe = c(uppe, (mul + mur)/2 + 2 * sqrt(vr))
      
      if (!is.null(muhat)) {
        muhatek = (mul - 2*sqrt(vl) + mur + 2*sqrt(vr))/2
        #check if muhatek follows the shape constraint
        lr = c(lc, rc)
        #sh=10 is for block ave
        #sh=9 is for block ord
        if (sh == 1 | sh == 10 | sh == 9) {
          if (length(lr) > 1) {
            bool = mul <= muhatek & mur >= muhatek
          } else if (length(lr) == 1) {
            if (at_min) {
              bool = mur >= muhatek
            } else if (at_max) {
              bool = mul <= muhatek 
            }
          }
        }
        if (sh == 2) {
          if (length(lr) > 1) {
            bool = mul >= muhatek & mur <= muhatek
          } else if (length(lr) == 1) {
            if (at_min) {
              bool = mur <= muhatek
            } else if (at_max) {
              bool = mul >= muhatek 
            }
          }
        }
        if (!bool) {
          if (!any(at_max) & !any(at_min)) {
             muhatek = (mul+mur)/2
          } else if (at_min) {
            muhatek = mul
          } else if (at_max)  {
            muhatek = mur
          }
        }
        muhate = c(muhate, muhatek)
      }
      
      if (!is.null(vsc)) {
        #varek = (mur + 2*sqrt(vr) - mul + 2*sqrt(vl))^2 / 16 
        rpt = mur + 2*sqrt(vr); lpt = mul - 2*sqrt(vl)
        varek = (rpt-lpt)^2/16
        vsce = c(vsce, varek)
        if (at_min & (ndek == 0)) {
          lwrek = -1e+5
          uppek = muhatek + 2*sqrt(varek)
        } else if (at_max & (ndek == 0)) {
          uppek = 1e+5
          lwrek = muhatek - 2*sqrt(varek)
        } else if (at_min & (ndek >= 1) & (ndek <= 10)) {
          varek_lwr = max(vsc)
          lwrek = muhatek - 2*sqrt(varek_lwr)
          uppek = muhatek + 2*sqrt(varek)
        } else if (at_max & (ndek >= 1) & (ndek <= 10)) {
          varek_upp = max(vsc)
          lwrek = muhatek - 2*sqrt(varek)
          uppek = muhatek + 2*sqrt(varek_upp)
        } else {
          lwrek = muhatek - 2*sqrt(varek)
          uppek = muhatek + 2*sqrt(varek)
        }
        lwre = c(lwre, lwrek)
        uppe = c(uppe, uppek)
      }
      
      ps = unique(c(lc, rc))
      ns = unique(c(n_lc, n_rc))
      pslst[[k]] = ps
      nslst[[k]] = ns
      domain.ord_em = c(domain.ord_em, dmem)
    }
    
    #yvec_all[obs_cells] = yvec
    #yvec_all[empty_cells] = ye
    
    #w_all[obs_cells] = w[obs_cells] 
    #w_all[obs_cells] = w
    #w_all[empty_cells] = we
    # 
    domain.ord_all[obs_cells] = domain.ord
    domain.ord_all[empty_cells] = domain.ord_em
    # 
    if (!is.null(muhat)) {
      muhat_all[obs_cells] = muhat
      muhat_all[empty_cells] = muhate
      muhat = muhat_all
    }
    
    if (!is.null(vsc)) {
      lwr_all[obs_cells] = lwr
      lwr_all[empty_cells] = lwre
      
      upp_all[obs_cells] = upp
      upp_all[empty_cells] = uppe
      
      upp = upp_all
      lwr = lwr_all
      
      vsc_all[obs_cells] = vsc
      vsc_all[empty_cells] = vsce
      vsc = vsc_all
    }
    
    #yvec = yvec_all
    #w = w_all
    domain.ord = domain.ord_all
    
  } else if (Nd >= 2) {
    vsc_all = lwr_all = upp_all = muhat_all = 1:M*0
    yvec_all = w_all = domain.ord_all = 1:M*0
    empty_grids = grid2[empty_cells, ]
    maxs = apply(grid2, 2, max)
    mins = apply(grid2, 2, min)
    at_max = (empty_grids == maxs)
    at_min = (empty_grids == mins)
    
    vsce = lwre = uppe = muhate = domain.ord_em = NULL
    ye = we = NULL
    nslst = pslst = vector('list', 2)
    nbors = NULL
    nbors_lst_0 = fn_nbors(empty_grids, grid2, mins, maxs)
    for(k in 1:nrow(empty_grids)){
      #new: if is new because of rm_id2
      nbors_k = nbors_lst_0[[k]]
      if (nrow(nbors_k) > 0) {
        ps = apply(nbors_k, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
        #temp
        ord = order(ps)
        ps2 = sort(ps)
        #print (ps)
        ns = nd[ps]
        
        pslst[[k]] = ps2
        nslst[[k]] = ns
        Ns = Nhat[ps]
        ys = yvec[which(obs_cells%in%ps)]
        
        ndek = nd_ne[k]
        #the way to choose l and r domainis are different in 2D case
        # # of neighbors may > 2
        row_id = as.numeric(rownames(grid2))
        D1 = Ds[1]; D2 = Ds[2]
        row_id = matrix(row_id, nrow=D1)
  
        ek = empty_cells[k]
        row_col_val_ek = which(row_id == ek, arr.ind=T)
        row_col_vals_nbors = NULL
        for(psi in ps){
          row_col_vals_nbors = rbind(row_col_vals_nbors, which(row_id == psi, arr.ind=T))
        }
        row_col_vals_nbors = as.data.frame(row_col_vals_nbors)
      
        d1_nbors = ps[which(row_col_vals_nbors$col == row_col_val_ek[2])]
        if (length(ps) > length(d1_nbors)) {
          d2_nbors = ps[!ps %in% d1_nbors]
        } else {d2_nbors = 0}
        #if (any(which(row_col_vals_nbors$row == row_col_val_ek[1]))){
        #  d2_nbors = ps[which(row_col_vals_nbors$row == row_col_val_ek[1])]
        #} else {
        #  d2_nbors = 0
        #}
        nbors_lst = vector('list', 2)
        nbors_lst[[1]] = d1_nbors
        nbors_lst[[2]] = d2_nbors
          
        mu_1d = muhat[which(obs_cells%in%d1_nbors)]
        mu_2d = muhat[which(obs_cells%in%d2_nbors)]
        #new:
        #if(sh[1] == 0){
        #  mu_1d = NULL
        #}
        #if(sh[2] == 0){
        #  mu_2d = NULL
        #}
        mu_12d = c(mu_1d, mu_2d)
        mu_max = mu_12d[which.max(mu_12d)]
        mu_min = mu_12d[which.min(mu_12d)]
        
        d12_nbors = c(d1_nbors, d2_nbors)
        max_ps = d12_nbors[which.max(mu_12d)]
        min_ps = d12_nbors[which.min(mu_12d)]
        
        v_max = vsc[which(obs_cells %in% max_ps)]
        v_min = vsc[which(obs_cells %in% min_ps)]
        
        muhatek = (mu_min - 2*sqrt(v_min) + mu_max + 2*sqrt(v_max))/2
        #print (paste('first ', muhatek))
        #check two shape constraints
        #at_maxk = ek == max(row_id)
        #at_mink = ek == min(row_id)
        at_maxk = c(row_col_val_ek[1] == D1, row_col_val_ek[2] == D2)
        at_mink = c(row_col_val_ek[1] == 1, row_col_val_ek[2] == 1)
          
        mul_2d = mur_2d = NULL
        row_ek = row_col_val_ek[1]; col_ek = row_col_val_ek[2]
        for(d in 1:2) {
          #print (c(k,d))
          #new: one dim has no constraint
          if (d == 1 & sh[d] == 0) {
            next 
          } else if (d == 2 & sh[d] == 0) {
            break
          }
          psd = nbors_lst[[d]]
          #if (length(psd) == 1)  {
          #  lps = rps = psd 
          #} else if (length(psd) > 1) {
            if (at_maxk[d]) {
              #lps = max(psd[psd < ek])
              lps = min(psd[psd < ek])
              #rps = lps
              rps = NULL
            }
            if (at_mink[d]) {
              #rps = min(psd[psd > ek])
              rps = max(psd[psd > ek])
              #lps = rps
              lps = NULL
            }
            if (!at_maxk[d] & !at_mink[d]) {
              lps = rps = NULL
              if (any(psd < ek)){
                lps = max(psd[psd < ek])
              }
              if (any(psd > ek)){
                rps = min(psd[psd > ek])
              }
            }
          #}
        
          if (!is.null(lps)) {
            mul = min(muhat[which(obs_cells%in%lps)])
          } else {mul = NULL}
          if (!is.null(rps)) {
            mur = max(muhat[which(obs_cells%in%rps)]) 
          } else {mur = NULL}
          
          lr = c(lps, rps)
          #temp for sh = 9
          if (sh[d] == 1 | sh[d] == 9) {
            if (length(lr) > 1) {
              bool = mul <= muhatek & mur >= muhatek
            } else if (length(lr) == 1) {
              if (d == 1) {
                if (row_ek == 1 | is.null(mul)) {
                  bool = mur >= muhatek
                } else if (row_ek == D1 | is.null(mur)) {
                  bool = mul <= muhatek 
                }
              }
              if (d == 2) {
                if (col_ek == 1) {
                  bool = mur >= muhatek
                } else if (col_ek == D2) {
                  bool = mul <= muhatek 
                }
              }
            }
          }
          if (sh[d] == 2) {
            if (length(lr) > 1) {
              bool = mul >= muhatek & mur <= muhatek
            } else if (length(lr) == 1) {
              if (d == 1) {
                if (row_ek == 1 | is.null(mul)) {
                  bool = mur <= muhatek
                } else if (row_ek == D1 | is.null(mur)) {
                  bool = mul >= muhatek 
                } 
              }
              if (d == 2) {
                if (col_ek == 1 | is.null(mul)) {
                  bool = mur <= muhatek
                } else if (col_ek == D2 | is.null(mur)) {
                  bool = mul >= muhatek 
                } 
              }
            }
          }
          #print (muhatek)
          if (!bool) {
            if (!is.null(mul) & !is.null(mur)) {
              # cat('mul:')
              # print (mul)
              # cat('mur:')
              # print (mur)
              # #will be used only when mul = mur? (most cases)
              muhatek = (mul+mur)/2
            } else if (!is.null(mul) & is.null(mur)) {
              muhatek = mul
            } else if (is.null(mul) & !is.null(mur))  {
              muhatek = mur
            }
          }
          #print (c(bool, muhatek))
        }
        
        muhate = c(muhate, muhatek)
        varek = (mu_max + 2*sqrt(v_max) - mu_min + 2*sqrt(v_min))^2 / 16 
        vsce = c(vsce, varek)
        lwrek = muhatek - 2*sqrt(varek)
        uppek = muhatek + 2*sqrt(varek)
        lwre = c(lwre, lwrek)
        uppe = c(uppe, uppek)
        
        nbors = rbind(nbors, nbors_k)
      } else {
        #temp: if there's no observed neighbor, use the previous imputation
        #new: zeros_ps: 1  2  3, won't find pslst[[k-1]] for 1
        #pslst[[k]] = pslst[[k-1]]
        #nslst[[k]] = nslst[[k-1]]
        ye = c(ye, rev(ye)[1])
        we = c(we, rev(we)[1])
        
        if (!is.null(muhat)) {
          muhate = c(muhate, rev(muhate)[1])
        }
        if (!is.null(vsc)) {
          vsce = c(vsce, rev(vsce)[1])
          lwre = c(lwre, rev(lwre)[1])
          uppe = c(uppe, rev(uppe)[1])
        }
      }
      #print (muhate)
    }
    #yvec_all[obs_cells] = yvec
    #yvec_all[empty_cells] = ye
    #w_all[obs_cells] = w[obs_cells] 
    #w_all[obs_cells] = w
    #w_all[empty_cells] = we
    #yvec = yvec_all
    #w = w_all
    
    domain.ord_all[obs_cells] = domain.ord
    #temp:
    domain.ord_em = obs_cells[which(obs_cells%in%ps[1:length(empty_cells)])]
    domain.ord_all[empty_cells] = domain.ord_em
    domain.ord = domain.ord_all
    
    if (!is.null(muhat)) {
      muhat_all[obs_cells] = muhat
      muhat_all[empty_cells] = muhate
      muhat = muhat_all
    }
    
    if (!is.null(vsc)){
      lwr_all[obs_cells] = lwr
      lwr_all[empty_cells] = lwre
      upp_all[obs_cells] = upp
      upp_all[empty_cells] = uppe
      vsc_all[obs_cells] = vsc
      vsc_all[empty_cells] = vsce
      
      upp = upp_all
      lwr = lwr_all
      vsc = vsc_all
    }
  }
  
  rslt = list(ps = ps, ns = ns, muhat = muhat, lwr = lwr, upp = upp, domain.ord = domain.ord, vsc = vsc)
  #rslt = list(ps = ps, ns = ns, muhat = muhat, lwr = lwr, upp = upp, domain.ord = domain.ord, yvec = yvec, w = w)
  return(rslt)
}


###############################
#subrountine to find nbors (>=2D)
###############################
fn_nbors = function(empty_grids, grid2, mins, maxs){
  nr = nrow(empty_grids)
  nb_lst = vector("list", length = nr)
  to_merge = list()
  iter = 1
  for(k in 1:nr){
    #ndek = nd_ne[k]
    #axes = apply(empty_grids[k,], 2, function(e) e + c(-1,1))
    nbors_k = NULL
    ptk = empty_grids[k, ]
    ptk = as.numeric(ptk)
    nc = ncol(grid2)
    for(i in 1:nc){
      pt1 = ptk[i]
      pt2 = ptk[-i] 
      axes = pt1 + c(-1,1) 
      pair_i = matrix(0, nrow=2, ncol=nc)
      pair_i[1, i] = axes[1]; pair_i[1, -i] = pt2
      pair_i[2, i] = axes[2]; pair_i[2, -i] = pt2
      nbors_k = rbind(nbors_k, pair_i)
    }
    rm_id = unique(c(which(apply(nbors_k, 1, function(e) any(e<mins))), which(apply(nbors_k, 1, function(e) any(e>maxs)))))
    if (length(rm_id) > 0) {
      nbors_k = nbors_k[-rm_id, , drop = FALSE]
    }
    #new: label the cells whose neighbors should be merged later
    # rm_id2 = NULL
    if (nrow(nbors_k) > 0) {
      for(h in 1:nrow(nbors_k)) {
        nbh = nbors_k[h, ]
        bool = any(apply(empty_grids, 1, function(e) all(e == nbh)))
        if (bool) {
          #rm_id2 = c(rm_id2, h) 
          ps = as.numeric(which(apply(empty_grids, 1, function(e) all(e == nbh))))
          to_merge[[iter]] = sort(c(k, ps))
          iter = iter + 1
        }
      }
    }
    nb_lst[[k]] = nbors_k
  }
  
  #need to be more efficient
  to_merge = unique(to_merge)
  nm = length(to_merge)
  if (nm > 1) {
    for(i in 1:(nm-1)) {
      to_merge_i = to_merge[[i]]
      vec = to_merge_i
      vec_ps = i
      for(j in (i+1):nm) {
        to_merge_j = to_merge[[j]]
        if (any(to_merge_j %in% to_merge_i)) {
          vec = sort(unique(c(vec, to_merge_j)))
          vec_ps = c(vec_ps, j)
        }
      }
      if (length(vec_ps) > 1) {
        to_merge[vec_ps] = lapply(to_merge[vec_ps], function(x) x = vec)
      }
    }
    to_merge = unique(to_merge)
  }
  
  if (nm > 0) {
    for(i in 1:length(to_merge)) {
      ps = to_merge[[i]]
      nbors_ps = unique(do.call(rbind, nb_lst[ps]))
      rm_id2 = NULL
      for(h in 1:nrow(nbors_ps)) {
        nbh = nbors_ps[h, ]
        bool = any(apply(empty_grids, 1, function(e) all(e == nbh)))
        if (bool) {
          rm_id2 = c(rm_id2, h) 
        }
      }
      if (length(rm_id2) > 0) {
        nbors_ps = nbors_ps[-rm_id2,,drop=FALSE]
      }
      nb_lst[ps] = lapply(nb_lst[ps], function(x) x = nbors_ps)
    }
  }
  return (nb_lst)
}

####################################################################################
#compute the variance for n=2...10 cells
#keep muhat, only compute weighted variance
#n=0 and 1: 100%; n=2,3 90%; n=4,5 70%; n=6,7, 50%; n=8,9 30% n=10, 10%
#Nd=2 only, Nd=1 coverage rate is good?
####################################################################################
impute_em3 = function(small_cells, M, yvec, Nd, nd, Nhat, w, domain.ord, grid2, 
                      Ds=NULL, sh=NULL, muhat=NULL, lwr=NULL, upp=NULL, vsc=NULL, new_obs_cells=NULL)
{
  ne = length(small_cells)
  #observed sample size for n=2...10 cells
  #the `small_cells` are not empty
  nd_ne = nd[small_cells]
  vsc_ne = vsc[which(new_obs_cells%in%small_cells)]
  
  lc = rc = 1:ne*0
  if (Nd >= 2) {
    vsc_all = lwr_all = upp_all = muhat_all = 1:M*0
    yvec_all = w_all = domain.ord_all = 1:M*0
    empty_grids = grid2[small_cells, ]
    maxs = apply(grid2, 2, max)
    mins = apply(grid2, 2, min)
    at_max = (empty_grids == maxs)
    at_min = (empty_grids == mins)
    
    vsce = lwre = uppe = muhate = domain.ord_em = NULL
    ye = we = NULL
    nslst = pslst = vector('list', 2)
    nbors = NULL
    nbors_lst_0 = fn_nbors(empty_grids, grid2, mins, maxs)
    for(k in 1:nrow(empty_grids)){
      #new: if is new because of rm_id2
      nbors_k = nbors_lst_0[[k]]
      if (nrow(nbors_k) > 0) {
        ps = apply(nbors_k, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
        #temp
        ord = order(ps)
        ps2 = sort(ps)
        #print (ps)
        ns = nd[ps]
        
        pslst[[k]] = ps2
        nslst[[k]] = ns
        Ns = Nhat[ps]
        ys = yvec[which(new_obs_cells%in%ps)]
        
        ndek = nd_ne[k]
        #the way to choose l and r domainis are different in 2D case
        # # of neighbors may > 2
        row_id = as.numeric(rownames(grid2))
        D1 = Ds[1]; D2 = Ds[2]
        row_id = matrix(row_id, nrow=D1)
        
        ek = small_cells[k]
        row_col_val_ek = which(row_id == ek, arr.ind=T)
        row_col_vals_nbors = NULL
        for(psi in ps){
          row_col_vals_nbors = rbind(row_col_vals_nbors, which(row_id == psi, arr.ind=T))
        }
        row_col_vals_nbors = as.data.frame(row_col_vals_nbors)
        
        d1_nbors = ps[which(row_col_vals_nbors$col == row_col_val_ek[2])]
        if (length(ps) > length(d1_nbors)) {
          d2_nbors = ps[!ps %in% d1_nbors]
        } else {d2_nbors = 0}
        #if (any(which(row_col_vals_nbors$row == row_col_val_ek[1]))){
        #  d2_nbors = ps[which(row_col_vals_nbors$row == row_col_val_ek[1])]
        #} else {
        #  d2_nbors = 0
        #}
        nbors_lst = vector('list', 2)
        nbors_lst[[1]] = d1_nbors
        nbors_lst[[2]] = d2_nbors
        
        mu_1d = muhat[which(new_obs_cells%in%d1_nbors)]
        mu_2d = muhat[which(new_obs_cells%in%d2_nbors)]
        mu_12d = c(mu_1d, mu_2d)
        mu_max = mu_12d[which.max(mu_12d)]
        mu_min = mu_12d[which.min(mu_12d)]
        
        d12_nbors = c(d1_nbors, d2_nbors)
        max_ps = d12_nbors[which.max(mu_12d)]
        min_ps = d12_nbors[which.min(mu_12d)]
        
        v_max = vsc[which(new_obs_cells %in% max_ps)]
        v_min = vsc[which(new_obs_cells %in% min_ps)]
        
        muhatek = (mu_min - 2*sqrt(v_min) + mu_max + 2*sqrt(v_max))/2
        #check two shape constraints
        #at_maxk = ek == max(row_id)
        #at_mink = ek == min(row_id)
        at_maxk = c(row_col_val_ek[1] == D1, row_col_val_ek[2] == D2)
        at_mink = c(row_col_val_ek[1] == 1, row_col_val_ek[2] == 1)
        
        muhate = c(muhate, muhatek)
        #new:
        varek_new = (mu_max + 2*sqrt(v_max) - mu_min + 2*sqrt(v_min))^2 / 16 
        #the variance by survey
        #the if else can be replaced by vectorized computation?
        vsc_nek = vsc_ne[k]
        if(ndek >= 1 & ndek <= 3){
          #varek = varek_new*.5 + vsc_nek*.5 (not good)
          #varek = varek_new*.8 + vsc_nek*.2
          #varek = varek_new*.9 + vsc_nek*.1
          varek = varek_new
        }else if(ndek == 4 | ndek == 5){
          #varek = varek_new*.5 + vsc_nek*.5 (not good)
          #varek = varek_new*.6 + vsc_nek*.4
          #varek = varek_new*.7 + vsc_nek*.3
          varek = varek_new
        }else if(ndek == 6 | ndek == 7){
          #varek = varek_new*.5 + vsc_nek*.5 (not good)
          #varek = varek_new*.4 + vsc_nek*.6
          #varek = varek_new*.5 + vsc_nek*.5
          varek = varek_new
        }else if(ndek == 8 | ndek == 9){
          #varek = varek_new*.2 + vsc_nek*.8
          #varek = varek_new*.1 + vsc_nek*.9
          #varek = varek_new*.3 + vsc_nek*.7
          varek = varek_new
        }else if(ndek == 10){
          #varek = varek_new*.1 + vsc_nek*.9
          #varek = varek_new*.1 + vsc_nek*.9
          varek = varek_new
        }
        
        lwrek = muhatek - 2*sqrt(varek)
        uppek = muhatek + 2*sqrt(varek)
        lwre = c(lwre, lwrek)
        uppe = c(uppe, uppek)
        vsce = c(vsce, varek_new)
        
        nbors = rbind(nbors, nbors_k)
      } else {
        #temp: if there's no observed neighbor, use the previous imputation
        #new: zeros_ps: 1  2  3, won't find pslst[[k-1]] for 1
        #pslst[[k]] = pslst[[k-1]]
        #nslst[[k]] = nslst[[k-1]]
        ye = c(ye, rev(ye)[1])
        we = c(we, rev(we)[1])
        
        if (!is.null(muhat)) {
          muhate = c(muhate, rev(muhate)[1])
        }
        if (!is.null(vsc)) {
          lwre = c(lwre, rev(lwre)[1])
          uppe = c(uppe, rev(uppe)[1])
          vsce = c(vsce, rev(vsce)[1])
        }
      }
      #print (c(k, length(lwre)))
    }
    
    # domain.ord_all[new_obs_cells] = domain.ord
    # #temp:
    # domain.ord_em = new_obs_cells[which(new_obs_cells%in%ps[1:length(small_cells)])]
    # domain.ord_all[small_cells] = domain.ord_em
    # domain.ord = domain.ord_all
    # 
    # # if (!is.null(muhat)) {
    #   muhat_all[new_obs_cells] = muhat
    #   muhat_all[small_cells] = muhate
    #   muhat = muhat_all
    # }
    # 
    if (!is.null(vsc)){
      lwr_all[new_obs_cells] = lwr
      #add the case when all cells are small cells
      if(!is.null(lwre)){
        lwr_all[small_cells] = lwre
        upp_all[small_cells] = uppe
        vsc_all[small_cells] = vsce
      }
      upp_all[new_obs_cells] = upp
      vsc_all[new_obs_cells] = vsc
 
      vsc = vsc_all  
      upp = upp_all
      lwr = lwr_all
    }
  }
  
  rslt = list(lwr = lwr, upp = upp, vsc = vsc)
  return(rslt)
}
