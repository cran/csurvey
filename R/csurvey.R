#############################
## main routine for the user#
#############################
csvy = function(formula, design, subset=NULL, nD=NULL, family=stats::gaussian(),
                 amat=NULL, level=0.95, n.mix=100L, test=TRUE,...) {
  cl = match.call()
  subset = substitute(subset)
  subset = eval(subset, model.frame(design), parent.frame())
  if (!is.null(subset)){
    if (any(is.na(subset)))
      stop("subset must not contain NA values")
    design = design[subset,]
  }
  if (is.character(family))
    family = get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family = family()
  if (is.null(family$family))
    stop("'family' not recognized!")
  labels = NULL
  mf = match.call(expand.dots = FALSE)
  m = match("formula", names(mf), 0L)
  mf = mf[c(1L, m)]
  mf[[1L]] = as.name("model.frame")
  mf = eval(mf, model.frame(design), parent.frame())
  ynm = names(mf)[1]
  mt = attr(mf, "terms")
  y = model.response(mf, "any")
  shapes_add = NULL
  xmat_add0 = NULL
  xmat_add = NULL
  xnms_add = NULL
  nums_add = NULL
  xid_add = 1
  relaxes = NULL
  block.ave.lst = vector("list", length = (ncol(mf)-1))
  block.ord.lst = vector("list", length = (ncol(mf)-1))
  cnms = colnames(mf)
  for (i in 2:ncol(mf)) {
    if (is.numeric(attributes(mf[,i])$shape)) {
      labels = c(labels, "additive")
      shapes_add = c(shapes_add, attributes(mf[,i])$shape)
      reli = attributes(mf[,i])$relax
      if (is.null(reli)) {
        reli = FALSE
      } else {
        reli = TRUE
      }
      relaxes = c(relaxes, reli)
      if (is.factor(mf[,i]) && is.character(levels(mf[,i]))) {
        xmat_add = cbind(xmat_add, as.numeric(mf[,i]))
      } else {
        xmat_add = cbind(xmat_add, mf[,i])
      }
      xmat_add0 = cbind(xmat_add0, mf[,i])
      xnms_add = c(xnms_add, attributes(mf[,i])$nm)
      xid_add = xid_add + 1
      
      block.ord.lst[[i-1]] = -1
      block.ave.lst[[i-1]] = -1
      if ((attributes(mf[,i])$categ == "block.ord")) {
        block.ord.lst[[i-1]] = attributes(mf[,i])$order
      }
      if ((attributes(mf[,i])$categ == "block.ave")) {
        block.ave.lst[[i-1]] = attributes(mf[,i])$order
      }
    } else if (is.null(attributes(mf[,i])$shape)) {
      block.ord.lst[[i-1]] = -1
      block.ave.lst[[i-1]] = -1
      reli = attributes(mf[,i])$relax
      if (is.null(reli)) {
        reli = FALSE
      } else {
        reli = TRUE
      }
      relaxes = c(relaxes, reli)
      shapes_add = c(shapes_add, 0)
      labels = c(labels, "additive")
      xmat_add = cbind(xmat_add, mf[,i])
      xmat_add0 = cbind(xmat_add0, mf[,i])
      xnms_add = c(xnms_add, cnms[i])
      xid_add = xid_add + 1
    }
  }
  
  xid_add = 2:ncol(mf)
  #xmat_add = mf[, -1, drop = FALSE]
  #how to handle ...?
  #fit = csvy.fit(design = design, fo = formula, subset = subset, family = family, M = nD, xm = xmat_add, xnms_add = xnms_add,
  #                sh = shapes_add, ynm = ynm, block.ave.lst = block.ave.lst, block.ord.lst = block.ord.lst,
  #                level = level, n.mix = n.mix, test = test,...)
  fit = eval(call("csvy.fit", design = design, family = family, M = nD, 
                   relaxes = relaxes, xm = xmat_add, xnms_add = xnms_add,
                   sh = shapes_add, ynm = ynm, amat = amat,
                   block.ave.lst = block.ave.lst, 
                   block.ord.lst = block.ord.lst, 
                   level = level, n.mix = n.mix, cl = cl, test = test))
  fit$survey.design = design
  fit$data = design$variables
  fit$na.action = attr(mf, "na.action")
  fit$call = cl
  fit$family = family
  fit$n.mix = n.mix
  fit$terms = mt
  fit$xmat_add = xmat_add 
  fit$xmat_add0 = xmat_add0
  fit$xnms_add = xnms_add
  fit$shapes_add = shapes_add
  fit$ynm = ynm
  class(fit) = c("csvy", "cgam", "svyglm", "glm")
  # classes = if (inherits(design, "svyrep.design")) {
  #   c("svrepglm", "svyglm")
  # } else {
  #   "svyglm"
  # }
  # structure(c(fit, list(call = cl, n.mix = n.mix, xmat_add = xmat_add, xmat_add0 = xmat_add0, xnms_add = xnms_add, na.action = attr(mf, "na.action"), 
  #                       shapes_add = shapes_add, ynm = ynm, family = family, terms = mt, survey.design = design, data = design$variables)),
  #           class = c("csvy", "cgam", classes))
  return (fit)
}

############
## csvy.fit#
############
csvy.fit = function(design, family=stats::gaussian(), M=NULL, relaxes=NULL, 
                    xm=NULL, xnms_add=NULL, sh=1, ynm=NULL, nd=NULL, 
                    ID=NULL, block.ave.lst=NULL, block.ord.lst=NULL, 
                    amat=NULL, level=0.95, n.mix=100L, cl=NULL, test=TRUE,...){
  #wp is population wt
  #Ds: observed domains?
  #make sure each column in xm is numeric:
  #xm = map_dbl(xm, .f = function(.x) as.numeric(.x))
  bool = apply(xm, 2, function(x) is.numeric(x))
  if(any(!bool)){
    xm = apply(xm, 2, function(x) as.numeric(x))
  }
  Ds = apply(xm, 2, function(x) length(unique(x)))
  xm.red = unique(xm)
  xvs2 = apply(xm.red, 2, function(x) sort(unique(x)))
  if(is.matrix(xvs2)){
    grid2 = expand.grid(as.data.frame(xvs2))
  }else if(is.list(xvs2)){
    grid2 = expand.grid(xvs2)
  }
  colnames(grid2) = xnms_add
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
  
  Tot_obs = length(obs_cells)
  #print (obs_cells)
  #print (Max_obs)
  if (is.null(M)) {
    #stop("nD (total number of domains) is not provided! \n")
    M = max(obs_cells)
  } else if (M %% 1 != 0 | M <= 1) {
    stop(paste("nD (total number of domains) must be a positive integer > 1!", "\n"))
    #M = max(obs_cells)
  } else if (M < Tot_obs) {
    stop(paste("There are at least", Tot_obs, "domains! nD (total number of domains) should be >= ",
               Tot_obs, sep = " ", "\n"))
  }
  
  #new: to be used in makeamat_2d_block
  M_0 = M
  #if (any(class(ds) == "survey.design")) {
  if (inherits(ds, "survey.design")) {
    weights = 1/(ds$prob)
  }
  #if (any(class(ds) == "svyrep.design")) {
  if (inherits(ds, "svyrep.design")) {
    weights = ds$pweights
  }
  #Nhat and stratawt will follow the order: 1 to M
  y = df[[ynm]]
  #test:
  if(family$family %in% c("binomial", "quasibinomial")){
    uvals = unique(y)
    if(length(uvals) > 2) {
      stop(paste("The response has more than 2 outcomes! family = binomial won't work!", "\n"))
    } else {
      if(is.factor(y)){
        ulvls = levels(y)
        success = ulvls[2]
        y = ifelse(y == success, 1, 0)
      } 
    }
  }
  n = length(y)
  ys = stratawt = vector("list", length=M)
  #add ybar_vec for cic of non-gaussian cases
  Nhat = nd = ybar_vec = 1:M*0
  uds_num = as.numeric(levels(uds))[uds] 
  #new: for binomial
  y_ch_id = NULL
  
  #check more!
  # if (length(Ds) == 1) {
  #   for(i in 1:M){
  #     udi = uds_num[i]
  #     ps = which(ID %in% udi)
  #     ysi = y[ps]
  #     wi = weights[ps]
  #     Nhi = sum(wi)
  # 
  #     Nhat[i] = Nhi
  #     nd[i] = length(wi)
  #     ys[[i]] = ysi
  #     #new: for binomial
  #     if (all(ysi == 0) | all(ysi == 1)) {
  #       y_ch_id = c(y_ch_id, i)
  #     }
  #     stratawt[[i]] = wi
  #     ybar_vec[i] = mean(ysi)
  #   }
  # }
  
  if (length(Ds) >= 1) {
    for(udi in uds_num){
      ps = which(ID %in% udi)
      ysi = y[ps]
      wi = weights[ps]
      Nhi = sum(wi)
      
      #check more
      #Nhat[which(1:M %in% udi)] = Nhi
      #nd[which(1:M %in% udi)] = length(wi)
      #ys[[which(1:M %in% udi)]] = ysi
      #print (which(1:M %in% udi)) #cannot be used in 1D case if uds_num not start from 1
      #print (which(uds_num %in% udi)) #cannot be used, cannot handle empty cells
      #print (udi)
      #test!
      if (min(uds_num) > 1) {
        Nhat[udi - min(uds_num) + 1] = Nhi
        nd[udi - min(uds_num) + 1] = length(wi)
        ys[[udi - min(uds_num) + 1]] = ysi
        #new: for binomial
        if (all(ysi == 0) | all(ysi == 1)) {
          y_ch_id = c(y_ch_id, udi- min(uds_num) + 1)
        }
        stratawt[[udi - min(uds_num) + 1]] = wi
        ybar_vec[udi - min(uds_num) + 1] = mean(ysi)
      } else {
        Nhat[udi] = Nhi
        nd[udi] = length(wi)
        ys[[udi]] = ysi
        #new: for binomial
        if (all(ysi == 0) | all(ysi == 1)) {
          #stop('check!')
          y_ch_id = c(y_ch_id, udi)
          #y_ch_id = c(y_ch_id, which(uds_num %in% udi))
        }
        
        #stratawt[[i]] = wi
        #ybar_vec[i] = mean(ysi)
        #stratawt[[which(1:M %in% udi)]] = wi
        #ybar_vec[which(1:M %in% udi)] = mean(ysi)
        stratawt[[udi]] = wi
        ybar_vec[udi] = mean(ysi)
      }
    }
  }
  
  #new: for binomial
  y_ch_id2 = NULL
  log_odds = NULL
  
  if (!is.null(y_ch_id) & family$family %in% c('binomial', 'quasibinomial')){
    n_ch_id = length(y_ch_id)
    y_ch_id2 = y_ch_id - 1
    
    #check!
    if (y_ch_id2[1] == 0) {
      nbor = 1
      while(all(ys[[nbor]] == 1) | all(ys[[nbor]] == 0) | is.null(ys[[nbor]])){
        nbor = nbor + 1
      }
      y_ch_id2[1] = nbor
    }
    
    ch = sapply(y_ch_id2, function(e) all(ys[[e]] == 1) | all(ys[[e]] == 0) | is.null(ys[[e]]))
    
    if (any(ch)) {
      ps = which(ch)
      for(i in 1:length(ps)) {
        psi = ps[i]
        nbor = y_ch_id2[psi]
        if (nbor == 1) {
          while(nbor %in% y_ch_id | is.null(ys[[nbor]])){
            nbor = nbor + 1
          }
        } else {
          while(nbor %in% y_ch_id | is.null(ys[[nbor]])) {
            nbor = nbor - 1
            #test more
            if (nbor == 1) {
              while(nbor %in% y_ch_id | is.null(ys[[nbor]])){
                nbor = nbor + 1
              }
            }
          }
        }
        y_ch_id2[psi] = nbor
      }
    }
    
    ###
    y_merge = rbind(y_ch_id2, y_ch_id)
    n_y_merge = ncol(y_merge)
    for(i in 1:n_y_merge){
      y_merge_i = c(y_merge[, i])
      den_i = 0
      num_i = 0
      for(j in y_merge_i){
        #print (j)
        ys_j = ys[[j]]
        den_i = length(ys_j) + den_i
        num_i = sum(ys_j) + num_i
      }
      p_i = num_i / den_i
      log_odds = c(log_odds, log(p_i / (1 - p_i)))
    }
  }
  
  #new:
  obs_cells = which(nd > 0)
  fo = as.formula(paste0("~", ynm))
  #use fo2 for svyglm
  xnms_add_fac = sapply(xnms_add, function(xnmi) paste0("factor(", xnmi, ")"))
  #xnms_add_fac = xnms_add
  #check more!
  #xnms_fo_fac = paste(xnms_add_fac, collapse = "+")
  xnms_fo_fac = paste(xnms_add_fac, collapse = "*")
  fo2 = formula(paste(ynm, "~", xnms_fo_fac))
  #new:
  fo2_null = formula(paste(ynm, "~", 1))
  xnms_fo = paste(xnms_add, collapse = "+")
  #new:
  if(family$family == "gaussian"){
    ans.unc = suppressWarnings(svyby(formula=fo, by=~ID, design=ds, FUN=svymean, 
                                     covmat=TRUE, multicore=FALSE, drop.empty.groups=FALSE))
    v1 = vcov(ans.unc)
    etahatu = yvecu = yvec = ans.unc[[ynm]]
    vsu = round(diag(v1), 5)
    ans.unc_null = svyglm(formula=fo2_null, design=ds, family=family)
    prior.weights = ans.unc_null$prior.weights
    ans.unc_cp_0 = ans.unc
  } else {
    fo_by = if (length(xnms_add) == 1) formula(paste("~", xnms_add)) else formula(paste("~", paste(xnms_add, collapse = "+")))
    ans.unc_cp_0 = suppressWarnings(svyby(formula = fo, by = fo_by, design = ds, FUN = svymean, keep.var = TRUE,
                                           keep.names = TRUE, verbose = FALSE, vartype = "se",
                                           drop.empty.groups = FALSE, covmat = TRUE,
                                           na.rm.by = FALSE, na.rm.all = FALSE))
    #ans.unc = svyglm(formula=fo2, design=ds, family=family, 
    #                 control = list(epsilon = 1e-5, maxit = 12))
    
    #ans.unc = tryCatch({svyglm(formula=fo2, design=ds, family=family, 
    #                           control = list(epsilon = 1e-7, maxit = 10))}, error = function(e) {e})
    #if(inherits(ans.unc, 'error')) {
    #  stop('check!')
    #} 
    
    ans.unc = suppressWarnings(svyglm(formula=fo2, design=ds, family=family, 
                                      control = list(epsilon = 1e-8, maxit = 10)))
    
    prior.weights = ans.unc$prior.weights
    #create newdata
    newdata = grid2
    empty_cells = which(nd == 0)
    # #predict.svyglm will ignore empty cells for 1D?
    if(length(Ds) > 1 & length(empty_cells) > 0) {
      newdata = newdata[-empty_cells, ,drop=FALSE]
    }
    
    if(length(empty_cells) == 0){
      p.unc = predict(ans.unc, newdata, se.fit = TRUE, type = 'link', vcov = TRUE)
      v1 = vcov(p.unc)
      etahatu = yvecu = yvec = as.data.frame(p.unc)$link
      vsu = round(diag(v1), 5)
    }
    
    if(length(empty_cells) > 0){
      tt = delete.response(terms(formula(ans.unc)))
      mf = model.frame(tt, data = newdata)
      mm = model.matrix(tt, mf)
      #modified from predict.svrepglm
      mm = mm[, -empty_cells]
      vv = mm %*% vcov(ans.unc) %*% t(mm)
      v1 = vv
      etahatu = yvecu = yvec = NULL
      etahatu = drop(mm %*% coef(ans.unc))
      yvecu = etahatu
      yvec = etahatu
      vsu = round(diag(v1), 5)
    }
    
    #temp fix:
    # bool1 = any(yvec > 5)
    # if (bool1) {
    #   ps1 = which(yvec > 5)
    #   etahatu[ps1] = yvecu[ps1] = yvec[ps1] = 5
    # }
    # bool2 = any(yvec < -5)
    # if (bool2) {
    #   ps2 = which(yvec < -5)
    #   etahatu[ps2] = yvecu[ps2] = yvec[ps2] = -5
    # }
    # 
  }
  
  z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
  hlu = z.mult*sqrt(vsu)
  lwru = yvecu - hlu
  uppu = yvecu + hlu
  #new: replace n=1 variance with the average
  #not to change vsu, which will be used for 
  ps = which(round(diag(v1), 7) == 0L)
  diag(v1)[ps] = mean(diag(v1)[-ps])
  
  w = Nhat/sum(Nhat)
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
  
  if(any(bool_ord == 0)){
    noord_cells = map(block.ord.lst, .f = function(.x) which(.x == 0))
    if(length(Ds) == 1){
      noord_cells = noord_cells[[1]]
    }
  }
  ne = length(empty_cells)
  nd_ne = nd[empty_cells]
  small_cells = which(nd >= 1 & nd <= 10)
  nsm = length(small_cells)
  zeros = ones = 0
  zeros_ps = NULL
  
  if(ne >= 1){
    if (any(nd == 0)) {
      zeros_ps = which(nd == 0)
      zeros = length(zeros_ps)
    }
    #zeros = ne - ones
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
  amat_0 = NULL #to be used in Nd > 1
  amat_ci = NULL #to be used for ci
  
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
                          zeros_ps=zeros_ps, noord_cells=noord_cells, xvs2=xvs2, grid2=grid2)
          amat_0 = makeamat(x=1:M_0, sh=sh, Ds=M_0, interp=TRUE, relaxes=FALSE, 
                            block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]], 
                            zeros_ps=NULL, noord_cells=noord_cells, xvs2=xvs2, grid2=grid2)
        }  
      } else if (sh == 0 & is.null(amat)) {
        #print (sh)
        stop ("User must define a constraint matrix!")
      }
    } else if (Nd >= 2) {
      #temp:
      if (any(sh > 0)) {
        #if (!relaxes) {
        #print (zeros_ps)
        #new: create another amat used for forcing monotonicity on ci
        amat_rslt = makeamat_2D(x=NULL, sh=sh, grid2=grid2, xvs2=xvs2,
                                zeros_ps=zeros_ps, noord_cells=noord_cells, 
                                Ds=Ds, interp=TRUE, 
                                block.ave.lst=block.ave.lst,
                                block.ord.lst=block.ord.lst, M_0=M_0)
        amat = amat_rslt$amat
        #amat_ci = amat_rslt$amat_ci #wrong, should include empty domains
        #to be used in impute_2 or impute_3
        amat_0_rslt = makeamat_2D(x=NULL, sh=sh, grid2=grid2, xvs2=xvs2,
                                  zeros_ps=NULL, noord_cells=noord_cells,
                                  Ds=Ds, interp=TRUE,
                                  block.ave.lst=block.ave.lst,
                                  block.ord.lst=block.ord.lst, M_0=M_0)
        amat_0 = amat_0_rslt$amat
        amat_ci = amat_0_rslt$amat_ci
      } else if (all(sh == 0) & is.null(amat)) {
        stop (paste("User must define a constraint matrix!", "\n"))
      }
    }
  }
  
  #new: for binomial, to avoid extreme etahatun
  if (ne == 0 & family$family %in% c("binomial", "quasibinomial")){
    if(!is.null(log_odds)){
      iter = 1
      for(i in 1:n_y_merge){
        y_merge_i = c(y_merge[, i])
        for(j in 1:length(y_merge_i)){
          ps_j = y_merge_i[j]
          yvec[ps_j] = log_odds[iter]
        }
        iter = iter + 1
      }
    }
    etahatu = yvecu = yvec
  } 
  
  if (ne >= 1 & family$family %in% c("binomial", "quasibinomial")){
    yvec_all = 1:M_0*0
    yvec_all[-empty_cells] = yvec
    if(!is.null(log_odds)){
      iter = 1
      for(i in 1:n_y_merge){
        y_merge_i = c(y_merge[, i])
        for(j in 1:length(y_merge_i)){
          ps_j = y_merge_i[j]
          yvec_all[ps_j] = log_odds[iter]
        }
        iter = iter + 1
      }
    }
    yvec = yvec_all[-empty_cells]
    etahatu = yvecu = yvec
  }
  
  #new: check irreducibility of amat
  if (!is.null(amat)) {
    #new:
    #if (family$family == "gaussian"){
    ans.polar = coneA(yvec, amat, w=w, msg=FALSE)
    etahat = round(ans.polar$thetahat, 10)
    #ans.polar = tryCatch({coneA(yvec, amat, w=w, msg=FALSE)}, error = function(e) {e})
    #if(inherits(ans.polar, 'error')) {
    #  stop('check!')
    #} else {
    #  etahat = round(ans.polar$thetahat, 10)
    #}
    #ans.polar = coneA(yvec, amat, w=w, msg=FALSE)
    #etahat = round(ans.polar$thetahat, 10)
    # } else {
    #   # umat = chol(v1)
    #   # uinv = solve(umat)
    #   # zvec = uinv %*% yvec
    #   # atil = amat %*% umat
    #   # ans.polar = coneA(zvec, atil, msg = FALSE)
    #   # phihat = round(ans.polar$thetahat, 10)
    #   # etahat = umat %*% phihat
    # 
    #   eans = eigen(v1, symmetric=TRUE)
    #   evecs = eans$vectors; evecs = round(evecs, 8)
    #   evals = eans$values
    #   sm = 1e-4
    #   if(any(evals < sm)){
    #     evals[which(evals < sm)] = sm
    #   }
    #   umat = evecs %*% diag(sqrt(evals))
    #   uinv = solve(umat) #solve will give singular error message when n is small
    #   atil = amat %*% umat
    #   zvec = uinv %*% yvec
    #   ans.polar = coneA(zvec, atil, msg=FALSE)
    #   phihat = round(ans.polar$thetahat, 10)
    #   etahat = umat %*% phihat
    # }
    
    #new
    #muhat = etahat 
    face = ans.polar$face
    #new: to be used in summary and anova
    edf = 1.5*ans.polar$df
    if (length(face) == 0){
      mat = wt %*% v1
    } else {
      dd = amat[face, ,drop=FALSE]
      wtinvdd = wtinv %*% t(dd)
      pmat_is = wtinvdd %*% solve(dd %*% wtinvdd) %*% dd
      mat = wt %*% (imat - pmat_is) %*% v1
    }
    #evec = yvec - muhat
    if(family$family == "gaussian"){
      evec = yvec - etahat
      CIC = t(evec) %*% wt %*% evec + 2*(sum(diag(mat)))
      CIC = as.numeric(CIC)
      #new: add CIC.un
      CIC.un = 2*(sum(diag(wt %*% v1)))
      CIC.un = as.numeric(CIC.un)
    }else{
      #check for empty cells
      evec = ybar_vec[which(nd != 0)] - c(family$linkinv(etahat))
      CIC = t(evec) %*% wt %*% evec + 2*(sum(diag(mat)))
      CIC = as.numeric(CIC)
      CIC.un = 2*(sum(diag(wt %*% v1)))
      CIC.un = as.numeric(CIC.un)
    }
  }
  
  acov = v1
  #if(is.integer(n.mix) & n.mix > 0){
  if(n.mix > 0){
    #get the constrained variance-covariance matrix
    lwr = upp = NULL
    acov = v1
    dp = -t(amat)
    dp = apply(dp, 2, function(e) e / (sum(e^2))^(.5))
    
    m_acc = M 
    sector = NULL
    times = NULL
    df.face = NULL
    iter = 1
    obs = 1:M
    #binomial can use mvrnorm?
    ysims = MASS::mvrnorm(n.mix, mu=etahat, Sigma=v1)
    for (iloop in 1:n.mix) {
      #ysim is the group mean
      #new:
      ysim = ysims[iloop, ]
      #new:
      #if (family$family == "gaussian"){
      ansi = coneA(ysim, amat, w=w, msg=FALSE)
      etahati = round(ansi$thetahat, 10)
      # } else {
      #   # umat = chol(v1)
      #   # uinv = solve(umat)
      #   # zsim = uinv %*% ysim
      #   # atil = amat %*% umat
      # 
      #   eans = eigen(v1, symmetric=TRUE)
      #   evecs = eans$vectors; evecs = round(evecs, 8)
      #   evals = eans$values
      #   sm = 1e-4
      #   if(any(evals < sm)){
      #     evals[which(evals < sm)] = sm
      #   }
      #   umat = evecs %*% diag(sqrt(evals))
      #   uinv = solve(umat) #solve will give singular error message when n is small
      #   atil = amat %*% umat
      #   zsim = uinv %*% ysim
      # 
      #   ansi = coneA(zsim, atil, msg=FALSE)
      #   phihati = round(ansi$thetahat, 10)
      #   etahati = umat %*% phihati
      #   #etahati = L %*% phihati
      # }
      
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
      imat = diag(M)
      sm_id = which((df.face[,2]/n.mix) < 1e-3)
      if (any(sm_id)) {
        df.face = df.face[-sm_id, ,drop=FALSE]
        sector = sector[-sm_id, ,drop=FALSE]
      }
      nsec = nrow(df.face)
      bsec = df.face
      bsec[,2] = bsec[,2] / sum(bsec[,2])
      
      acov = matrix(0, nrow=M, ncol=M)
      #check more:
      #if(family$family == 'gaussian') {
      wtinv = diag(1/w)
      # } else {
      #   #umat = chol(v1)
      #   #wtinv = chol2inv(umat)
      # 
      #   eans = eigen(v1, symmetric=TRUE)
      #   evecs = eans$vectors; evecs = round(evecs, 8)
      #   evals = eans$values
      #   sm = 1e-4
      #   if(any(evals < sm)){
      #     evals[which(evals < sm)] = sm
      #   }
      #   umat = evecs %*% diag(sqrt(evals))
      #   uinv = solve(umat)
      #   wtinv = uinv %*% t(uinv)
      # }
      
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
    lwr = etahat - hl
    upp = etahat + hl
  }  else {
    z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
    #z.mult = 2
    vsc = diag(v1)
    hl = z.mult*sqrt(vsc)
    lwr = etahat - hl
    upp = etahat + hl
  }
  
  vscu = diag(v1)
  hlu = z.mult*sqrt(vscu)
  lwru = etahatu - hlu
  uppu = etahatu + hlu
  
  if (ne >= 1) {
    #etahatu = yvecu
    #new: if there's empty cells, just augment muhatu to include NA's
    etahatu_all = 1:(M+zeros)*0
    etahatu_all[zeros_ps] = NA
    etahatu_all[obs_cells] = etahatu
    etahatu = etahatu_all
    
    lwru_all = 1:(M+zeros)*0
    lwru_all[zeros_ps] = NA
    lwru_all[obs_cells] = lwru
    lwru = lwru_all
    
    uppu_all = 1:(M+zeros)*0
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
    #revise later: replace lwru with lwr etc
    ans_im = impute_em2(empty_cells, obs_cells, M=(M+zeros), etahatu, Nd, nd, Nhat, w, 
                        domain.ord, grid2, Ds, sh, etahat, lwr=lwr, upp=upp, vsc=vsc, amat_0)
    etahat = ans_im$muhat
    lwr = ans_im$lwr
    upp = ans_im$upp
    vsc = ans_im$vsc
    domain.ord = ans_im$domain.ord
  }
  
  sig1 = vsc
  sig2 = sig1
  vsc_mix = sig1
  
  hl = z.mult*sqrt(vsc_mix)
  lwr = etahat - hl
  upp = etahat + hl
  
  #new: learn from cgam; nasa's idea
  #for now, only handles one predictor
  #new: don't constrain endpoints to follow monotonicity; make a copy of old lwr and upp
  lwr_0 = lwr
  upp_0 = upp
  nci = length(lwr)
  
  #new: check monotonicity of lwr and upp
  sh_0 = sh
  #relabel sh == 1 and 9?
  if(any(sh == 9)){
    sh_0[which(sh == 9)] = 1
  }
  #remove sh = 0, there is no amat for sh = 0
  Ds_0 = Ds
  if(any(sh == 0)){
    rm_id = which(sh == 0)
    sh_0 = sh_0[-rm_id]
    Ds_0 = Ds[-rm_id]
  }
  
  #print (dim(amat_0))
  #print (length(lwr))
  
  check_ps_lwr = round(c(amat_0 %*% lwr), 4) 
  check_ps_upp = round(c(amat_0 %*% upp), 4) 
  
  not_monotone_ci = any(check_ps_lwr < 0) | any(check_ps_upp < 0)
  
  #check more
  wvec = nd/sum(nd)
  if(any(wvec == 0)){
    wvec[which(wvec == 0)] = 1e-4
  }
  
  vsc_mix_0 = vsc_mix
  if (length(sh_0) == 1 & n.mix > 0 & not_monotone_ci) {
    #lwr = fix_monotone_1D(lwr, sh_0, nd)
    #upp = fix_monotone_1D(upp, sh_0, nd)
    #new:
    #wvec1 = fix_monotone_coneA_1D(lwr, sh_0, nd)
    #wvec2 = fix_monotone_coneA_1D(upp, sh_0, nd)
    #print (wvec2)
    lwr = coneA(lwr, amat_0, w = wvec)$thetahat
    upp = coneA(upp, amat_0, w = wvec)$thetahat
    vsc_mix_0 = ((upp - lwr) / 4)^2
  }
  
  #new:
  check_ps_etahat = round(c(amat_0 %*% etahat), 4) 
  if (length(sh_0) > 1 & n.mix > 0 & any(check_ps_etahat < 0) ) {
    #wvec = fix_monotone_coneA_2D(amat_ci, amat = amat_0, ci = etahat, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    etahat = coneA(etahat, amat_0, w = wvec)$thetahat
  }
  
  if (length(sh_0) > 1 & n.mix > 0 & not_monotone_ci) {
    #etahat = fix_monotone_2D(amat_ci, amat = amat_0, ci = etahat, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    #lwr = fix_monotone_2D(amat_ci, amat = amat_0, ci = lwr, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    #upp = fix_monotone_2D(amat_ci, amat = amat_0, ci = upp, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    #wvec1 = fix_monotone_coneA_2D(amat_ci, amat = amat_0, ci = lwr, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    #wvec2 = fix_monotone_coneA_2D(amat_ci, amat = amat_0, ci = upp, grid = grid2, Ds = Ds_0, nd, sh = sh_0, sh_all = sh, xvs2)
    
    lwr = coneA(lwr, amat_0, w = wvec)$thetahat
    upp = coneA(upp, amat_0, w = wvec)$thetahat
    #new:
    vsc_mix_0 = ((upp - lwr) / 4)^2
  }
  
  # vscu = diag(v1)
  # hlu = z.mult*sqrt(vscu)
  # lwru = etahatu - hlu
  # uppu = etahatu + hlu
  
  #new: test of H_0: theta in V and H_1: theta in C
  #do lsq fit with Atheta = 0 to get fit under H0
  #use the 1st option to compute T1 and T2
  bval = NULL
  pval = NULL
  if(test){
    m = qr(amat)$rank
    #vmat = qr.Q(qr(t(amat)), complete = TRUE)[, -(1:(qr(t(amat))$rank)), drop = FALSE]
    vec = amat %*% yvecu
    tst = amat %*% v1 %*% t(amat)
    #print (ginv(tst))
    T1_hat = t(vec) %*% MASS::ginv(amat %*% v1 %*% t(amat)) %*% vec
    
    eans = eigen(v1, symmetric=TRUE)
    evecs = eans$vectors; evecs = round(evecs, 8)
    #new: change negative eigenvalues to be a small positive value
    evals = eans$values
    sm = 1e-7
    
    #need test more
    #neg_ps = which(evals < sm)
    #if(length(neg_ps) > 0){
    #  evals[neg_ps] = sm
    #}
    L = evecs %*% diag(sqrt(evals))
    Linv = solve(L) #solve will give singular error message when n is small
    atil = amat %*% L
    Z_s = Linv %*% yvecu
    
    theta_hat = coneA(Z_s, atil, msg=FALSE)$thetahat
    T2_hat = t(Z_s-theta_hat) %*% (Z_s-theta_hat)
    
    bval = (T1_hat-T2_hat)/T1_hat
    
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
      pval = 1 - sum(pbeta(bval, m:0/2, 0:m/2)*mdist)
    } else {pval = 1}
  }
  #set h_grid = NULL for now
  h_grid = NULL
  muhat = c(family$linkinv(etahat))
  muhat.un = c(family$linkinv(etahatu))
  muhat_all = 1:n*0
  etahat_all = 1:n*0
  #new: for null deviance of gaussian case
  #ybar_all = 1:n*0
  for (udi in uds){
    udi = as.numeric(udi)
    ps = which(ID %in% udi)
    #should use muhat, not etahat
    muhat_all[ps] = muhat[which(uds%in%udi)]
    etahat_all[ps] = etahat[which(uds%in%udi)]
    #ybar_all[ps] = ybar_vec[which(uds%in%udi)]
  }
  #weights or prior.weights?
  #prior.weights seems correct
  deviance = sum(family$dev.resids(y, muhat_all, prior.weights))
  #CIC = t(evec) %*% wt %*% evec + 2*(sum(diag(mat)))
  #same as CIC:
  #for gaussian: -2loglike = nlog(t(evec) %*% wt %*% evec)
  #for binomial: yvec is not 0-1, and it doesn't work, then the penalty won't match
  #logLike_CIC = sum(family$dev.resids(yvec, muhat, diag(wt))) + 2*(sum(diag(mat)))
  logLike_CIC = deviance + 2*(sum(diag(mat)))
  ans = new.env()
  #revise later: about domain.ord
  ans$linear.predictors = etahat_all
  ans$etahat = c(etahat)
  ans$linear.predictors.un = ans.unc$linear.predictors
  ans$etahatu = c(yvecu)
  ans$fitted.values = muhat_all
  ans$muhat = muhat
  ans$fitted.values.un = ans.unc$fitted.values
  ans$muhatu = muhat.un
  ans$cov.un = v1 #unconstrained cov
  ans$cov.unscaled = vsc_mix #constrained cov with imputation
  if(n.mix == 0L){
    ans$acov = v1
  } else {
    ans$acov = acov #mixture variance-covariance without imputation
  }
  if(family$family == "gaussian"){
    ans$null.deviance = ans.unc_null$deviance
    #ans$null.deviance = sum(family$dev.resids(y, mean(y), prior.weights))
    #ans$null.deviance = sum(family$dev.resids(y, ybar_all, prior.weights))
  }else{
    ans$null.deviance = ans.unc$null.deviance
  }
  #ans$df.null = n - 1
  ans$df.null = M - 1
  ans$deviance = deviance
  #ans$df.residual = n - edf
  ans$df.residual = M - edf
  ans$edf = edf 
  ans$lwr = c(lwr)
  ans$upp = c(upp)
  ans$lwru = c(lwru)
  ans$uppu = c(uppu)
  ans$y = y
  ID = as.numeric(levels(ID))[ID] 
  ans$domain = ID
  ans$grid_ps = pskeep
  ans$amat = amat
  ans$Ds = Ds
  ans$domain.ord = domain.ord
  ans$nd = nd
  ans$grid = grid2
  ans$ec = empty_cells
  ans$hkeep = hkeep
  ans$h_grid = h_grid
  ans$zeros_ps = zeros_ps
  ans$pval = pval
  ans$bval = bval
  ans$CIC = CIC
  ans$CIC.un = CIC.un
  ans$logLike_CIC = logLike_CIC
  #new: inherit attributes if svyby is used
  if (family$family == "gaussian") {
    ans.unc_cp = data.frame(matrix(nrow = M_0, ncol = ncol(ans.unc)))
    colnames(ans.unc_cp) = colnames(ans.unc)
    #ans.unc_cp$ID = factor(1:M_0)
    #ans.unc_cp = ans.unc_cp_0
    ans.unc_cp[[ynm]] = as.numeric(etahat)
    #ans.unc_cp$se = vsc_mix^(.5)
    ans.unc_cp$se = vsc_mix_0^(.5) #new
    attr(ans.unc_cp, "class") = c("svyby", "data.frame")
    attr(ans.unc_cp, "svyby") = attr(ans.unc, "svyby")
    attr(ans.unc_cp, "var") = acov
  } else {
    ans.unc_cp = data.frame(matrix(nrow = M_0, ncol = ncol(ans.unc_cp_0)))
    colnames(ans.unc_cp) = colnames(ans.unc_cp_0)
    ans.unc_cp = ans.unc_cp_0
    ans.unc_cp[[ynm]] = as.numeric(etahat)
    #ans.unc_cp$se = vsc_mix^(.5)
    ans.unc_cp$se = vsc_mix_0^(.5) #new
    attr(ans.unc_cp, "class") = c("svyby", "data.frame")
    attr(ans.unc_cp, "svyby") = attr(ans.unc_cp_0, "svyby")
    attr(ans.unc_cp, "var") = acov
  }
  ans$w = w
  ans$xm.red = xm.red
  ans$ne = ne
  ans$empty_cells = empty_cells
  ans$small_cells = small_cells
  ans$obs_cells = obs_cells
  ans$zeros = zeros
  ans$Nd = Nd
  ans$Nhat = Nhat
  ans$amat_0 = amat_0
  ans$amat_ci = amat_ci
  #very little difference between weights and prior.weights?
  ans$weights = weights
  ans$prior.weights = prior.weights
  ans$ans.unc = ans.unc
  ans$ans.unc_cp = ans.unc_cp
  ans$offset = ans.unc$offset
  ans$contrasts = ans.unc$contrasts
  return (ans)
}

####################################
#inherit from cgam
#apply plotpersp to a csvy.fit object
#####################################
# plotpersp = function(object,...) {
# #x1nm = deparse(substitute(x1))
# #x2nm = deparse(substitute(x2))
#  UseMethod("plotpersp", object)
# }

#--------------------------------------------------
#ci can be ci or upp
#need to include block ordering
#--------------------------------------------------
fix_monotone_1D = function(ci, sh, nd,...){
  check_ps_ci = round(diff(ci), 6)
  bool_ci1 = any(check_ps_ci < 0) & sh == 1
  bool_ci2 = any(check_ps_ci > 0) & sh == 2
  not_monotone_ci = ifelse(bool_ci1 | bool_ci2, TRUE, FALSE)
  nci = length(ci)
  
  nrep = 0
  while(not_monotone_ci & nrep < 20) {
    nrep = nrep + 1
    if (sh == 1) {
      check_id_ci = which(check_ps_ci < 0)
      n_check_ci = length(check_id_ci)
      for(k in 1:n_check_ci){
        i = check_id_ci[k]
        if (any(which(ci[(i + 1):nci] < ci[i]))){
          ps_left = i
          check_ps = which(ci[(i + 1):nci] < ci[i])
          ps_right = max(((i+1):nci)[check_ps])
          nd_lr = nd[ps_left:ps_right]
          if (any(nd_lr > 0)){
            ws = nd_lr / sum(nd_lr)
            ci[ps_left:ps_right] = sum(ws * ci[ps_left:ps_right])
          }   
        }
      }
    }
    
    if (sh == 2){
      check_id_ci = which(check_ps_ci > 0)
      n_check_ci = length(check_id_ci)
      for(k in n_check_ci:1){
        #the domain that violates the constraint is the ith domain
        i = check_id_ci[k]
        if (any(which(ci[1:i] < ci[i + 1]))){
          check_ps = which(ci[1:i] < ci[i + 1])
          ps_left = min((1:i)[check_ps])
          ps_right = i + 1
          nd_lr = nd[ps_left:ps_right]
          if (any(nd_lr > 0)){
            ws = nd_lr / sum(nd_lr)
            ci[ps_left:ps_right] = sum(ws * ci[ps_left:ps_right])
          }       
        }
      }
    }
    
    check_ps_ci = round(diff(ci), 6)
    bool_ci1 = any(check_ps_ci < 0) & sh == 1
    bool_ci2 = any(check_ps_ci > 0) & sh == 2
    not_monotone_ci = ifelse(bool_ci1 | bool_ci2, TRUE, FALSE)
  }
  return (ci)
}

#--------------------------------------------------
fix_monotone_2D = function(amat_ci, amat, ci, grid, Ds, nd, sh, sh_all, xvs2,...){
  if(all(sh == 1)) {
    ci = fix_monotone_2D_0(amat, ci, grid, Ds, nd)$ci
  }
  if(all(sh == 2)){
    ci = fix_monotone_2D_0(-amat, -ci, grid, Ds, nd)$ci
    ci = -ci
  }
  if(any(sh == 1) & any(sh == 2)){
    nD = length(Ds)
    decr_ps = which(sh_all == 2) #sh_all include sh = 0, xvs2 include sh = 0
    
    xvs2_0 = xvs2
    for(i in decr_ps){
      if(is.matrix(xvs2_0)){
        xvs2_0[,i] = rev(xvs2_0[,i])
      }else if(is.list(xvs2_0)){
        xvs2_0[[i]] = rev(xvs2[[i]])
      }
    }
    
    if(is.matrix(xvs2_0)){
      grid_rev = expand.grid(as.data.frame(xvs2_0))
    }else if(is.list(xvs2_0)){
      grid_rev = expand.grid(xvs2_0)
    }
    
    ps_rev = apply(grid_rev, 1, function(elem) which(apply(grid, 1, function(gi) all(gi == elem))))
    ci_rev = ci[ps_rev]
    
    amat_rev = amat
    sh_track = NULL
    for(i in 1:nD){
      amati = amat_ci[[i]]
      shi = sh[i]
      sh_track = c(sh_track, rep(shi, nrow(amati)))
    }
    amat_rev[which(sh_track == 2), ] = -amat[which(sh_track == 2), , drop = FALSE]
    
    nd_rev = nd[ps_rev]
    Ds_rev = Ds #?
    ci_rev = fix_monotone_2D_0(amat_rev, ci_rev, grid_rev, Ds, nd_rev)$ci
    ci = ci_rev[ps_rev]
  }
  return (ci)
}

#--------------------------------------------------
fix_monotone_2D_0 = function(amat, ci, grid, Ds, nd,...){
  sm = 1e-5
  check_ps_ci = round(c(amat %*% ci), 5) #4th or 5th digit?
  not_monotone_ci = any(check_ps_ci < 0)
  nrep = 0
  upp_limit = 1000 # should be 2^D1 * 2^D2?
  clusters_merged = list()
  
  while(not_monotone_ci & nrep < upp_limit){
    nrep = nrep + 1
    check_id_ci = which(check_ps_ci < 0)
    amat_not_mono = amat[check_id_ci, ,drop = FALSE]
    grid_ci = cbind(grid, ci, t(amat_not_mono)) #no use, for visual check only
    #different from 1D case
    ps_left_cands = which(apply(amat_not_mono, 2, function(e) any(-1 %in% e)))
    ps_left_cands = sort(ps_left_cands)
    n_lefts = length(ps_left_cands)
    
    rows_lefts = NULL
    clusters = vector('list', length = n_lefts) #each cluster will be weight-averaged
    for(k in 1:n_lefts){
      ps_left_k = ps_left_cands[k]
      colk = amat_not_mono[, ps_left_k, drop = TRUE]
      row_ps = which(colk == -1)
      rows_lefts = c(rows_lefts, row_ps)
      rowk = amat_not_mono[row_ps, ,drop = FALSE] #need to move along rowk to find 1
      ps_right_k = apply(rowk, 1, function(es) max(which(es == 1)))
      clusters[[k]] = c(ps_left_k, ps_right_k)
    }
    
    #make a copy
    clusters_0 = clusters
    #merge elements in clusters which have common domains
    #clusters_merged = list()
    iter = 1
    nc = length(clusters)
    
    #----------------------------------------------------------------------------------------------------
    #found here: https://stackoverflow.com/questions/47322126/merging-list-with-common-elements
    #----------------------------------------------------------------------------------------------------
    i = rep(1:length(clusters_0), lengths(clusters_0))
    j = factor(unlist(clusters_0))
    tab = sparseMatrix(i = i, j = as.integer(j), x = TRUE, dimnames = list(NULL, levels(j)))
    connects = Matrix::tcrossprod(tab, boolArith = TRUE)
    group = clusters(graph_from_adjacency_matrix(connects, mode = "undirected"))$membership
    clusters_merged = tapply(clusters_0, group, function(x) sort(unique(unlist(x))))
    
    nc_merged = length(clusters_merged)
    for(k in 1:nc_merged){
      ps_lr = clusters_merged[[k]]
      #print (ps_lr)
      nd_lr = nd[ps_lr]
      #print (nd_lr)
      if (any(nd_lr > 0)){
        ws = nd_lr / sum(nd_lr)
        ci[ps_lr] = sum(ws * ci[ps_lr])
        #temp
        if(any(nd_lr == 0)){
          nd_lr[which(nd_lr == 0)] = round(mean(nd_lr[which(nd_lr != 0)]))
          nd[ps_lr] = nd_lr
        }
      }   
      
      #temp: difficult case, the merged cluster consists of zero cells only
      if(all(nd_lr == 0)){
        ci_left_most = ci[ps_lr[1]]
        ci[ps_lr] = ci_left_most
      }
    }
    check_ps_ci = round(c(amat %*% ci), 5)
    not_monotone_ci = any(check_ps_ci < 0)
  }
  #print (nrep)
  #rslt = list(ci = ci, nrep = nrep, clusters_merged = clusters_merged)
  rslt = list(ci = ci, nrep = nrep)
  return (rslt)
}



#########################################
#new plotpersp
#########################################
plotpersp.csvy = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL, data=NULL, ci=c("none", "lwr", "up", "both"),
                          transpose=FALSE, main=NULL, categ=NULL, categnm=NULL, 
                          surface = c("C","U"), type = c("link", "response"),
                          col = "white", cex.main=.8, xlab=NULL, 
                          ylab=NULL, zlab=NULL, zlim=NULL, box=TRUE,
                          axes=TRUE, th=NULL, ltheta=NULL, 
                          ticktype="detailed", nticks=5, palette=NULL, NCOL=NULL,...) {
  x1nm = deparse(substitute(x1))
  x2nm = deparse(substitute(x2))
  type = match.arg(type)
  surface = match.arg(surface)
  ci = match.arg(ci)
  #print (ci)
  if (!inherits(object, "csvy")) {
    warning("calling plotpersp(<fake-csvy-object>) ...")
  }
  xnms = object$xnms_add
  shp = object$shapes_add
  Ds = object$Ds
  Nd = length(Ds)
  xmat = object$xmat_add
  bool = apply(xmat, 2, function(x) is.numeric(x))
  if(any(!bool)){
    xmat = apply(xmat, 2, function(x) as.numeric(x))
  }
  ynm = object$ynm
  family = object$family
  etahat = object$etahat
  lwr = object$lwr
  upp = object$upp
  grid = object$grid
  main = ifelse(is.null(main), "Constrained Fit", main)
  
  switch(surface, U = {
    etahat = object$etahatu
    lwr = object$lwru
    upp = object$uppu
    main = "Unconstrained Fit"
  }, C = )
  
  switch(type, response = {
    etahat = family$linkinv(etahat)
    lwr = family$linkinv(lwr)
    upp = family$linkinv(upp)
  }, link = )
  
  obs = 1:Nd
  if (x1nm == "NULL" | x2nm == "NULL") {
    if (length(xnms) >= 2) {
      if (is.null(categ)) {
        x1nm = xnms[1]
        x2nm = xnms[2]
        x1_id = 1
        x2_id = 2
      } else {
        if (!is.character(categ)) {
          stop("categ must be a character argument!")
        } else if (any(grepl(categ, xnms))) {
          id = which(grepl(categ, xnms))
          x12nms = xnms[-id]
          x12ids = obs[-id]
          x1nm = x12nms[1]
          x2nm = x12nms[2]
          x1_id = x12ids[1]
          x2_id = x12ids[2]
        }  
      }
    } else {stop ("Number of non-parametric predictors must >= 2!")}
  } else {
    #new: use the data as the environment for x1 and x2
    x1 = data[[x1nm]]; x2 = data[[x2nm]]
    if (all(xnms != x1nm)) {
      if (length(x1) != nrow(xmat)) {
        stop(cat("Number of observations in the data set is not the same as the number of elements in x1!"), "\n")
      }
      bool = apply(xmat, 2, function(x) all(x1 == x))
      if (any(bool)) {
        x1_id = obs[bool]
        #change x1nm to be the one in formula
        x1nm = xnms[bool]
      } else {
        stop (cat(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the model!"), "\n")
      }
    } else {
      x1_id = obs[xnms == x1nm]
    }
    
    if (all(xnms != x2nm)) {
      if (length(x2) != nrow(xmat)) {
        stop (cat("Number of observations in the data set is not the same as the number of elements in x2!"), "\n")
      }
      bool = apply(xmat, 2, function(x) all(x2 == x))
      if (any(bool)) {
        x2_id = obs[bool]
        x2nm = xnms[bool]
      } else {
        stop (cat(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the model!", "\n"))
      }
    } else {
      x2_id = obs[xnms == x2nm]
    }
  }
  #xmat keeps the order of x's in the original data frame, but x1 and x2 may change depending on x1_id and x2_id
  
  #x12_id = c(2, 3)
  x12_id = c(x1_id, x2_id)
  # x1_id = x12_id[1]
  # x2_id = x12_id[2]
  if (Nd > 2) {
    x3_ids = obs[-x12_id] 
  }
  shp_12 = shp[x12_id]
  
  x1nm = xnms[x1_id]
  x2nm = xnms[x2_id]
  if (Nd > 2) {
    x3nms = xnms[x3_ids]
  }
  
  x1 = xmat[, x1_id]
  x2 = xmat[, x2_id]
  if (Nd > 2) {
    x3s = xmat[, x3_ids, drop = FALSE] #need change
    x3_vals = apply(x3s, 2, median)
  }
  
  if (Nd > 2) {
    if(is.null(categ)) {
      ps = 1:nrow(grid) #track rows with the median values of columns that are not x1 or x2
      nx3 = length(x3_vals)
      for(k in 1:nx3){
        x3_val_k = x3_vals[k]
        ps_k = which(grid[, x3_ids[k]] %in% x3_val_k)
        ps = intersect(ps, ps_k)
      }
      surf.etahat = matrix(etahat[ps], Ds[x1_id], Ds[x2_id])
      surf.lwr = matrix(lwr[ps], Ds[x1_id], Ds[x2_id])
      surf.upp = matrix(upp[ps], Ds[x1_id], Ds[x2_id])
    } else {
      if (!is.character(categ)) {
        stop("categ must be a character argument!")
      } else if (any(grepl(categ, x3nms))) {
        x3_id = which(xnms %in% categ)
        x3nm = xnms[x3_id]
        if (x3nm %in% c(x1nm, x2nm)) {
          stop("categ must be different than x1 and x2!")
        }
        x3p = 1:Ds[x3_id]
        nsurf = length(x3p)
        surf.etahat = vector('list', nsurf)
        surf.lwr = vector('list', nsurf)
        surf.upp = vector('list', nsurf)
        
        #test more
        ps = 1:nrow(grid) #track rows with the median values of columns that are not x1 or x2
        x3_others_id = which(!x3nms %in% x3nm)
        x3_vals_others = x3_vals[x3_others_id]
        nx3 = length(x3_others_id)
        if(nx3 >= 1){#used when there are more than 3 x's
          for(k in 1:nx3){
            x3_other_val_k = x3_vals_others[k]
            colk = x3_ids[x3_others_id[k]]
            ps_k = which(grid[, colk] %in% x3_other_val_k)
            ps = intersect(ps, ps_k)
          }
        }
        
        for(k in 1:nsurf){
          x3_val_k = x3p[k]
          ps_k = intersect(ps, which(grid[, x3_id] %in% x3_val_k))
          surf.etahat[[k]] = matrix(etahat[ps_k], Ds[x1_id], Ds[x2_id])
          surf.lwr[[k]] = matrix(lwr[ps_k], Ds[x1_id], Ds[x2_id])
          surf.upp[[k]] = matrix(upp[ps_k], Ds[x1_id], Ds[x2_id])
        }
      } else {
        stop(cat(categ, "is not an exact character name defined in the csurvey fit!", "\n"))
      }
      
      if(is.null(palette)){
        palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", 
                    "seagreen1", "palegreen", "azure", "whitesmoke")
      }
      if (nsurf > 1 && nsurf < 11) {
        col = palette[1:nsurf]
      } else {
        #new: use rainbow
        col = topo.colors(nsurf)
      }
      width = nsurf^.5
      wd = round(width)
      if (width > wd) {
        wd = wd + 1
      }  
      if ((wd^2 - nsurf) >= wd ) {
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
    }
  } else if (Nd == 2) {
    surf.etahat = matrix(etahat, Ds[x1_id], Ds[x2_id])
    surf.lwr = matrix(lwr, Ds[x1_id], Ds[x2_id])
    surf.upp = matrix(upp, Ds[x1_id], Ds[x2_id])
  } else if (Nd == 1) {
    stop('persp plot only works when there are >= 2 predictors!')
  }
  
  t_col = function(color, percent = 80, name = NULL) {
    rgb.val = col2rgb(color)
    ## Make new color using input color as base and alpha set by trobjectparency
    t.col = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                maxColorValue = 255,
                alpha = (100-percent)*255/100,names = name)
    ## Save the color
    invisible(t.col)
  }
  col.upp = t_col("green", perc = 90, name = "lt.green")
  col.lwr = t_col("pink", perc = 80, name = "lt.pink")
  
  if (is.null(xlab)) {
    #xlab = deparse(x1nm)
    xlab = x1nm
  }
  if (is.null(ylab)) {
    #ylab = deparse(x2nm)
    ylab = x2nm
  }
  if (is.null(zlab)) {
    if("response" %in% type){
      zlab = paste("Est mean of", ynm)
    } 
    if("link" %in% type){
      zlab = paste("Est systematic component of", ynm)
    }
  }
  
  ang = NULL
  if (is.null(th)) {
    if (all(shp_12 == 1)) {
      ang = -40
    } else if (all(shp_12 == 2)){
      ang = 40
    } else if (shp_12[1] == 1 & shp_12[2] == 2) {
      ang = -50
    } else if (shp_12[1] == 2 & shp_12[2] == 1) {
      ang = -230
    } else {
      ang = -40
    }
  } else {
    ang = th
  }
  
  x1p = 1:Ds[x1_id]
  x2p = 1:Ds[x2_id]
  
  if (is.list(surf.upp)) {
    zlim0 = vector('list', length = nsurf)
    for(i in 1:nsurf){
      if (ci != "none") {
        rg = sapply(surf.upp, max) - sapply(surf.lwr, min)
        zlim0[[i]] = c(min(surf.lwr[[i]]) - rg[i]/8, max(surf.upp[[i]]) + rg[i]/8)
      } else {
        rg = sapply(surf.etahat, max) - sapply(surf.etahat, min)
        zlim0[[i]] = c(min(surf.etahat[[i]])  - rg[i]/10, max(surf.etahat[[i]]) + rg[i]/10)
      }     
    }
    #zlim0 = c(min(surf.lwr[[1]]) - rg/5, max(surf.upp[[1]]) + rg/5)
    #print (zlim0)
  } else {
    rg = max(surf.upp) - min(surf.lwr)
    zlim0 = c(min(surf.lwr) - rg/5, max(surf.upp) + rg/5)
  }
  
  if (is.list(surf.etahat)) {#this means categ is not null
    if(is.null(NCOL)){
      old.par = par(mfrow = c(fm[1], fm[2]))
    } else {
      nr = fm[1]*fm[2]/NCOL
      old.par = par(mfrow = c(nr, NCOL))
    }
    on.exit(par(old.par))
    par(mar = c(4, 1, 1, 1))
    par(cex.main = cex.main)
    for(i in 1:nsurf){
      if (is.null(categnm)) {
        persp(x1p, x2p, surf.etahat[[i]], zlim = zlim0[[i]], col=col[i], xlab = xlab,
              ylab = ylab, zlab = zlab, main = paste0(x3nm, " = ", x3p[i]), 
              theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes, nticks=nticks)
        if (ci == "lwr" | ci == "both") {
          par(new = TRUE)
          persp(x1p, x2p, surf.lwr[[i]], zlim = zlim0[[i]], col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
        }
        if (ci == "up" | ci == "both") {
          par(new = TRUE)
          persp(x1p, x2p, surf.upp[[i]], zlim = zlim0[[i]], col=col.upp, theta = ang, box=FALSE, axes=FALSE)
        }
      } else {
        if (length(categnm) == nsurf) {
          persp(x1p, x2p, surf.etahat[[i]], zlim = zlim0[[i]], col=col[i], xlab = xlab,
                ylab = ylab, zlab = zlab, main = categnm[i], 
                theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes, nticks=nticks)
          if (ci == "lwr" | ci == "both") {
            par(new = TRUE)
            persp(x1p, x2p, surf.lwr[[i]], zlim = zlim0[[i]], col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
          }
          if (ci == "up" | ci == "both") {
            par(new = TRUE)
            persp(x1p, x2p, surf.upp[[i]], zlim = zlim0[[i]], col=col.upp, theta = ang, box=FALSE, axes=FALSE)
          }
        } else if (length(categnm) == 1) {
          persp(x1p, x2p, surf.etahat[[i]], zlim = zlim0[[i]], col=col[i], xlab = xlab,
                ylab = ylab, zlab = zlab, main = paste0(categnm, " = ", x3p[i]), 
                theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes, nticks=nticks)
          if (ci == "lwr" | ci == "both") {
            par(new = TRUE)
            persp(x1p, x2p, surf.lwr[[i]], zlim = zlim0[[i]], col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
          }
          if (ci == "up" | ci == "both") {
            par(new = TRUE)
            persp(x1p, x2p, surf.upp[[i]], zlim = zlim0[[i]], col=col.upp, theta = ang, box=FALSE, axes=FALSE)
          }
        }
      }
    }
  } else {
    #old.par = ifelse(ci == "both", par(mfrow=c(1, 2)),  par(mfrow=c(1, 1)))
    if (ci == "both") {old.par = par(mfrow=c(1, 2))} else {old.par = par(mfrow=c(1, 1))}
    persp(x1p, x2p, surf.etahat, zlim = zlim0, 
          col=col, xlab = xlab, ylab = ylab, zlab = zlab, main = main,
          theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes, nticks=nticks)
    if (ci == "up") {
      par(new = TRUE)
      persp(x1p, x2p, surf.upp, zlim = zlim0, col=col.upp, theta = ang, box=FALSE, axes=FALSE)
    }
    #par(new = FALSE)
    
    if (ci == "lwr") {
      par(new = TRUE)
      persp(x1p, x2p, surf.lwr, zlim = zlim0, col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
    }
    
    if (ci == "both") {
      par(new = TRUE)
      persp(x1p, x2p, surf.upp, zlim = zlim0, col=col.upp, theta = ang, box=FALSE, axes=FALSE)
      par(new=FALSE)
      
      persp(x1p, x2p, surf.etahat, zlim = zlim0, 
            col=col, xlab = xlab, ylab = ylab, zlab = zlab, main = main,
            theta = ang, ltheta = -135, ticktype = ticktype, box=box, axes=axes, nticks=nticks)
      par(new = TRUE)
      persp(x1p, x2p, surf.lwr, zlim = zlim0, col=col.lwr, theta = ang, box=FALSE, axes=FALSE)
    }
    
    par(new=FALSE)
    on.exit(par(old.par))
  }
  #rslt = list(surf.etahat = surf.etahat, surf.upp = surf.upp, surf.lwr = surf.lwr, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, theta = ang, 
  #            ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype)
  #invisible(rslt)
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
# constr = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 0
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# incr = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 1
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

#relaxed increasing
# relax.incr = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 1
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   attr(x, "relax") = TRUE
#   #class(x) = "additive"
#   x
# }

# decr = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 2
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

#relaxed decreasing
# relax.decr = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 2
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   attr(x, "relax") = TRUE
#   #class(x) = "additive"
#   x
# }

# conv = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 3
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# conc = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 4
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# incr.conv = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 5
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# decr.conv = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 6
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# incr.conc = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 7
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

# decr.conc = function(x, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 8
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "additive"
#   #class(x) = "additive"
#   x
# }

block.Ord = function(x, order = NULL, numknots = 0, knots = 0, space = "E")
{
  cl = match.call()
  pars = match.call()[-1]
  attr(x, "nm") = deparse(pars$x)
  attr(x, "shape") = 9
  attr(x, "numknots") = numknots
  attr(x, "knots") = knots
  attr(x, "space") = space
  attr(x, "categ") = "block.ord"
  attr(x, "order") = order
  #class(x) = "additive"
  x
}

# block.Ave = function(x, order = NULL, numknots = 0, knots = 0, space = "E")
# {
#   cl = match.call()
#   pars = match.call()[-1]
#   attr(x, "nm") = deparse(pars$x)
#   attr(x, "shape") = 10
#   attr(x, "numknots") = numknots
#   attr(x, "knots") = knots
#   attr(x, "space") = space
#   attr(x, "categ") = "block.ave"
#   attr(x, "order") = order
#   #class(x) = "additive"
#   x
# }

######
#1D
######
makeamat = function(x, sh, Ds = NULL, suppre = FALSE, interp = FALSE, 
                    relaxes = FALSE, h = NULL, block.ave = NULL, 
                    block.ord = NULL, zeros_ps = NULL, noord_cells = NULL,
                    D_index = 1, xvs2 = NULL, grid2 = NULL) {
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
    #print (amat)
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
            c1 = min(obs[abs(x - xu[i]) < sm]); c2 = min(obs[abs(x - xu[i + 1]) < sm])
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
      
      #if(ncol(grid2) == 1){
      rm_id = unique(sort(c(zeros_ps, noord_ps)))
      #} else {
      #  rm_id = unique(sort(noord_ps))
      #}
      if(length(rm_id)>0){
        block.ord_nz = block.ord[-rm_id]
      }
      
      ubck = unique(block.ord_nz)
      #use values in block.ord as an ordered integer vector
      ubck = sort(ubck)
      
      nbcks = length(table(block.ord_nz))
      szs = as.vector(table(block.ord_nz))
      
      #new:
      if(is.list(xvs2)){
        ps = sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
      }else if(is.matrix(xvs2)){
        #D_index=1
        #ps = lapply(xvs2, function(.x) which(as.vector(grid2) %in% .x)[1])
        #ps = as.matrix(grid2)
        #test! wrong to use x
        ps = 1:M
      }
      
      #if(length(noord_ps) > 0){
      #  ps = ps[-noord_ps]
      #}
      if(length(rm_id) > 0){
        ps = ps[-rm_id]
      }
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
        
        #ps_bck1 = ps[which(block.ord_nz == ubck1)]
        #ps_bck2 = ps[which(block.ord_nz == ubck2)]
        #nr = l1*l2; nc = M
        #for each element in block1, the value for each element in block 2 will be the same
        #for each element in block1, the number of comparisons is l2
        amat_bcki = NULL
        rm_l = length(zeros_ps) + length(noord_cells)
        #rm_l = length(zeros_ps) 
        #rm_l = 0
        #test more
        #l_noord_ps = length(which(grid2[, D_index] %in% noord_ps))
        #rm_l = length(zeros_ps) + l_noord_ps
        
        #tmp:
        #if(D_index == 1){
        #  M.ord = length(block.ord)
        #}else{
        #M.ord = nc2
        #  M.ord = nrow(grid2)
        #}
        M.ord = length(block.ord)
        
        #amatj = matrix(0, nrow = l2, ncol = M.ord)
        amatj = matrix(0, nrow = l2, ncol = M.ord-rm_l)
        bck_1k = bck_1[1]
        #bck_1k = ps_bck1[1]
        amatj[, bck_1k] = -1
        row_pointer = 1
        for(j in 1:l2){
          bck_2j = bck_2[j]
          #bck_2j = ps_bck2[j]
          amatj[row_pointer, bck_2j] = 1
          row_pointer = row_pointer + 1
        }
        amatj0 = amatj
        amat_bcki = rbind(amat_bcki, amatj)
        
        if(l1 >= 2){
          for(k in 2:l1) {
            bck_1k = bck_1[k]; bck_1k0 = bck_1[k-1]
            #bck_1k = ps_bck1[k]; bck_1k0 = ps_bck1[k-1]
            amatj = amatj0
            #set the value for the kth element in block 1 to be -1
            #set the value for the (k-1)th element in block 1 back to 0
            #keep the values for block 2
            amatj[, bck_1k] = -1; amatj[, bck_1k0] = 0
            amatj0 = amatj
            amat_bcki = rbind(amat_bcki, amatj)
          }
        }
        amat = rbind(amat, amat_bcki)
      }
      
      if(length(noord_cells)>0){
        nr = nrow(amat); nc = ncol(amat)
        amat_tmp = matrix(0, nrow = nr, ncol = (M.ord-rm_l+length(noord_cells)))
        if(length(zeros_ps) > 0){
          ps = which(block.ord[-zeros_ps] != 0)
        }else{
          ps = which(block.ord != 0)
        }
        amat_tmp[, ps] = amat
        amat = amat_tmp
      }
      return (amat)
    } else if ((length(block.ave) > 1) & all(block.ord == -1)) {
      block.ave_nz = block.ave
      #rm_id = unique(sort(c(zeros_ps, noord_cells)))
      
      noord_ps = which(block.ave == 0)
      rm_id = unique(sort(c(zeros_ps, noord_ps)))
      
      if(length(rm_id)>0){
        block.ave_nz = block.ave[-rm_id]
      }
      #if(length(zeros_ps)>0){
      #  block.ave_nz = block.ave[-zeros_ps]
      #}
      
      #if(length(noord_cells)>0){
      #  block.ave_nz = block.ave_nz[-noord_cells]
      #}
      
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
        #print (sz1)
        #print (sz2)
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
makeamat_2D = function(x=NULL, sh, grid2=NULL, xvs2=NULL, zeros_ps=NULL, noord_cells=NULL,
                       Ds = NULL, suppre = FALSE, interp = FALSE, 
                       relaxes = FALSE, block.ave.lst = NULL, block.ord.lst = NULL, M_0 = NULL) {
  #n = length(x)
  #xu = sort(unique(x))
  #n1 = length(xu)
  sm = 1e-7
  ms = NULL
  Nd = length(Ds)
  ne = length(zeros_ps)
  #new: for ci
  amat_ci = vector('list', length = Nd)
  if (Nd == 2) {
    if (any(sh > 0)) {
      D1 = Ds[1]
      D2 = Ds[2]
      obs = 1:(D1*D2)
      #new: group sh = 9 | 10 with other shapes
      if (sh[1] <= 10) {
        #new:
        if(sh[1] == 0){
          amat1 = NULL
        }else{
          if (length(zeros_ps) > 0) {
            amat0_lst = list()
            grid3 = grid2[-zeros_ps, ,drop=FALSE]
            row_id = as.numeric(rownames(grid3))
            cts = row_id[which(row_id%%D1 == 0)]
            #temp: check if zeros_ps contains some point %% D1 = 0
            if (any((zeros_ps %% D1) == 0)) {
              zeros_ps_at_cts = zeros_ps[which((zeros_ps %% D1) == 0)]
              cts_add = row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
              cts = unique(sort(c(cts, cts_add)))
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
                if(sh[1] < 9){
                  amat0i = makeamat(1:D1i, sh[1])
                } else {
                  #new:
                  zeros_ps_i = zerosi%%D1
                  if(any(zeros_ps_i %in% 0)){
                    zeros_ps_i[which(zeros_ps_i %in% 0)] = D1
                  } else {
                    zeros_ps_i = zerosi%%D1
                  }
                  #not sure....
                  amat0i = makeamat(1:D1i, sh[1], block.ave = block.ave.lst[[1]], 
                                    block.ord = block.ord.lst[[1]], 
                                    zeros_ps = zeros_ps_i, 
                                    noord_cells = noord_cells[[1]], 
                                    D_index = 1, xvs2=xvs2, grid2=grid2)
                }
                amat0_lst[[i]] = amat0i
              } else if (length(zerosi) >= 1 & (ed == st)) {
                #next
                amat0i = 0
                amat0_lst[[i]] = amat0i
              } else if (length(zerosi) == 0) {
                if(sh[1] < 9){
                  amat0i = makeamat(1:D1, sh[1])
                }else{
                  #new:
                  amat0i = makeamat(1:D1i, sh[1], block.ave = block.ave.lst[[1]], 
                                    block.ord = block.ord.lst[[1]], 
                                    zeros_ps = NULL, 
                                    noord_cells = noord_cells[[1]], 
                                    D_index = 1,
                                    xvs2=xvs2, grid2=grid2)
                }
                amat0_lst[[i]] = amat0i
              }
              #amat0_lst[[i]] = amat0i
            }
          } else {
            #if(sh[1] < 9){
            amat0 = makeamat(1:D1, sh=sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]],
                             zeros_ps=NULL, noord_cells = noord_cells[[1]], D_index = 1, xvs2 = xvs2, grid2 = grid2)
            #} else {
            #  amat0 = makeamat_2d_block(x=1:D1, sh=sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]], 
            #        zeros_ps=NULL, noord_cells = noord_cells[[1]], D_index = 1, xvs2 = xvs2, grid2 = grid2, M_0 = M_0)
            #}
            amat0_lst = rep(list(amat0), D2)
          }
          amat0_lst = Filter(Negate(is.null), amat0_lst)
          amat1 = bdiag(amat0_lst)
          amat1 = as.matrix(amat1)
          amat_ci[[1]] = amat1
        }
      }
      #2D
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
          amat2 = makeamat_2d_block(x=NULL, sh[2], block.ave = block.ave.lst[[2]], 
                                    block.ord = block.ord.lst[[2]], 
                                    zeros_ps = zeros_ps, 
                                    noord_cells = noord_cells[[2]], 
                                    D_index=2, xvs2=xvs2, grid2=grid2, M_0=M_0)
        }
      }
      amat = rbind(amat1, amat2)
      #new: for ci
      amat_ci[[2]] = amat2
    } 
  } else if (Nd >= 3) {
    M = length(Ds)
    D1 = Ds[1]
    cump = cumprod(Ds)
    nc = cump[M]
    amat1 = amat2 = amatm = NULL
    
    # 1D: group sh=9 or sh=10 with other shapes in 1D
    if(sh[1] == 0){
      amat1 = NULL
    } else {
      if(length(zeros_ps)>0){
        amat0_lst = list()
        grid3 = grid2[-zeros_ps, ,drop=FALSE]
        row_id = as.numeric(rownames(grid3))
        cts = row_id[which(row_id%%D1 == 0)]
        if(any((zeros_ps %% D1) == 0)){
          zeros_ps_at_cts = zeros_ps[which((zeros_ps %% D1) == 0)]
          cts_add = row_id[sapply(zeros_ps_at_cts, function(elem) max(which(row_id < elem)))]
          #cts = sort(c(cts, cts_add))
          cts = unique(sort(c(cts, cts_add)))
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
            #amat0i = makeamat(1:D1i, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            #new:
            zeros_ps_i = zerosi%%D1
            if(any(zeros_ps_i %in% 0)){
              zeros_ps_i[which(zeros_ps_i %in% 0)] = D1
            } #else {
            #zeros_ps_i = zerosi%%D1
            #}
            amat0i = makeamat(1:D1i, sh[1], block.ave = block.ave.lst[[1]], block.ord = block.ord.lst[[1]], 
                              zeros_ps = zeros_ps_i, noord_cells = noord_cells[[1]], xvs2=xvs2, grid2=grid2)
            amat0_lst[[i]] = amat0i
          } else if (length(zerosi) >= 1 & (ed == st)) {
            amat0i = 0
            amat0_lst[[i]] = amat0i
          } else if (length(zerosi) == 0) {
            #amat0i = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]])
            amat0i = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]], 
                              zeros_ps = zeros_ps, noord_cells = noord_cells[[1]], xvs2=xvs2, grid2=grid2)
            amat0_lst[[i]] = amat0i
          }
        }
        amat0_lst = Filter(Negate(is.null), amat0_lst)
        amat1 = bdiag(amat0_lst)
        amat1 = as.matrix(amat1)
      }else{
        n1 = nc/D1
        amat1_0 = makeamat(1:D1, sh[1], block.ave=block.ave.lst[[1]], block.ord=block.ord.lst[[1]],
                           zeros_ps = zeros_ps, noord_cells = noord_cells[[1]], xvs2=xvs2, grid2=grid2)
        amat1_0_lst = rep(list(amat1_0), n1)
        amat1 = bdiag(amat1_0_lst)
        amat1 = as.matrix(amat1)
      }
    }
    amat = amat1
    #new: for ci
    amat_ci[[1]] = amat1
    #2D
    for (i in 2:(M-1)) {
      Di = Ds[i]
      cump = cumprod(Ds[1:i])
      nci = cump[i]
      gap = cump[i-1]
      ni = nc/nci
      
      if (sh[i] %in% c(9,10)) {
        amati = makeamat_2d_block(x=NULL, sh[i], block.ave = block.ave.lst[[i]], 
                                  block.ord = block.ord.lst[[i]], 
                                  zeros_ps = zeros_ps, 
                                  noord_cells = noord_cells[[i]], 
                                  D_index=i, xvs2=xvs2, grid2=grid2, M_0=M_0)
      } else {
        #new: sh[i] == 0
        if(sh[i] == 0){
          amati = NULL
        }
        if(sh[i] > 0) {
          if ((sh[i] %in% c(1,2))) {
            nr = gap*(Di - 1)
          } 
          if(sh[i] > 2){
            if (Di < 3) {
              stop ("For monotonicity + convexity constraints, number of domains must be >= 3!")
            } else {
              nr = gap*(Di - 2)
            }
          }
          amati_0 = makeamat_2d(sh=sh[i], nr2=nr, nc2=nci, gap=gap, D2=Di, Ds=Ds)
          
          if(length(zeros_ps) == 0){
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
                    #test:
                    if(row_ps > 0){
                      amati_0_cp[row_ps, (col_ps+gap)] = 1
                      amati_0_cp[row_ps, col_ps] = 0
                    }
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
        }
      }  
      amat = rbind(amat, amati)
      #new: for ci
      amat_ci[[i]] = amati
    }
    
    #3D
    cump = cumprod(Ds[1:(M-1)])
    gap = cump[M-1]
    Dm = Ds[M]
    if((sh[M] %in% c(9,10))){
      amatm = makeamat_2d_block(x=NULL, sh[M], block.ave = block.ave.lst[[M]], 
                                block.ord = block.ord.lst[[M]], 
                                zeros_ps = zeros_ps, 
                                noord_cells = noord_cells[[M]], 
                                D_index=M, xvs2=xvs2, grid2=grid2, M_0=M_0)
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
              #test:
              if(row_ps > 0){
                amatm_cp[row_ps, (col_ps+gap)] = 1
                amatm_cp[row_ps, col_ps] = 0
              }
            }
          }
          amatm_cp = amatm_cp[-rm_id, ]
          amatm_cp = amatm_cp[,-zeros_ps]
          amatm = amatm_cp
        }
      }
    }
    amat = rbind(amat, amatm)
    #new: for ci
    amat_ci[[M]] = amatm
  }
  rslt = list(amat = amat, amat_ci = amat_ci)
  return (rslt)
  #return (amat)
}

############################################
#local function to be called in makeamat_2D
############################################
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


################################################################
#local function to be called in makeamat_2D for block ordering
################################################################
makeamat_2d_block = function(x=NULL, sh, block.ave = NULL, block.ord = NULL, zeros_ps = NULL, 
                             noord_cells = NULL, D_index = 1, xvs2 = NULL, 
                             grid2 = NULL, M_0 = NULL) {
  #block average is not added yet
  if (all(block.ave == -1) & (length(block.ord) > 1)) {
    ps = sapply(xvs2[[D_index]], function(.x) which(grid2[, D_index] %in% .x)[1])
    #nbcks is the total number of blocks
    block.ord_nz = block.ord
    noord_ps = which(block.ord == 0)
    #this rm_id is not for empty cells; it is for some domain we don't want to compare ... with
    rm_id = unique(sort(noord_ps))
    #test:
    #z_rm_id = NULL
    #if(length(zeros_ps) > 0){
    #  for(e in zeros_ps){
    #    z_rm_id = c(z_rm_id, max(which(ps <= e)))
    #  }
    #}
    #if(ncol(grid2) == 1){
    #rm_id = unique(sort(c(z_rm_id, noord_ps)))
    #} else {
    #  rm_id = unique(sort(noord_ps))
    #}
    if(length(rm_id)>0){
      block.ord_nz = block.ord[-rm_id]
    }
    
    ubck = unique(block.ord_nz)
    #use values in block.ord as an ordered integer vector
    ubck = sort(ubck)
    
    nbcks = length(table(block.ord_nz))
    szs = as.vector(table(block.ord_nz))
    
    if(length(rm_id) > 0){
      ps = ps[-rm_id]
    }
    #amat dimension: l1*l2...*lk by M
    amat = NULL
    for(i in 1:(nbcks-1)) {
      ubck1 = ubck[i]
      ubck2 = ubck[i+1]
      bck_1 = which(block.ord_nz == ubck1)
      bck_2 = which(block.ord_nz == ubck2)
      l1 = length(bck_1)
      l2 = length(bck_2)
      
      ps_bck1 = ps[which(block.ord_nz == ubck1)]
      ps_bck2 = ps[which(block.ord_nz == ubck2)]
      #for each element in block1, the value for each element in block 2 will be the same
      #for each element in block1, the number of comparisons is l2
      amat_bcki = NULL
      
      #tmp:
      if(D_index == 1){
        M.ord = length(block.ord) 
      }else{
        M.ord = M_0
        #M.ord = nrow(grid2)
      }
      
      amatj = matrix(0, nrow = l2, ncol = M.ord)
      #amatj = matrix(0, nrow = l2, ncol = M.ord-rm_l)
      #bck_1k = bck_1[1]
      bck_1k = ps_bck1[1]
      amatj[, bck_1k] = -1
      row_pointer = 1
      for(j in 1:l2){
        #bck_2j = bck_2[j]
        bck_2j = ps_bck2[j]
        amatj[row_pointer, bck_2j] = 1
        row_pointer = row_pointer + 1
      }
      amatj0 = amatj
      amat_bcki = rbind(amat_bcki, amatj)
      
      if(l1 >= 2){
        for(k in 2:l1) {
          #bck_1k = bck_1[k]; bck_1k0 = bck_1[k-1]
          bck_1k = ps_bck1[k]; bck_1k0 = ps_bck1[k-1]
          amatj = amatj0
          #set the value for the kth element in block 1 to be -1
          #set the value for the (k-1)th element in block 1 back to 0
          #keep the values for block 2
          amatj[, bck_1k] = -1; amatj[, bck_1k0] = 0
          amatj0 = amatj
          amat_bcki = rbind(amat_bcki, amatj)
        }
      }
      amat = rbind(amat, amat_bcki)
    }
    
    if(D_index > 1){
      iter = 1
      amat0 = amat 
      nr = nrow(amat0)
      nc = ncol(amat0)
      ps = sapply(xvs2[[D_index]], function(.x) which(grid2[,D_index] %in% .x)[1])
      ps1_ed = ps[2] - 1
      
      if(length(noord_cells) > 0){
        ps = ps[-noord_cells]
      }
      
      while(iter < ps1_ed) {
        ps1 = ps + iter
        amati = matrix(0, nrow=nr, ncol=nc)
        amati[, ps1] = amat0[, ps]
        amat = rbind(amat, amati)
        iter = iter + 1
      }
    }
    
    #test more...
    #handle empty cells
    if (length(zeros_ps) > 0) {
      compared_pairs = apply(amat, 1, function(e) which(e != 0))
      #check if or not any empty cell has been compared; if yes, then delete the rows
      rm_rows = which(apply(compared_pairs, 2, function (e) any(e %in% zeros_ps)))
      if (length(rm_rows) > 0) {
        amat = amat[-rm_rows, ,drop = FALSE]
      }
      #columns for empty cells must be removed too; there are no thetahat for these domains when run coneA
      amat = amat[, -zeros_ps]
    }
  }
  return (amat)
}

#########################
#empty cell imputation 
##########################
impute_em2 = function(empty_cells, obs_cells, M, yvec, Nd, nd, Nhat, w, domain.ord, grid2, Ds=NULL, sh=NULL, muhat=NULL, lwr=NULL, upp=NULL, vsc=NULL, amat_0=NULL)
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
    #domain.ord_all[which(nd > 0)] = domain.ord
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
    nbors_lst_0 = fn_nbors_3(empty_cells, amat_0, muhat)
    for(k in 1:length(empty_cells)){
      #new: if is new because of rm_id2
      nbors_k_maxs = nbors_lst_0$nb_lst_maxs[[k]]
      nbors_k_mins = nbors_lst_0$nb_lst_mins[[k]]
      em_k = empty_cells[k]
      
      #revised: nbors_k_mins and nbors_k_maxs must be in obs_cells
      #if(length(nbors_k_mins) > 0 | length(nbors_k_maxs) > 0){
      if(length(nbors_k_mins) > 0 & any(nbors_k_mins %in% obs_cells) | length(nbors_k_maxs) > 0 & any(nbors_k_maxs %in% obs_cells)){
        max_ps = min_ps = NULL
        if(length(nbors_k_maxs)>0){
          mu_max = max(muhat[obs_cells %in% nbors_k_maxs])
          max_ps = (nbors_k_maxs[which(muhat[obs_cells %in% nbors_k_maxs] == mu_max)])[1]
        }else{
          #mu_max = muhat[em_k]
          #max_ps = NULL
          mu_min = min(muhat[obs_cells %in% nbors_k_mins])
          min_ps = (nbors_k_mins[which(muhat[obs_cells %in% nbors_k_mins] == mu_min)])[1]
          mu_max = mu_min
          max_ps = min_ps
          dmem = obs_cells[which(obs_cells%in%min_ps)]
        }
        
        if(length(nbors_k_mins)>0){
          mu_min = min(muhat[obs_cells %in% nbors_k_mins])
          min_ps = (nbors_k_mins[which(muhat[obs_cells %in% nbors_k_mins] == mu_min)])[1]
        } else{
          #mu_min = muhat[em_k]
          #min_ps = NULL
          mu_max = max(muhat[obs_cells %in% nbors_k_maxs])
          max_ps = (nbors_k_maxs[which(muhat[obs_cells %in% nbors_k_maxs] == mu_max)])[1]
          mu_min = mu_max
          min_ps = max_ps
          dmem = obs_cells[which(obs_cells%in%max_ps)]
        }
        
        if(!is.null(max_ps)){
          v_max = vsc[which(obs_cells %in% max_ps)]
        }
        if(!is.null(min_ps)){
          v_min = vsc[which(obs_cells %in% min_ps)]
        }
        
        if(!is.null(max_ps) & !is.null(min_ps)){
          varek_new = (mu_max + 2*sqrt(v_max) - mu_min + 2*sqrt(v_min))^2 / 16 
          muhatek = (mu_min - 2*sqrt(v_min) + mu_max + 2*sqrt(v_max))/2
          #?mu_min
          dmem = obs_cells[which(obs_cells%in%min_ps)]
        }
        
        varek = varek_new
        vsce = c(vsce, varek)
        muhate = c(muhate, muhatek)
        domain.ord_em = c(domain.ord_em, dmem)
        
        if (length(nbors_k_mins) == 0) {
          lwrek = -1e+5
          uppek = muhatek + 2*sqrt(varek)
        } else if (length(nbors_k_maxs) == 0) {
          uppek = 1e+5
          lwrek = muhatek - 2*sqrt(varek)
        } else if (length(nbors_k_mins) > 0 & length(nbors_k_maxs) > 0) {
          uppek = muhatek + 2*sqrt(varek)
          lwrek = muhatek - 2*sqrt(varek)
        }
        
        uppe = c(uppe, uppek)
        lwre = c(lwre, lwrek)
        
      } else {
        #ye = c(ye, rev(ye)[1])
        #we = c(we, rev(we)[1])
        #need to test more
        if(!is.null(rev(domain.ord_em)[1])){
          domain.ord_em = c(domain.ord_em, rev(domain.ord_em)[1])
        }
        if (!is.null(muhat)) {
          muhate = c(muhate, rev(muhate)[1])
        }
        if (!is.null(vsc)) {
          vsce = c(vsce, rev(vsce)[1])
          lwre = c(lwre, rev(lwre)[1])
          uppe = c(uppe, rev(uppe)[1])
        }
      }
    }
    
    domain.ord_all[obs_cells] = domain.ord
    #temp:
    #domain.ord_em = obs_cells[which(obs_cells%in%ps[1:length(empty_cells)])]
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
  
  rslt = list(muhat = muhat, lwr = lwr, upp = upp, domain.ord = domain.ord, vsc = vsc)
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
####################################################################################
impute_em3 = function(small_cells, M, yvec, Nd, nd, Nhat, w, domain.ord, grid2, 
                      Ds=NULL, sh=NULL, muhat=NULL, lwr=NULL, upp=NULL, 
                      vsc=NULL, new_obs_cells=NULL, amat_0=NULL)
{
  ne = length(small_cells)
  #observed sample size for n=2...10 cells
  #the `small_cells` are not empty
  nd_ne = nd[small_cells]
  #vsc_ne = vsc[which(new_obs_cells%in%small_cells)]
  vsc_ne = vsc[small_cells]
  
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
    nbors_lst_0 = fn_nbors_2(small_cells, amat_0, muhat)
    for(k in 1:length(small_cells)){
      #new: if is new because of rm_id2
      #nbors_k_maxs = nbors_lst_0$nb_lst_maxs[[k]]
      #nbors_k_mins = nbors_lst_0$nb_lst_mins[[k]]
      max_ps = (nbors_lst_0$nb_lst_maxs[[k]])[1]
      min_ps = (nbors_lst_0$nb_lst_mins[[k]])[1]
      sm_k = small_cells[k]
      
      if(length(max_ps) > 0){      
        mu_max = muhat[max_ps]
        v_max = vsc[max_ps]
      }
      if(length(min_ps) > 0){
        mu_min = muhat[min_ps]
        v_min = vsc[min_ps]
      }
      #new:
      vsc_nek = vsc_ne[k]
      if(length(min_ps) > 0 & length(max_ps) > 0){
        varek_new = (mu_max + 2*sqrt(v_max) - mu_min + 2*sqrt(v_min))^2 / 16 
      }else{
        varek_new = vsc_nek
      }
      #the variance by survey
      #the if else can be replaced by vectorized computation?
      vsce = c(vsce, varek_new)
      #print(c(k, length(vsce)))
    }
    
    if (!is.null(vsc)){
      vsc_all = vsc
      vsc_all[small_cells] = vsce
      vsc = vsc_all  
    }
  }
  #rslt = list(lwr = lwr, upp = upp, vsc = vsc)
  rslt = list(vsc = vsc)
  return(rslt)
}


###############################
#subrountine to find nbors (>=2D)
###############################
fn_nbors_2 = function(small_cells, amat, muhat){
  nr = length(small_cells)
  nb_lst_mins = nb_lst_maxs = vector("list", length = nr)
  #to_merge = list()
  #iter = 1
  for(k in 1:nr){ 
    ptk = small_cells[k]
    col_ptk = amat[, ptk]
    rows_ps = which(col_ptk != 0)
    amat_sub = amat[rows_ps, ,drop=FALSE]
    #nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0)) %>% as.vector() %>% unique() 
    nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0))
    nbors_k = unique(as.vector(nbors_k))
    nbors_k = nbors_k[-which(nbors_k == ptk)]
    #1: min candidates; -1: max candidates
    sgn_nbors_k = col_ptk[which(col_ptk != 0)]
    nbors_k_mins = nbors_k[sgn_nbors_k == 1]
    nbors_k_maxs = nbors_k[sgn_nbors_k == -1]
    nbors_k_min = nbors_k_max = integer(0L)
    if(length(nbors_k_maxs) > 0){
      max_mu = max(muhat[nbors_k_maxs])
      nbors_k_max = nbors_k_maxs[which(muhat[nbors_k_maxs] == max_mu)]
    }
    if(length(nbors_k_mins) > 0){
      min_mu = min(muhat[nbors_k_mins])
      nbors_k_min = nbors_k_mins[which(muhat[nbors_k_mins] == min_mu)]
    }
    
    nb_lst_mins[[k]] = nbors_k_min
    nb_lst_maxs[[k]] = nbors_k_max
  }
  
  rslt = list(nb_lst_maxs=nb_lst_maxs, nb_lst_mins=nb_lst_mins)
  return (rslt)
}

##############################################################
#subrountine to find nbors (>=2D)
#new neighbor routine for >= 2D empty cells
##############################################################
fn_nbors_3 = function(empty_cells, amat_0, muhat){
  nr = length(empty_cells)
  nb_lst_mins = nb_lst_maxs = vector("list", length = nr)
  to_merge_mins = to_merge_maxs = list()
  iter = 1
  for(k in 1:nr){
    nbors_k = NULL
    ptk = empty_cells[k]
    col_ptk = amat_0[, ptk]
    rows_ps = which(col_ptk != 0)
    amat_sub = amat_0[rows_ps, ,drop=FALSE]
    #nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0)) %>% as.vector() %>% unique() 
    nbors_k = apply(amat_sub, 1, function(.x) which(.x != 0))
    nbors_k = unique(as.vector(nbors_k))
    nbors_k = nbors_k[-which(nbors_k == ptk)]
    #1: min candidates; -1: max candidates
    sgn_nbors_k = col_ptk[which(col_ptk != 0)]
    nbors_k_mins = nbors_k[sgn_nbors_k == 1]
    nbors_k_maxs = nbors_k[sgn_nbors_k == -1]
    
    #new: need to remove empty_cells from nbors_k_maxs and nbors_k_mins
    nbors_k_mins = setdiff(nbors_k_mins, empty_cells)
    nbors_k_maxs = setdiff(nbors_k_maxs, empty_cells)
    
    if(length(nbors_k_maxs) > 0){
      nbors_k_maxs = nbors_k_maxs[which(muhat[nbors_k_maxs] == max(muhat[nbors_k_maxs]))]
    }
    if(length(nbors_k_mins) > 0){
      nbors_k_mins = nbors_k_mins[which(muhat[nbors_k_mins] == min(muhat[nbors_k_mins]))]
    }
    
    #new: label the cells whose neighbors should be merged later
    rm_id2 = NULL
    if (length(nbors_k_maxs) > 0) {
      for(h in 1:length(nbors_k_maxs)){
        nbh = nbors_k_maxs[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          ps = which(empty_cells == nbh)
          to_merge_maxs[[iter]] = sort(c(k, ps))
          iter = iter + 1
        }
      }
    }
    
    if (length(nbors_k_mins) > 0) {
      for(h in 1:length(nbors_k_mins)){
        nbh = nbors_k_mins[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          ps = which(empty_cells == nbh)
          to_merge_mins[[iter]] = sort(c(k, ps))
          iter = iter + 1
        }
      }
    }
    
    nb_lst_mins[[k]] = nbors_k_mins 
    nb_lst_maxs[[k]] = nbors_k_maxs
  }
  
  #need to be more efficient
  to_merge_maxs = unique(to_merge_maxs)
  to_merge_mins = unique(to_merge_mins)
  nm_maxs = length(to_merge_maxs)
  nm_mins = length(to_merge_mins)
  if (nm_maxs > 1) {
    for(i in 1:(nm_maxs-1)){
      to_merge_i = to_merge_maxs[[i]]
      vec = to_merge_i
      vec_ps = i
      for(j in (i+1):nm_maxs){
        to_merge_j = to_merge_maxs[[j]]
        if(any(to_merge_j %in% to_merge_i)){
          vec = sort(unique(c(vec, to_merge_j)))
          vec_ps = c(vec_ps, j)
        }
      }
      if (length(vec_ps) > 1) {
        to_merge_maxs[vec_ps] = lapply(to_merge_maxs[vec_ps], function(x) x = vec)
      }
    }
    to_merge_maxs = unique(to_merge_maxs)
  }
  
  if (nm_maxs > 0) {
    for(i in 1:length(to_merge_maxs)) {
      ps = to_merge_maxs[[i]]
      #nbors_ps_maxs = unlist(nb_lst_maxs[ps]) %>% unique()
      nbors_ps_maxs = unique(unlist(nb_lst_maxs[ps]))
      rm_id2_maxs = NULL
      for(h in 1:length(nbors_ps_maxs)) {
        nbh = nbors_ps_maxs[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_maxs = c(rm_id2_maxs, h)
        }
      }
      if (length(rm_id2_maxs) > 0) {
        nbors_ps_maxs = nbors_ps_maxs[-rm_id2_maxs]
      }
      nb_lst_maxs[ps] = lapply(nb_lst_maxs[ps], function(x) x = nbors_ps_maxs)
      
      ###
      for(h in 1:length(nbors_ps_maxs)) {
        nbh = nbors_ps_maxs[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_maxs = c(rm_id2_maxs, h)
        }
      }
      if (length(rm_id2_maxs) > 0) {
        nbors_ps_maxs = nbors_ps_maxs[-rm_id2_maxs]
      }
      nb_lst_maxs[ps] = lapply(nb_lst_maxs[ps], function(x) x = nbors_ps_maxs)
    }
  }
  
  if (nm_mins > 0) {
    for(i in 1:length(to_merge_mins)) {
      ps = to_merge_mins[[i]]
      #nbors_ps_mins = unlist(nb_lst_mins[ps]) %>% unique()
      nbors_ps_mins = unique(unlist(nb_lst_mins[ps]))
      rm_id2_mins = NULL
      for(h in 1:length(nbors_ps_mins)) {
        nbh = nbors_ps_mins[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_mins = c(rm_id2_mins, h)
        }
      }
      if (length(rm_id2_mins) > 0) {
        nbors_ps_mins = nbors_ps_mins[-rm_id2_mins]
      }
      nb_lst_mins[ps] = lapply(nb_lst_mins[ps], function(x) x = nbors_ps_mins)
      
      ###
      for(h in 1:length(nbors_ps_mins)) {
        nbh = nbors_ps_mins[h]
        bool = any(empty_cells %in% nbh)
        if (bool) {
          rm_id2_mins = c(rm_id2_mins, h)
        }
      }
      if (length(rm_id2_mins) > 0) {
        nbors_ps_mins = nbors_ps_mins[-rm_id2_mins]
      }
      nb_lst_mins[ps] = lapply(nb_lst_mins[ps], function(x) x = nbors_ps_mins)
    }
  }
  
  rslt = list(nb_lst_maxs=nb_lst_maxs, nb_lst_mins=nb_lst_mins)
  return (rslt)
}


############
#fitted.csvy
############
fitted.csvy = function(object,...)
{
  ans = object$muhat
  ans
}


##############
#confint.csvy
##############
confint.csvy = function(object, parm = NULL, level = 0.95, type = c("link", "response"),...) {
  type = match.arg(type)
  n.mix = object$n.mix
  #if (n.mix %% 1 != 0 | n.mix <= 0){
  #  stop("confint function only works when n.mix is a positive integer!")
  #} else {
  z.mult = qnorm((1 - level) / 2, lower.tail = FALSE)
  # vs = vcov(object)
  vs = object$cov.unscaled
  #lwr = object$etahat - z.mult * sqrt(vs)
  #upp = object$etahat + z.mult * sqrt(vs)
  lwr = object$lwr
  upp = object$upp
  
  lwr = switch(type,
                link = lwr,
                response = object$family$linkinv(lwr)
  )
  upp = switch(type,
                link = upp,
                response = object$family$linkinv(upp)
  )
  
  ci.bands = as.data.frame(do.call(cbind, list(lwr, upp)))
  # learn from confint.svyglm:
  a = (1 - level) / 2
  a = c(a, 1 - a)
  pct = paste0(format(100 * a,
                       trim = TRUE,
                       scientific = FALSE, digits = 3
  ), "%")
  colnames(ci.bands) = pct
  return(ci.bands)
  #}
}

###############
#print method
###############
print.csvy = function(x,...){
  print(x$survey.design, varnames = FALSE, design.summaries = FALSE,...)
  NextMethod()
}

###################################################################################
#extract mixture variance-covariance matrix
###################################################################################
vcov.csvy = function(object,...){
  #v = object$acov_cp
  v = object$acov
  #temp: with empty cells, off-diagonal values in v are not imputed
  if (all(object$nd > 0)) {
    dimnames(v) = list(names(coef(object)), names(coef(object)))
  } else {
    nv = dim(v)[1]
    vec_nm = as.character(1:nv)
    dimnames(v) = list(vec_nm, vec_nm)
  }
  return(v)
}


#############################################
#summary for csvy
#I don't remember why I changed this summary function
#use the one in _0517 instead
#############################################
# summary.csvy = function(object, ...) {
#   family = object$family
#   CIC = object$CIC
#   CIC.un = object$CIC.un
#   deviance = object$deviance
#   null.deviance = object$null.deviance
#   bval = object$bval
#   pval = object$pval
#   # edf
#   df.residual = object$df.residual
#   df.null = object$df.null
#   edf = object$edf
#   tms = object$terms
#   rslt = NULL
#   # learn from summary.svyglm
#   # not work for gaussian; svyby
#   #dispersion = svyvar(resid(object,"pearson"), object$survey.design, na.rm = TRUE)
#   #disperson = survey:::summary.svrepglm(object)$disperson
#   if(inherits(object, "svrepglm")){
#     presid = resid(object,"pearson")
#     dispersion = sum(object$survey.design$pweights*presid^2, na.rm = TRUE) / sum(object$survey.design$pweights)
#   } else if (inherits(object, "svyglm")){
#     dispersion = svyvar(resid(object, "pearson"), object$survey.design, na.rm = TRUE)
#   }
#   if (!is.null(pval)) {
#     # rslt = data.frame("edf" = round(resid_df_obs, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
#     # same as cgam:
#     rslt = data.frame("edf" = round(edf, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
#     # check?
#     rownames(rslt) = rev((attributes(tms)$term.labels))[1]
#     structure(list(
#       call = object$call, coefficients = rslt, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
#       null.deviance = null.deviance, deviance = deviance, family = family, dispersion = dispersion
#     ), class = "summary.csvy")
#   } else {
#     structure(list(
#       call = object$call, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
#       null.deviance = null.deviance, deviance = deviance, family = family, dispersion = dispersion
#     ), class = "summary.csvy")
#     # warning("summary function for a csvy object only works when test = TRUE!")
#   }
# }


summary.csvy <- function(object,...) {
  family <- object$family
  CIC <- object$CIC
  CIC.un <- object$CIC.un
  deviance <- object$deviance
  null.deviance <- object$null.deviance
  bval <- object$bval
  pval <- object$pval
  #edf
  df.residual <- object$df.residual
  df.null <- object$df.null
  edf <- object$edf
  tms <- object$terms
  rslt <- NULL
  if (!is.null(pval)) {
    #rslt = data.frame("edf" = round(resid_df_obs, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
    #same as cgam:
    rslt <- data.frame("edf" = round(edf, 4), "mixture of Beta" = round(bval, 4), "p.value" = round(pval, 4))
    #check?
    rownames(rslt) <- rev((attributes(tms)$term.labels))[1]
    structure(list(call = object$call, coefficients = rslt, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
                   null.deviance = null.deviance, deviance = deviance, family = family), class = "summary.csvy")
    
    #structure(list(call = object$call, coefficients = rslt, CIC = CIC, df.residual = df.residual,
    #               null.deviance = null.deviance, deviance = deviance, family = family), class = "summary.csvy")
    #ans = list(call = x$call, coefficients = rslt, CIC = CIC, 
    #           null.deviance = null.deviance, deviance = deviance, 
    #           family = family)
    #class(ans) = "summary.csvy"
    #ans
  } else {
    structure(list(call = object$call, CIC = CIC, CIC.un = CIC.un, df.null = df.null, df.residual = df.residual,
                   null.deviance = null.deviance, deviance = deviance, family = family), class = "summary.csvy")
    #warning("summary function for a csvy object only works when test = TRUE!")
  }
}

#####################
#print.summary.csvy#
#####################
print.summary.csvy = function(x,...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Null deviance: ", round(x$null.deviance, 4), "", "on", x$df.null, "", "degrees of freedom", "\n")
  cat("Residual deviance: ", round(x$deviance, 4), "", "on", x$df.residual, "", "observed degrees of freedom", "\n")
  if (!is.null(x$coefficients)) {
    cat("\n")
    cat("Approximate significance of constrained fit: \n")
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  }
  if (!is.null(x$CIC)) {
    cat("CIC (constrained estimator): ", round(x$CIC, 4))
    cat("\n")
    cat("CIC (unconstrained estimator): ", round(x$CIC.un, 4))
  }
}


##############
#predict.csvy#
##############
predict.csvy = function(object, newdata=NULL, type=c("link","response"), 
                        se.fit=TRUE, level=0.95, n.mix=100,...)
{
  family = object$family
  type = match.arg(type)
  if (!inherits(object, "csvy")) { 
    warning("calling predict.csvy(<fake-csvy-object>) ...")
  }
  #match newdata with grid
  xnms_add = object$xnms_add
  grid2 = object$grid
  colnames(grid2) = xnms_add
  #add colnames in the object?
  if (missing(newdata) || is.null(newdata)) {
    newdata = grid2
  }
  if (!is.data.frame(newdata)) {
    stop ("newdata must be a data frame!")	
  }
  #learnt from predict.svyglm:
  #if (type=="terms")
  #  return(predterms(object,se=se.fit,...))
  #make it simpler?
  if(ncol(newdata) == 1){
    ps = apply(newdata, 1, function(elem) grid2[which(apply(grid2, 1, function(gi) all(gi == elem))),])
  }
  if(ncol(newdata) > 1){
    ps0 = sapply(colnames(newdata), function(elem) which(sapply(colnames(grid2), function(gi) all(gi == elem))))
    #order the column in case the user doesn't define the variables as the order used in the formula
    newdata = newdata[, as.numeric(ps0)]
    ps = apply(newdata, 1, function(elem) which(apply(grid2, 1, function(gi) all(gi == elem))))
  }
  
  etahat = object$etahat 
  vsc_mix = object$cov.unscaled
  
  #interval = match.arg(interval)
  if(!se.fit) {
    fit = etahat[ps]
    fit = switch(type, link = fit, response = object$family$linkinv(fit))
    return (fit)
  } else {
    #when n.mix = 0, then there is no vsc_mix
    if(is.null(vsc_mix)) {
      ynm = object$ynm
      amat = object$amat
      nd = object$nd
      M = length(nd)
      muhat = object$muhat
      #add in csvy
      #revise later: etahat.s
      etahatu = object$etahatu
      w = object$w
      v1 = object$cov.un
      domain.ord = object$domain.ord
      #move the imputation code to be in predict.csvy
      xm.red = object$xm.red
      ne = object$ne
      xvs2 = apply(xm.red, 2, function(x) 1:length(unique(x)))
      zeros_ps = empty_cells = object$empty_cells
      small_cells = object$small_cells
      obs_cells = object$obs_cells
      nsm = length(small_cells)
      zeros = object$zeros
      Nhat = object$Nhat
      Ds = object$Ds
      sh = object$shapes_add
      amat_0 = object$amat_0
      Nd = object$Nd
      
      lwr = upp = NULL
      acov = v1
      dp = -t(amat)
      dp = apply(dp, 2, function(e) e / (sum(e^2))^(.5))
      
      M = M - zeros
      m_acc = M
      sector = NULL
      times = NULL
      df.face = NULL
      iter = 1
      obs = 1:M
      #ysims = MASS::mvrnorm(n.mix, mu=muhat, Sigma=v1)
      n.mix = 100
      #binomial can use mvrnorm?
      #etahat returned from object has been imputed
      if(length(empty_cells) > 0){
        etahat = etahat[-empty_cells]
      }
      ysims = MASS::mvrnorm(n.mix, mu=etahat, Sigma=v1)
      for (iloop in 1:n.mix) {
        #ysim is the group mean
        #new:
        ysim = ysims[iloop, ]
        ansi = coneA(ysim, amat, w=w, msg=FALSE)
        etahati = round(ansi$thetahat, 10)
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
        imat = diag(M)
        sm_id = which((df.face[,2]/n.mix) < 1e-3)
        if (any(sm_id)) {
          df.face = df.face[-sm_id, ,drop=FALSE]
          sector = sector[-sm_id, ,drop=FALSE]
        }
        nsec = nrow(df.face)
        bsec = df.face
        bsec[,2] = bsec[,2] / sum(bsec[,2])
        
        acov = matrix(0, nrow=M, ncol=M)
        wtinv = diag(1/w)
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
      lwr = loweta = etahat - hl
      upp = uppeta = etahat + hl
      
      vscu = diag(v1)
      hlu = z.mult*sqrt(vscu)
      lwru = lowetau = etahatu - hlu
      uppu = uppetau = etahatu + hlu
      
      if(is.matrix(xvs2)){
        grid2 = expand.grid(as.data.frame(xvs2))
      }else if(is.list(xvs2)){
        grid2 = expand.grid(xvs2)
      }
      
      if (ne >= 1) {
        #muhatu = yvecu
        #new: if there's empty cells, just augment muhatu to include NA's
        etahatu_all = 1:(M+zeros)*0
        etahatu_all[zeros_ps] = NA
        etahatu_all[obs_cells] = etahatu
        etahatu = etahatu_all
        
        ans_im = impute_em2(empty_cells, obs_cells, M=(M+zeros), etahatu, Nd, nd, Nhat, w,
                            domain.ord, grid2, Ds, sh, etahat, lwr, upp, vsc, amat_0)
        #ans_im = with(object, impute_em2(empty_cells, obs_cells, M=(M+zeros), etahatu, Nd, nd, Nhat, w,
        #                    domain.ord, grid2, Ds, sh, etahat, lwr, upp, vsc, amat_0))
        etahat = ans_im$muhat
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
          ans_im2 = impute_em2(small_cells, new_obs_cells, M=(M+zeros), etahatu,
                               Nd, nd, Nhat, w, domain.ord, grid2, Ds, sh, etahat,
                               lwr, upp, vsc, amat_0)
          #ans_im2 = with(object, impute_em2(small_cells, new_obs_cells, M=(M+zeros), etahatu,
          #         Nd, nd, Nhat, w, domain.ord, grid2, Ds, sh, etahat, lwr, upp, vsc, amat_0))
          sig2 = ans_im2$vsc
        }
        if(Nd >= 2){
          ans_im2 = impute_em3(small_cells, M=(M+zeros), etahatu, Nd, nd, Nhat, w, domain.ord,
                               grid2, Ds, sh, etahat, lwr, upp, vsc, new_obs_cells, amat_0)
          sig2 = ans_im2$vsc
        }
      }
      
      vsc_mix = sig1
      vsc_mix_imp = NULL
      if (Nd >= 1 & nsm >= 1) {
        imp_ps = sort(c(zeros_ps, small_cells))
        nd_imp_ps = nd[imp_ps]
        l_imp_ps = length(nd_imp_ps)
        
        sig1_imp = sig1[imp_ps]
        sig2_imp = sig2[imp_ps]
        
        vsc_mix_imp = 1:l_imp_ps*0
        #tmp:
        if(Nd > 1){
          nbors_lst_0 = fn_nbors_2(small_cells, amat_0, muhat)
        }
        
        k_sm = 1
        for(k in 1:l_imp_ps){
          ndek = nd_imp_ps[k]
          s1k = sig1_imp[k]
          s2k = sig2_imp[k]
          #tmp
          if(Nd > 1 & ndek > 0){
            max_ps = nbors_lst_0$nb_lst_maxs[[k_sm]]
            min_ps = nbors_lst_0$nb_lst_mins[[k_sm]]
            #max_ps = min_ps = NULL
            if(length(max_ps)>0){
              mu_max = muhat[max_ps]
            }
            if(length(min_ps)>0){
              mu_min = muhat[min_ps]
            }
            k_sm = k_sm + 1
          } else if (Nd == 1 | ndek == 0) {
            min_ps = max_ps = 1
          }
          
          #when min_ps or max_ps = NULL: kth cell is in a 2D case, small cell, at the edge of the
          #grid and there is only one constraint in amat
          if(length(min_ps) == 0 | length(max_ps) == 0) {
            varek = s1k
          } else {
            if(ndek == 0){
              #varek = s1k
              varek = s2k
            } else if (ndek == 1) {
              varek = s1k*.2 + s2k*.8
            } else if(ndek == 2 | ndek == 3){
              varek = s1k*.3 + s2k*.7
            }else if(ndek == 4 | ndek == 5){
              varek = s1k*.6 + s2k*.4
            }else if(ndek == 6 | ndek == 7){
              varek = s1k*.6 + s2k*.4
            }else if(ndek == 8 | ndek == 9 | ndek == 10){
              varek = s1k*.8 + s2k*.2
            }#else if(ndek == 10){
            #  varek = s1k*.8 + s2k*.2
            #}
          }
          vsc_mix_imp[k] = varek
        }
        vsc_mix[imp_ps] = vsc_mix_imp
      }
    }
    #eta = etahat[ps]; se = vsc[ps]
    #learnt from predict.svyglm:
    z.mult = qnorm((1 - level)/2, lower.tail=FALSE)
    fit = etahat[ps]
    se.fit = sqrt(vsc_mix[ps])
    lwr = etahat[ps] - z.mult*se.fit
    upp = etahat[ps] + z.mult*se.fit
    
    #new: make a copy; not change endpoints' c.i.
    lwr_0 = lwr
    upp_0 = upp
    #new: learn from cgam; nasa's idea
    #for now, only handles one predictor
    sh = object$shapes_add
    if (length(sh) == 1) {
      nci = length(lwr)
      if (sh == 1) #incr
      {
        for(i in (nci - 1):1)
        {
          if(upp[i] > upp[i + 1])
          {
            upp[i] = upp[i + 1]
          }
        }
        for(i in 1:(nci - 1))
        {
          if(lwr[i + 1] < lwr[i])
          {
            lwr[i + 1] = lwr[i]
          }
        } 
        
        #new: don't change endpoints' ci!
        lwr[1] = lwr_0[1]
        upp[1] = upp_0[1]
        lwr[nci] = lwr_0[nci]
        upp[nci] = upp_0[nci]
      }
      
      if (sh == 2) #decr
      {
        for(i in (nci - 1):1)
        {
          if(lwr[i] < lwr[i + 1])
          {
            lwr[i] = lwr[i + 1]
          }
        }
        
        for(i in 1:(nci - 1))
        {
          if(upp[i + 1] > upp[i])
          {
            upp[i + 1] = upp[i]
          }
        }
        
        #new: don't change endpoints' ci!
        lwr[1] = lwr_0[1]
        upp[1] = upp_0[1]
        lwr[nci] = lwr_0[nci]
        upp[nci] = upp_0[nci]
      }
    }
    
    fit = switch(type, link = fit, response = object$family$linkinv(fit))
    lwr = switch(type, link = lwr, response = object$family$linkinv(lwr))
    upp = switch(type, link = upp, response = object$family$linkinv(upp))
    se.fit = switch(type, link = se.fit, response = object$family$linkinv(se.fit))
    
    #no unconstrained:...
    if(type=='link'){
      ans = list(fit = fit, lwr = lwr, upp = upp, se.fit = se.fit)
    }
    if(type=='response'){
      ans = list(fit = fit, lwr = lwr, upp = upp)
    }
  }
  return (ans)
}


#---------------------------------------------------------------------------------------------------------------------------------------
#routines modified from methods of svyby or svyglm, svyrepglm
#---------------------------------------------------------------------------------------------------------------------------------------
# dotchart.csvy = function(x,...,pch = 19){
#   xx = x$ans.unc_cp
#   dotchart(xx,...,pch = pch)
# }

#---------------------------------------------------------------------------------------------------------------------------------------
deff.csvy = function(object,...) {
  x = object$ans.unc_cp
  deff(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
barplot.csvy = function(height, beside = TRUE,...){
  x = height$ans.unc_cp
  #survey:::barplot.svyby(xx, beside = TRUE,...)
  barplot(x, beside = TRUE,...)
}

#plot.csvy = barplot.csvy
#---------------------------------------------------------------------------------------------------------------------------------------
SE.csvy = function(object,...){
  x = object$ans.unc_cp
  SE(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
coef.csvy = function(object, ...) {
  x = object$ans.unc_cp
  coef(x,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
svycontrast.csvy = function(stat, contrasts,...) {
  x = stat$ans.unc_cp
  svycontrast(x, contrasts,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#fix
ftable.csvy = function(x,...) {
  xx = x$ans.unc_cp
  ftable(xx,...)
}

#---------------------------------------------------------------------------------------------------------------------------------------
# confint.csvyby = function(object, parm, level = 0.95, df = Inf,...){
#   x = object$ans.unc_cp
#   confint(x, parm, level, df,...)
# }

#---------------------------------------------------------------------------------------------------------------------------------------
# deff.csvyby = function(object,...){
#   x = object$ans.unc_cp
#   if(missing(x) | is.null(x)){
#     stop("only works when the fit is a class of svyby!")
#   }
#   deff(x,...)
# }

