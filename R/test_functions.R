



# RMST.method.test
# input:
# time, status, group:    vectors of the same length
# OR
# formula, event, data:   formula, event variable name and data set

# hyp_mat:                partitionized matrix as list, i.e. list of length L with matrices as entries
#                         alternatively, "center", "Dunnett", "Tukey" as 'type' in ?contr_mat can be used
# hyp_vec:                partitionized vector as list (default = zero vector)
# tau

# output:
# Object of class "GFDrmst"
# method: character of used method
# test_stat: vector of test statistics
# p.value: vector of adjusted p values
# res: list with hypothesis matrices, estimators, confidence intervals if the
#       matrices are vectors,  test statistics, critical values, decisions
# alpha: level of significance

RMST.asymptotic.test <- function(time = NULL, status = NULL, group = NULL,
                                 formula = NULL, event = NULL, data = NULL,
                                 hyp_mat, hyp_vec = NULL,
                                 tau, stepwise = FALSE, alpha = 0.05,
                                 Nres = 4999, seed = 1){

  set.seed(seed)

  if(is.null(time) | is.null(status) | is.null(group)){
    dat2 <- formula2input(formula, event, data)
    time   <- dat2$time
    status <- dat2$event
    group  <- dat2$group
  }

  k <- length(unique(group))
  group <- factor(group, labels = 1:k)
  # some possible values for hyp_mat
  if(is.character(hyp_mat)){
    if(hyp_mat %in% c("center", "Dunnett", "Tukey")){
      Y <- GFDmcv::contr_mat(k, type = hyp_mat)
      X <- data.frame(t(Y))
      colnames(X) <- rownames(Y)
      hyp_mat <- as.list(X)
      hyp_mat <- lapply(hyp_mat, t)
    }else{
      if(hyp_mat== "crossed factorial"){
        if(is.null(formula)) stop("Formula needed to realize the crossed factorial design.")
        formula2 <- paste0(formula, "*", event)
        dat <- model.frame(formula2, data)
        nadat <- names(dat)
        nf <- ncol(dat) - 2
        fl <- NA
        for (aa in 1:nf) {
          fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
        }
        levels <- list()
        for (jj in 1:nf) {
          levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
        }
        lev_names <- expand.grid(levels)
        if (nf == 1) {
          hyp_mat <- list(diag(fl) - matrix(1/fl, ncol = fl,
                                            nrow = fl))
        }else {
          lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),]
          formula3 <- as.formula(formula)
          nr_hypo <- attr(terms(formula3), "factors")
          fac_names_original <- colnames(nr_hypo)
          formula1 <- as.formula( gsub("\\+", "*", formula))
          nr_hypo <- attr(terms(formula1), "factors")
          fac_names <- colnames(nr_hypo)
          perm_names <- t(attr(terms(formula1), "factors")[-1, ])

          hyp_mat <- HC(fl, perm_names, fac_names)[[1]]
          names(hyp_mat) <- HC(fl, perm_names, fac_names)[[2]]
          hyp_mat <- hyp_mat[fac_names_original]

        }
      }else{ stop("hyp_mat is an invalid character input") }
    }}else{if(is.null(hyp_vec)){
      warning("hyp_vec is chosen as zero vector since it is unspecified!")
    }}
  if(is.null(hyp_vec)){
    hyp_vec <- lapply(hyp_mat, function(x) rep(0, nrow(x)))
    #warning("hyp_vec is chosen as zero vector since it is unspecified!")
  }
  if(is.matrix(hyp_mat)) hyp_mat <- list(hyp_mat)
  if(!is.list(hyp_vec)) hyp_vec <- list(hyp_vec)
  L <- length(hyp_mat)

  # check whether all hypotheses have a solution in [0,tau]^k
  solution <- logical(length(hyp_mat))
  for(i in 1:length(hyp_mat)){
    solution[i] <- existence(hyp_mat[[i]], hyp_vec[[i]], tau)
  }
  if(any(!solution)){
    stop(paste0("Hypothesis ", which(!solution),
                " does not have a possible solution in [0,tau]^k. "))
  }

  rmst <- numeric(k)
  var <- numeric(k)
  teststats <- numeric(L)

  attr(teststats, "names") <- paste0("W_n(H_", 1:L, ", c_", 1:L, ")")

  status <- factor(status, levels = c(0,1))
  #time <- factor(time, exclude = c(NA, NaN))
  values <- data.frame(time=time, status=status, group=group)

  for(i in 1:k) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }

  Sigma_hat <- diag(var)

  # calculate the test statistics
  for(i in 1:L){
    vec <- hyp_mat[[i]] %*% rmst - hyp_vec[[i]]
    teststats[i] <- t(vec) %*% MASS::ginv(hyp_mat[[i]] %*% Sigma_hat %*% t(hyp_mat[[i]])) %*% vec
    if(teststats[i] == 0) teststats[i] <- ifelse(all(vec == 0), 0, Inf)
  }

  namen <- names(hyp_mat)
  if(is.null(namen)) namen <- 1:length(hyp_mat)

  # calculate the (adjusted) p-values
  # if all matrices are vectors
  if(all(sapply(hyp_mat, function(x) nrow(x)) == 1)){
    # calculate the scaling matrices for all contrast matrices
    D <- diag(sapply(hyp_mat, function(vec){
      MASS::ginv(sqrt(vec %*% Sigma_hat %*% t(vec))) #using the g-inverse to avoid problems with zero variances
    }), nrow = L, ncol = L)
    # calculate the quantiles
    H <- global_mat(hyp_mat, k)
    # equicoordinate normal quantiles
    sigma <- (D%*%H%*%Sigma_hat%*%t(H)%*%D)
    # since errors occur when any of diag(sigma) = 0, we do not consider these components
    my_index <- (diag(sigma) != 0)
    if(sum(my_index)>0){
      quant <- mvtnorm::qmvnorm(1-alpha, sigma=sigma[my_index, my_index], tail = "both.tails")$quantile
    }else{ quant <- 0 }
    pvalues <- numeric(L)

    #out <- list()
    res <- matrix(list() , nrow = L, ncol = 8)
    colnames(res) <- c("hyp_matrix", "estimator", "lwr_conf", "upr_conf", "test_stat", "critical value", "adj_pvalue", "decision")
    rownames(res) <- names(hyp_mat)

    for(l in 1:L){
      if(sum(my_index) > 0){
        bound <- rep(sqrt(teststats[l]), sum(my_index))
        pvalues[l] <- (1-mvtnorm::pmvnorm(lower = -bound, upper = bound, sigma=sigma[my_index, my_index]))
        # calculate the simultaneous confidence intervals
        conf.int <- c(hyp_mat[[l]]%*%rmst) + sqrt(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]]))) * quant * t(c(-1,1))
        # determine the smallest and largest possible value
        conf.int[1] <- max(conf.int[1], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau))
        conf.int[2] <- min(conf.int[2], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0))
        attr(conf.int, "conf.level") <- 1-alpha
      }else{
        pvalues[l] <- as.numeric(hyp_mat[[l]]%*%rmst == hyp_vec[[l]])
        conf.int <- t(rep(hyp_mat[[l]]%*%rmst, 2))
        attr(conf.int, "conf.level") <- 1-alpha
      }
      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
      res[[l,3]] <- conf.int[1]
      res[[l,4]] <- conf.int[2]
      res[[l,5]] <- teststats[l]
      res[[l,6]] <- quant^2
    }

    if(stepwise & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMST.asymptotic.test(time, status, group, hyp_mat = hyp_mat[ind],
                                 hyp_vec = hyp_vec[ind], tau = tau,
                                 stepwise = FALSE, alpha = alpha,
                                 Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the confidence interval and critical value
      res <- res[,-c(3,4,6)]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

    if(stepwise) method <- "Multiple asymptotic RMST Wald-type tests with stepwise extension"
    if(!stepwise) method <- "Multiple asymptotic RMST Wald-type tests"

    out <- list(method = method,
                test_stat = teststats,
                p.value = pvalues,
                res = res,
                alpha = alpha)
    class(out) <- "GFDrmst"

    return(out)
  }else{
    # generate Nres multivariate random vectors
    random_numbers <- mvtnorm::rmvnorm(Nres , sigma=Sigma_hat)
    # determine values for approximating the limiting distribution
    random_values <- t(sapply(hyp_mat, function(mat) apply(random_numbers, 1, function(z) t(mat%*%z)%*%MASS::ginv(mat%*%Sigma_hat%*%t(mat))%*%mat%*%z)))
    pvalues <- p.value(data_mat =  random_values, teststat=teststats)

    res <- matrix(list() , nrow = L, ncol = 6)
    colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")
    rownames(res) <- names(hyp_mat)
    cv <- crit_values2(random_values, alpha = alpha)

    for(l in 1:L){
      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
      res[[l,3]] <- teststats[l]
      res[[l,4]] <- cv[l]

    }

    if(stepwise & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMST.asymptotic.test(time, status, group, hyp_mat = hyp_mat[ind],
                                 hyp_vec = hyp_vec[ind], tau = tau,
                                 stepwise = FALSE, alpha = alpha,
                                 Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the critical value
      res <- res[,-4]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

    if(stepwise) method <- "Multiple asymptotic RMST Wald-type tests with stepwise extension"
    if(!stepwise) method <- "Multiple asymptotic RMST Wald-type tests"

    out <- list(method = method,
                test_stat = teststats,
                p.value = pvalues,
                res = res,
                alpha = alpha)
    class(out) <- "GFDrmst"

    return(out)
  }
}


RMST.groupwise.test <- function(time = NULL, status = NULL, group = NULL,
                                formula = NULL, event = NULL, data = NULL,
                                hyp_mat, hyp_vec = NULL,
                                tau, stepwise = FALSE, alpha = 0.05,
                                Nres = 4999, seed = 1){
  set.seed(seed)

  if(is.null(time) | is.null(status) | is.null(group)){
    dat2 <- formula2input(formula, event, data)
    time   <- dat2$time
    status <- dat2$event
    group  <- dat2$group
  }

  k <- length(unique(group))
  group <- factor(group, labels = 1:k)
  # some possible values for hyp_mat
  if(is.character(hyp_mat)){
    if(hyp_mat %in% c("center", "Dunnett", "Tukey")){
      Y <- GFDmcv::contr_mat(k, type = hyp_mat)
      X <- data.frame(t(Y))
      colnames(X) <- rownames(Y)
      hyp_mat <- as.list(X)
      hyp_mat <- lapply(hyp_mat, t)
    }else{
      if(hyp_mat== "crossed factorial"){
        if(is.null(formula)) stop("Formula needed to realize the crossed factorial design.")
        formula2 <- paste0(formula, "*", event)
        dat <- model.frame(formula2, data)
        nadat <- names(dat)
        nf <- ncol(dat) - 2
        fl <- NA
        for (aa in 1:nf) {
          fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
        }
        levels <- list()
        for (jj in 1:nf) {
          levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
        }
        lev_names <- expand.grid(levels)
        if (nf == 1) {
          hyp_mat <- list(diag(fl) - matrix(1/fl, ncol = fl,
                                            nrow = fl))
        }else {
          lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),]
          formula3 <- as.formula(formula)
          nr_hypo <- attr(terms(formula3), "factors")
          fac_names_original <- colnames(nr_hypo)
          formula1 <- as.formula( gsub("\\+", "*", formula))
          nr_hypo <- attr(terms(formula1), "factors")
          fac_names <- colnames(nr_hypo)
          perm_names <- t(attr(terms(formula1), "factors")[-1, ])

          hyp_mat <- HC(fl, perm_names, fac_names)[[1]]
          names(hyp_mat) <- HC(fl, perm_names, fac_names)[[2]]
          hyp_mat <- hyp_mat[fac_names_original]

        }
      }else{ stop("hyp_mat is an invalid character input") }
    }}else{if(is.null(hyp_vec)){
      warning("hyp_vec is chosen as zero vector since it is unspecified!")
    }}
  if(is.null(hyp_vec)){
    hyp_vec <- lapply(hyp_mat, function(x) rep(0, nrow(x)))
    #warning("hyp_vec is chosen as zero vector since it is unspecified!")
  }
  if(is.matrix(hyp_mat)) hyp_mat <- list(hyp_mat)
  if(!is.list(hyp_vec)) hyp_vec <- list(hyp_vec)
  L <- length(hyp_mat)
  n_vec <- table(group)

  # check whether all hypotheses have a solution in [0,tau]^k
  solution <- logical(length(hyp_mat))
  for(i in 1:length(hyp_mat)){
    solution[i] <- existence(hyp_mat[[i]], hyp_vec[[i]], tau)
  }
  if(any(!solution)){
    stop(paste0("Hypothesis ", which(!solution),
                " does not have a possible solution in [0,tau]^k. "))
  }

  rmst <- numeric(k)
  var <- numeric(k)
  teststats <- numeric(L)

  attr(teststats, "names") <- paste0("W_n(H_", 1:L, ", c_", 1:L, ")")

  status <- factor(status, levels = c(0,1))
  #time <- factor(time, exclude = c(NA, NaN))
  values <- data.frame(time=time, status=status, group=group)

  for(i in 1:k) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }

  Sigma_hat <- diag(var)

  # calculate the test statistics
  for(i in 1:L){
    vec <- hyp_mat[[i]] %*% rmst - hyp_vec[[i]]
    teststats[i] <- t(vec) %*% MASS::ginv(hyp_mat[[i]] %*% Sigma_hat %*% t(hyp_mat[[i]])) %*% vec
    if(teststats[i] == 0) teststats[i] <- ifelse(all(vec == 0), 0, Inf)
  }

  namen <- names(hyp_mat)
  if(is.null(namen)) namen <- 1:length(hyp_mat)

  erg <- replicate(Nres , expr = {

    # generate the groupwiese bootstrap sample

    # calculate the rmst and variance estimator for the bootstrap sample
    bs_rmst <- numeric(k)
    bs_var <- numeric(k)

    for(i in 1:k) {
      # generate a bootstrap sample in group i
      values2 <- values[group == i,][sample(1:n_vec[i],replace = TRUE),]
      #values2 <- sort_data(values2)
      #calculate for each group
      temp <- RMST(values2$time, values2$status, tau)
      bs_rmst[i] <- temp["rmst"]
      bs_var[i] <- temp["var_rmst"]
    }

    # calculate the groupwise BS test statistics for all matrices in hyp_mat
    ts_groupwiseBS <- sapply(hyp_mat, function(c_mat){
      vec <- c_mat %*% (bs_rmst- rmst)
      tbs <- t(vec) %*% MASS::ginv(c_mat %*% diag(bs_var) %*% t(c_mat)) %*% vec
      if(tbs == 0) tbs <- ifelse(all(vec == 0), 0, Inf)
      tbs
    })

    ts_groupwiseBS
  }, simplify = TRUE)
  dim(erg) <- c(L,Nres)

  # (adjusted) p values
  pvalues <- p.value(erg, teststat = teststats)

  # if all matrices are vectors
  if(all(sapply(hyp_mat, function(x) nrow(x)) == 1)){
    res <- matrix(list() , nrow = L, ncol = 8)
    rownames(res) <- names(hyp_mat)
    colnames(res) <- c("hyp_matrix", "estimator", "lwr_conf", "upr_conf", "test_stat", "critical value", "adj_pvalue", "decision")

    quant <- crit_values2(erg, alpha = alpha)

    for(l in 1:L){
      # determine simultaneous 1-alpha confidence intervals
      conf.int <- NULL

      conf.int <- c(hyp_mat[[l]]%*%rmst) + sqrt(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]])) * quant[l]) * c(-1,1)
      # determine the smallest and largest possible value
      conf.int[1] <- max(conf.int[1], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau))
      conf.int[2] <- min(conf.int[2], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0))
      attr(conf.int, "conf.level") <- 1-alpha

      # save the needed values for res
      res[l,3:4] <- conf.int
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
      res[[l,5]] <- teststats[l]
      res[[l,6]] <- quant[l]
    }
    if(stepwise  & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMST.groupwise.test(time, status, group, hyp_mat = hyp_mat[ind],
                                hyp_vec = hyp_vec[ind], tau = tau,
                                stepwise = FALSE, alpha = alpha,
                                Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the confidence interval and critical value
      res <- res[,-c(3,4,6)]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")


  }else{
    res <- matrix(list() , nrow = L, ncol = 6)
    colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")
    rownames(res) <- names(hyp_mat)

    quant <- crit_values2(erg, alpha = alpha)

    for(l in 1:L){

      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
    }

    res[,3] <- teststats
    res[,4] <- quant

    if(stepwise  & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMST.groupwise.test(time, status, group, hyp_mat = hyp_mat[ind],
                                hyp_vec = hyp_vec[ind], tau = tau,
                                stepwise = FALSE, alpha = alpha,
                                Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the critical value
      res <- res[,-4]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

  }

  if(stepwise) method <- "Multiple groupwise RMST Wald-type tests with stepwise extension"
  if(!stepwise) method <- "Multiple groupwise RMST Wald-type tests"

  out <- list(method = method,
              test_stat = teststats,
              p.value = pvalues,
              res = res,
              alpha = alpha)
  class(out) <- "GFDrmst"
  return(out)
}


RMST.permutation.test <- function(time = NULL, status = NULL, group = NULL,
                                  formula = NULL, event = NULL, data = NULL,
                                  hyp_mat, hyp_vec = NULL,
                                  tau, stepwise = FALSE, alpha = 0.05,
                                  Nres = 4999, seed = 1){
  set.seed(seed)

  if(is.null(time) | is.null(status) | is.null(group)){
    dat2 <- formula2input(formula, event, data)
    time   <- dat2$time
    status <- dat2$event
    group  <- dat2$group
  }

  k <- length(unique(group))
  group <- factor(group, labels = 1:k)
  # some possible values for hyp_mat
  if(is.character(hyp_mat)){
    if(hyp_mat %in% c("center", "Dunnett", "Tukey")){
      Y <- GFDmcv::contr_mat(k, type = hyp_mat)
      X <- data.frame(t(Y))
      colnames(X) <- rownames(Y)
      hyp_mat <- as.list(X)
      hyp_mat <- lapply(hyp_mat, t)
    }else{
      if(hyp_mat== "crossed factorial"){
        if(is.null(formula)) stop("Formula needed to realize the crossed factorial design.")
        formula2 <- paste0(formula, "*", event)
        dat <- model.frame(formula2, data)
        nadat <- names(dat)
        nf <- ncol(dat) - 2
        fl <- NA
        for (aa in 1:nf) {
          fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
        }
        levels <- list()
        for (jj in 1:nf) {
          levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
        }
        lev_names <- expand.grid(levels)
        if (nf == 1) {
          hyp_mat <- list(diag(fl) - matrix(1/fl, ncol = fl,
                                            nrow = fl))
        }else {
          lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),]
          formula3 <- as.formula(formula)
          nr_hypo <- attr(terms(formula3), "factors")
          fac_names_original <- colnames(nr_hypo)
          formula1 <- as.formula( gsub("\\+", "*", formula))
          nr_hypo <- attr(terms(formula1), "factors")
          fac_names <- colnames(nr_hypo)
          perm_names <- t(attr(terms(formula1), "factors")[-1, ])

          hyp_mat <- HC(fl, perm_names, fac_names)[[1]]
          names(hyp_mat) <- HC(fl, perm_names, fac_names)[[2]]
          hyp_mat <- hyp_mat[fac_names_original]

        }
      }else{ stop("hyp_mat is an invalid character input") }
    }}else{if(is.null(hyp_vec)){
      warning("hyp_vec is chosen as zero vector since it is unspecified!")
    }}
  if(is.null(hyp_vec)){
    hyp_vec <- lapply(hyp_mat, function(x) rep(0, nrow(x)))
    #warning("hyp_vec is chosen as zero vector since it is unspecified!")
  }
  if(is.matrix(hyp_mat)) hyp_mat <- list(hyp_mat)
  if(!is.list(hyp_vec)) hyp_vec <- list(hyp_vec)
  L <- length(hyp_mat)
  n_vec <- table(group)

  # check whether all hypotheses have a solution in [0,tau]^k
  solution <- logical(length(hyp_mat))
  for(i in 1:length(hyp_mat)){
    solution[i] <- existence(hyp_mat[[i]], hyp_vec[[i]], tau)
  }
  if(any(!solution)){
    stop(paste0("Hypothesis ", which(!solution),
                " does not have a possible solution in [0,tau]^k. "))
  }

  rmst <- numeric(k)
  var <- numeric(k)
  teststats <- numeric(L)

  attr(teststats, "names") <- paste0("W_n(H_", 1:L, ", c_", 1:L, ")")

  status <- factor(status, levels = c(0,1))
  #time <- factor(time, exclude = c(NA, NaN))
  values <- data.frame(time=time, status=status, group=group)

  for(i in 1:k) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }

  Sigma_hat <- diag(var)

  # calculate the test statistics
  for(i in 1:L){
    vec <- hyp_mat[[i]] %*% rmst - hyp_vec[[i]]
    teststats[i] <- t(vec) %*% MASS::ginv(hyp_mat[[i]] %*% Sigma_hat %*% t(hyp_mat[[i]])) %*% vec
    if(teststats[i] == 0) teststats[i] <- ifelse(all(vec == 0), 0, Inf)
  }

  namen <- names(hyp_mat)
  if(is.null(namen)) namen <- 1:length(hyp_mat)

  perm_values <- values[,c("time", "status")]

  erg <- replicate(Nres , expr = {

    # generate the permuted sample
    perm_values$group <- sample(group)

    # calculate the perm. test statistics
    rmst_perm <- numeric(k)
    var_perm <- numeric(k)
    teststats_perm <- numeric(L)

    for(i in 1:k) {
      values2 <- perm_values[perm_values$group == i,]
      #calculate for each group
      temp <- RMST(values2$time, values2$status, tau)
      rmst_perm[i] <- temp["rmst"]
      var_perm[i] <- temp["var_rmst"]
    }

    Sigma_hat_perm <- diag(var_perm)

    # calculate the test statistics
    for(i in 1:L){
      vec <- hyp_mat[[i]] %*% rmst_perm
      teststats_perm[i] <- t(vec) %*% MASS::ginv(hyp_mat[[i]] %*% Sigma_hat_perm %*% t(hyp_mat[[i]])) %*% vec
      if(teststats_perm[i] == 0) teststats_perm[i] <- ifelse(all(vec == 0), 0, Inf)
    }

    teststats_perm
  }, simplify = TRUE)
  dim(erg) <- c(L,Nres)

  # unadjusted p values
  pvalues <- rowMeans(erg > teststats)
  # Bonferroni adjustment
  if(!stepwise) pvalues <- p.adjust(pvalues, method = "bonferroni")
  # Holm adjustment
  if(stepwise) pvalues <- p.adjust(pvalues, method = "holm")

  # if all matrices are vectors
  if(all(sapply(hyp_mat, function(x) nrow(x)) == 1)){
    res <- matrix(list() , nrow = L, ncol = 8)
    colnames(res) <- c("hyp_matrix", "estimator", "lwr_conf", "upr_conf",
                       "test_stat", "critical value", "adj_pvalue", "decision")
    rownames(res) <- names(hyp_mat)

    quant <- numeric(L)
    for(l in 1:L){
      quant[l] <- crit_values2(t(erg[l,]), alpha = alpha/L)
    }

    for(l in 1:L){
      # determine simultaneous 1-alpha confidence intervals
      conf.int <- NULL

      conf.int <- c(hyp_mat[[l]]%*%rmst) + sqrt(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]])) * quant[l]) * c(-1,1)
      # determine the smallest and largest possible value
      conf.int[1] <- max(conf.int[1], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau))
      conf.int[2] <- min(conf.int[2], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0))
      attr(conf.int, "conf.level") <- 1-alpha

      # save the needed values for res
      res[l,3:4] <- conf.int
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
      res[[l,5]] <- teststats[l]
      res[[l,6]] <- quant[l]
    }
    if(stepwise & (L > 1)){# delete the confidence interval and critical value
      res <- res[,-c(3,4,6)]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")


  }else{
    res <- matrix(list() , nrow = L, ncol = 6)
    colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")
    rownames(res) <- names(hyp_mat)

    quant <- numeric(L)
    for(l in 1:L){
      quant[l] <- crit_values2(t(erg[l,]), alpha = alpha/L)
    }

    for(l in 1:L){

      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmst
    }

    res[,3] <- teststats
    res[,4] <- quant

    if(stepwise & (L > 1)){
      # delete the critical value
      res <- res[,-4]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

  }

  if(stepwise) method <- "Permutation RMST Wald-type tests with Holm correction"
  if(!stepwise) method <- "Permutation RMST Wald-type tests with Bonferroni correction"

  out <- list(method = method,
              test_stat = teststats,
              p.value = pvalues,
              res = res,
              alpha = alpha)
  class(out) <- "GFDrmst"
  return(out)
}

# one function for all methods
RMST.test <- function(
    time = NULL,
    status = NULL,
    group = NULL,
    formula = NULL,
    event = NULL,
    data = NULL,
    hyp_mat,
    hyp_vec = NULL,
    tau,
    method = c("groupwise", "permutation", "asymptotic"),
    stepwise = FALSE,
    alpha = 0.05,
    Nres = 4999,
    seed = 1
){
  method <- match.arg(method)
  switch(method,
         groupwise   = RMST.groupwise.test(time, status, group, formula, event,
                                           data, hyp_mat, hyp_vec, tau, stepwise,
                                           alpha, Nres, seed),
         permutation = RMST.permutation.test(time, status, group, formula, event,
                                             data, hyp_mat, hyp_vec, tau, stepwise,
                                             alpha, Nres, seed),
         asymptotic  = RMST.asymptotic.test(time, status, group, formula, event,
                                            data, hyp_mat, hyp_vec, tau, stepwise,
                                            alpha, Nres, seed))
}
