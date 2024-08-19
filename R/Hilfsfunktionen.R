

# is there a solution mu for the hypothesis H %*% mu = c in [0,tau]^k
existence <- function(H, c, tau){
  # set c >= 0
  H[c < 0,] <- -H[c < 0,]
  c <- abs(c)

  r <- nrow(H)
  k <- ncol(H)

  Objective.in <- c(rep(0,k),rep(1,r))
  zero_mat_r <- matrix(0, nrow=r, ncol=k)
  zero_mat_k <- matrix(0, nrow=k, ncol=r)
  Const.mat <- rbind(cbind(H,diag(r)),
                     cbind(diag(k),zero_mat_k),
                     cbind(diag(k),zero_mat_k),
                     cbind(zero_mat_r,diag(r)))
  Const.rhs <- c(c, numeric(k), rep(tau,k), numeric(r))
  Const.dir <- c(rep("=",r),rep(">=",k),rep("<=",k),rep(">=",r))

  Optimum <- lp(direction="min",Objective.in,Const.mat,Const.dir,Const.rhs)
  Optimum$objval == 0
}

#No main effect of factor a
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_A <- function(n_group_a, n_group_b) {
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  jbb <- matrix(1,1, n_group_b)/n_group_b
  ha <- pa %x% jbb

  return(ha)
}

#No main effect of factor B
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_B <- function(n_group_a, n_group_b) {
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  jaa <- matrix(1, 1, n_group_a)/n_group_a
  ha <- jaa %x% pb

  return(ha)
}


#null_mat_fac(2, 2)

#No interaction effect
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_AB <- function(n_group_a, n_group_b){
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  hab <- pa %x% pb

  return(hab)
}



# if a formula is given, we need to extract the time, status, group vector from the formula + data
formula2input <- function(formula, event = "event", data){
  formula2 <- paste0(formula, "*", event)
  dat <- model.frame(formula2, data)
  subject <- 1:nrow(dat)
  n_all <- length(subject)
  formula <- as.formula(formula)
  nf <- ncol(dat) - 2
  nadat <- names(dat)
  if (anyNA(data[, nadat])) {
    stop("Data contains NAs!")
  }
  names(dat) <- c("time", nadat[2:(1 + nf)], "event")
  dat2 <- data.frame(dat, subject = subject)
  nadat2 <- nadat[-c(1, nf + 2)]
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
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    dat2$group <- rep(1:length(n), n)
  }else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),
    ]
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    dat2$group <- rep(1:length(n), n)
  }
  return(dat2[,c("time", "event", "group")])
}

# show which groups correspond to which factor values
formula2groups <- function(formula, event = "event", data){
  formula2 <- paste0(formula, "*", event)
  dat <- model.frame(formula2, data)
  nf <- ncol(dat) - 2
  nadat <- names(dat)
  if (anyNA(data[, nadat])) {
    stop("Data contains NAs!")
  }
  names(dat) <- c("time", nadat[2:(1 + nf)], "event")
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)
  if(nf > 1) lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),]
  lev_names <- cbind(1:nrow(lev_names), lev_names)
  colnames(lev_names) <- c("Group", colnames(dat)[2:(nf + 1)])
  return(lev_names)
}

#Sort object of data_gen
#Input:
# values:   matrix - created by data_gen

sort_data <- function(values) values[order(values[, 1]), ]

# ## a faster function than table(), which is not really faster...
# # x: numeric vector of survival times
# # y: numeric vector of 0 (= censored) and 1 (= uncensored)
# # output
# # table of x, y
# fast_table <- function (x,y){
#   x_uniq <- unique(x)
#   out <- matrix(NA, ncol=2, nrow=length(x_uniq))
#   colnames(out) <- c(0,1)
#   rownames(out) <- x_uniq
#
#   censored <- (y==0)
#   x_cens <- x[censored]
#   x_uncens <- x[!censored]
#   out[,1] <- sapply(x_uniq, function(z) sum(x_cens==z))
#   out[,2] <- sapply(x_uniq, function(z) sum(x_uncens==z))
#   out <- out[order(x_uniq),]
#
#   out
# }

# function: simple_surv() determines basic survival quantities in a faster
#           manner than survfit()
# input: time - survival time
#        cens - indicator for censoring (1=not censored, 0=censored)
# output: matrix of dimensions (n, 4) to speed up later calculations

simple_surv <- function(time, cens) {
  n.all <- length(time)
  tab <- table(time, cens)
  d <- tab[,1] # number of withdrawals in time points
  w <- tab[,2] # number of events in time points
  n <- c(n.all, n.all - cumsum(d + w)) # calculate risk set
  n <- n[-length(n)]
  s <- cumprod(1 - (w /  n)) # calculate Kaplan-Meier-Estimator
  cbind(as.numeric(row.names(tab)), s, w, n)
}

# function: RMST() estimates Restricted Mean Survival Time and its variance
#           for a sample of tupel (time, cens), given a time point tau
# input: time - survival time
#        cens - indicator for censoring (1=censored, 0=not censored)
#        tau - restriction time point
# output: rmst - estimated restricted mean survival time for tau
#         var_rmst - estimated variance of the rmst

RMST <- function(time, cens, tau){
  #n <- length(time) # number of observation in the beginning of study
  survtab <- simple_surv(time, cens) # fit convenient model for quantities

  t_i <- survtab[,1] <= tau # identify time points t_i <= tau
  t_sel <- survtab[,1][t_i] # select relevent time points
  S_km <- survtab[,2][t_i] # calculate Kaplan-Meier estimator

  w.factor <- diff(c(0, t_sel,tau)) # width of the area under S(t) from t_0=0
  rmst <- sum(w.factor * c(1, S_km)) # calculate area under S(t) up to t_i=tau

  area <- rev(cumsum(rev(diff(c(t_sel,tau)) * S_km)))^2 # calculate areas under S(t) from t_i to tau for
  d <- survtab[,3][t_i] # determine number of events
  Y <- survtab[,4][t_i] # determine number of individuals under risk
  var_rmst <- sum(( d/Y/(Y-d)) * area, na.rm = TRUE) # calculate final variance without
  # the factor n_i, because this factor
  # is eliminated by the factor in front
  # of the variances (i.e. n/n_i) and in
  # front of the test statistic (i.e. n^{1/2})
  # We add na.rm = TRUE to the sum to easily handle
  # the case Y = d (in the case the last summand is 0/0)
  c(rmst = rmst, var_rmst = var_rmst)
}


# calculate the p value
p.value <- function(data_mat, teststat){
  data_order <- t( apply(data_mat, 1, sort) )
  pvalues <- numeric(nrow(data_mat))
  B <- ncol(data_mat)

  for(l in 1:length(teststat)){
    if(teststat[l] < data_order[l,1]){
      pvalues[l] <- 1
    }else{
      beta <- mean(data_mat[l,] >= teststat[l])
      if(beta < 1){
        x <- round((1-beta)*B)
        quants <- data_order[,x]
      }else{
        quants <- -Inf
      }
      pvalues[l] <- mean(apply(data_mat > quants, 2, any))
    }
  }
  pvalues
}


# half partition search, quicker for higher degree of dependence
crit_values2 <- function(data_mat, alpha = 0.05){
  # the input is a matrix data_mat of dimenson n_obs x dim_obs containing the B (resampling) observations X_b of dimension dim_obs
  # for our specific purpose X_i is the vector of the test statistics based on the different contrast vectors, i.e. dim_obs = r in the pdf
  # each column represent one observation (i.e. one resampling iteration)
  # the output are the critical values for each dimension (large X_b lead to a rejection)
  # such that the overall level is alpha
  n <- length(data_mat[1,])
  dimension <- length(data_mat[,1])

  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t( apply(data_mat, 1, sort) )

  # worst case1
  #for each point outside the box only one coordinate leads to a rejection.
  # Thus, for each coordinate (alpha/dim) * number of obs are outside the box  (Bonferoni correction)
  j_low <- ceiling( (alpha/dimension) * n ) - 1
  A <- data_mat / data_order[,n - j_low]
  A[is.na(A)] <- 1
  alpha_low <- mean( apply(A, 2, max) > 1 ) # count the points being outside the box (where the box borders are given by the critical values)

  # worst case1
  # something like totally linear dependence (in the dimension two)
  j_high <- ceiling( alpha * n )
  A <- data_mat / data_order[,n - j_low]
  A[is.na(A)] <- 1
  alpha_high <- mean( apply(A, 2, max) > 1 ) # count the points being outside the box (where the box borders are given by the critical values)

  # now we search for values j_low and j_high = j_low + 1 , such that alpha_low <= alpha and alpha_high > alpha
  while( j_high - j_low > 1){
    j_mid <- ceiling( j_low + (j_high - j_low)/2) # approx. middle between j_low and j_high
    A <- data_mat / data_order[,n - j_mid]
    A[is.na(A)] <- 1
    alpha_sim <- mean( apply(A, 2, max) > 1 )
    ifelse(alpha_sim <= alpha, j_low <- j_mid, j_high <- j_mid)
  }
  data_order[,n - j_low ] # critical values
}

# function for getting the global hypothesis matrix from a list of single hypothesis matrices
global_mat <- function(c_mat, k) matrix(unlist(c_mat), ncol = k, byrow=TRUE)

# This function is adapted from the HC function in the GFDsurv package
# https://CRAN.R-project.org/package=GFDsurv
HC <- function (fl, perm_names, names)
{
  nf <- length(fl)
  P <- function(x) {
    P <- diag(x) - matrix(1/x, ncol = x, nrow = x)
    return(P)
  }
  One <- function(x) {
    I <- t(1/x * matrix(1, ncol = 1, nrow = x))
    return(I)
  }
  tmp <- 0
  for (i in 1:nf) {
    tmp <- c(tmp, choose(nf, i))
    nh <- sum(tmp)
  }
  Z <- 0:(nf - 1)
  position <- rep(0, nh)
  for (i in 1:nh) {
    position[i] <- position[i] + sum(2^Z[which(perm_names[i,
    ] == 1)])
  }
  Vektor <- rep(NA, nh)
  for (j in 1:nh) {
    Vektor[position[j]] <- j
  }
  fac_names <- names[Vektor]
  kp <- function(A) {
    kp <- A[[1]]
    for (i in 2:length(A)) {
      kp <- kp %x% A[[i]]
    }
    return(kp)
  }
  A <- list()
  for (i in 1:nf) {
    A[[i]] <- One(fl[i])
  }
  A_alt <- A
  B <- list()
  for (i in 1:nf) {
    B[[i]] <- P(fl[i])
  }
  n1 <- c(0, 1)
  liste <- vector("list", nf)
  for (i in 1:nf) {
    liste[[i]] <- n1
  }
  G <- expand.grid(liste)
  G <- G[2:(dim(G)[[1]] - 1), ]
  C <- list()
  for (i in 1:dim(G)[1]) {
    index <- which(G[i, ] == 1)
    for (j in 1:length(index)) {
      A[[index[j]]] <- B[[index[j]]]
    }
    C[[i]] <- A
    A <- A_alt
  }
  hypo <- vector("list", nh)
  for (i in 1:(nh - 1)) {
    C_tmp <- C[[i]]
    hypo[[i]] <- kp(C_tmp)
  }
  hypo[[nh]] <- kp(B)
  return(list(hypo, fac_names))
}
