#' Robust EM algorithm for Spline Regression Mixture Model (SRM) with Irregular observations
#' @param X observation times, long form [n*m, 1], i.e. n subjects each with m observations, observation time can be irregular
#' @param Y observed outcomes, long form [n*m, 1]
#' @param id subject indicators to tell which data belongs to which subject [n*m, 1]
#' @param splineOrder B-Spline order
#' @param splineType select between c("spline", "B-spline") indicating either natural cubic spline or b-splines
#' @param nknots number of internal knots selected equally
#' @param verbose Screen output?
#' @param MaxIter Maximum iteration
#' @import splines MASS
#' @examples
#' test = generateToyDataset(n1 = 40, n2 = 30, n3 = 30)
#' X  = rep(test$x, 100)
#' id = rep(seq(100), each = 100)
#' Y  = matrix(t(test$Y), ncol = 1)
#' output = emSRMir(X, Y, id, 3, splineType = "B-spline", 3)
#' @export
emSRMir <- function(X, Y, id, splineOrder, splineType = c("spline", "B-spline"), nknots, verbose = TRUE, MaxIter = 1000) {

  splineType <- match.arg(splineType)

  # n <- nrow(Y)
  # m <- ncol(Y)
  n = length(unique(id))

  M <- splineOrder
  dimBeta <- M + nknots

  # Construction of the design matrix
  if (splineType == "B-spline") {
    Xstack = bs( X, intercept = T, df = nknots+splineOrder, degree = splineOrder-1)
  } else if (splineType == "spline") {
    Xstack = ns( X, intercept = T, df = nknots+2+2)
  }

  x.list = split(as.data.frame(Xstack),  id)
  x.list = lapply(x.list, data.matrix)
  Y.dat  = apply(data.matrix(Y), 2, scale) # in case input Y is a vector
  y.list = split(as.data.frame(Y.dat),  id)
  y.list = lapply(y.list, data.matrix)

  # Threshold for testing the convergence
  Epsilon <- 1e-6

  # ------ Step 1 initialization ----- #
  Beta <- 1
  K <- n
  # -----------(28) ---------------- #
  gama <- 1e-6
  # Ystack <- t(matrix(data = t(as.numeric(repmat(t(Y), n, 1))), nrow = m)) # [nxn x m];
  # dij <- apply( (Ystack - repmat(Y, n, 1)) ^ 2, 1, sum)
  # dmin <- min(dij[dij > 0])
  XX.wait = crossprod(Xstack)
  Xy.wait = crossprod(Xstack, t(t(Y.dat)))
  Y2.wait = colSums(Y.dat^2)
  SSR.all = Y2.wait - diag(t(Xy.wait) %*% ginv(XX.wait) %*% Xy.wait)
  e.sigma = (SSR.all)/(nrow(Xstack) - dimBeta)
  Q <- e.sigma

  # Ytild <- matrix(data = t(Y), ncol = 1) # long form now
  Ytild <- Y.dat

  # Initialize the mixing proportions
  Alphak <- 1 / K * ones(K, 1)
  Pik <- 1 / K * ones(K, 1)

  # Initialize the regression parameters and the variances
  Betak <-  zeros(dimBeta, n)
  Sigmak2 <- zeros(K, 1)

  for (k in 1:K) {
    x.in    = x.list[[k]]
    y.in    = y.list[[k]]
    N.in    = dim(x.in)[k]
    XX      = crossprod(x.in) # faster than t(x.in) %*% x.in
    A.mat   = ginv(XX)
    Xy      = crossprod(x.in, t(t(y.in))) # t(x.in) %*% t(t(y.in))
    Y2      = colSums(y.in^2)
    SSR0    = diag(Y2 - t(Xy) %*% A.mat %*% Xy)
    betak   = crossprod(A.mat, Xy)
    Betak[, k] <- betak
    Dk      = unlist(lapply(split((Y.dat - Xstack %*% betak)^2, id), sum))
    Dk      <- sort(Dk)
    Sigmak2[k] <- Dk[ceiling(sqrt(K))]

    # 1/m added/removed recently; or Dk(end)
    # --------------------------- #

  }


  #--------- Step 3 (4) --------#

  # Compute the posterior cluster probabilites (responsibilities) for the
  # initial guess of the model parameters
  #############################
  #                           #
  #       E-Step              #
  #                           #
  #############################
  PikFik <- zeros(n, K)
  log_fk_xij <- zeros(nrow(Xstack), K)
  log_Pik_fk_Xi <- zeros(n, K)
  log_Pik_Fik <- zeros(n, K)

  # E-Step
  for (k in 1:K) {

    pik <- Pik[k]
    betak <- Betak[, k]
    sigmak2 <- Sigmak2[k]
    #fik = normpdf(X,muk,sigmak); #Gaussian density
    z <- ((Ytild - Xstack %*% betak) ^ 2) / sigmak2
    log_fk_xij[, k] <- -0.5 * (log(2 * pi) + log(sigmak2)) - 0.5 * z  # [nxm x 1] : univariate Gaussians
    # log-lik for the expected n_k curves of cluster k
    log_fk_Xi <- unlist(lapply(split(log_fk_xij[,k], id), sum)) #  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
    log_Pik_fk_Xi[, k] <- log(pik) + log_fk_Xi # [n x K]

    log_Pik_Fik[, k] <- log_Pik_fk_Xi[, k]
    #PikFik(:,k) = pik * exp(log_fk_Xi);
  }
  # Posterior = PikFik./(sum(PikFik,2)*ones(1,K));
  log_Prosterior <- lognormalize(log_Pik_fk_Xi)
  Posterior <- exp(log_Prosterior)
  Tauik <- Posterior


  #############################
  #                           #
  # main Robust EM-MxReg loop #
  #                           #
  #############################

  stored_J <- c() # to store the maximized penalized log-likelihood criterion
  pen_loglik_old <- -Inf
  iter <- 1 # iteration number
  converged <- FALSE
  stored_K <- c() # To store the estimatde number of clusters at each iteration

  while (!converged && (iter <= MaxIter)) {

    stored_K <- c(stored_K, K)

    # Print the value of the optimized criterion
    # pen_loglik = sum(log(sum(PikFik,2)),1 ) + Beta*n*sum(Alphak.*log(Alphak));
    # pen_loglik = (sum(log_Prosterior(:) .* Posterior(:)) - sum(log_Prosterior(:) .* log_Pik_Fik(:)))+ Beta*n*sum(Alphak.*log(Alphak));
    pen_loglik <- sum(logsumexp(log_Pik_Fik, 1)) + Beta * n * sum(Alphak * log(Alphak))

    if (verbose) {
      message("EM Iteration: ", iter - 1, " | Number of clusters K: ", K, " | Penalized log-likelihood: "  , pen_loglik)
    }


    #############################
    #                           #
    #       M-Step              #
    #                           #
    #############################
    for (k in 1:K) {

      tauik <- Tauik[, k]
      # ------Step 4 (25) ---------- #
      # Update of the regression coefficients
      Wk <- matrix(rep(tauik, table(id)), ncol = 1) # cluster_weights(:)% [mn x 1]
      # meme chose
      # temp =  repmat(tauik,1,m)';% [m x n]

      # cluster_weights = cluster_weights(:);
      wYk <- sqrt(Wk) * Ytild # fuzzy cluster k
      wXk <- sqrt(Wk %*% ones(1, dimBeta)) * Xstack # [(n*m)*(M+nknots)]x
      # maximization w.r.t betak: Weighted least squares
      # betak  = inv(phik'*phik + 1e-4*eye(dimBeta))*phik'*Yk;
      betak <- solve(crossprod(wXk, wXk) + 1e-4 * diag(dimBeta), crossprod(wXk, wYk))
      Betak[, k] <- betak

      # ------ Cooected with step 5 (13) ---------- #
      # mixing proportions : alphak_EM
      pik <- sum(tauik) / n # alpha_k^EM
      Pik[k] <- pik
    }
    # ------- step 5 (13) ------- #
    AlphakOld <- Alphak

    Alphak <- Pik + Beta * Alphak * (log(Alphak) - sum(Alphak * log(Alphak)))

    # ------- step 6 (24)  ------- #
    # update beta
    E <- sum(AlphakOld * log(AlphakOld))
    # eta <- min(1, 0.5 ^ floor((median(table(id)) / 2) - 1))
    eta <- min(1, 0.5 ^ floor((max(table(id)) / 2) - 1))
    pik_max <- max(Pik)
    alphak_max <- max(AlphakOld)
    Beta <- min(sum(exp(-eta * n * abs(Alphak - AlphakOld))) / K, (1 - pik_max) / (-alphak_max * E))

    # ------- step 7 --------- #
    # Kold = K;
    # update the number of clusters K
    small_klas <- which(Alphak < 1 / n)
    # ------- step 7  (14) ------- #
    K <- K - length(small_klas)

    # discard the small clusters
    if (length(small_klas) > 0) {
      Pik <- Pik[-small_klas]
      Alphak <- Alphak[-small_klas]
      log_fk_xij <- log_fk_xij[, -small_klas]
      log_Pik_fk_Xi <- log_Pik_fk_Xi[, -small_klas]
      log_Pik_Fik <- log_Pik_Fik[, -small_klas]
      PikFik <- PikFik[, -small_klas]
      log_Prosterior <- log_Prosterior[, -small_klas]
      Posterior <- Posterior[, -small_klas]
      Sigmak2 <- Sigmak2[-small_klas]
      Betak <- Betak[, -small_klas]
    }
    # ------- step 7  (15) normalize the Pik and Alphak ------- #
    Pik <- Pik / sum(Pik)
    Alphak <- Alphak / sum(Alphak)
    # ------- step 7 (16)  normalize the posterior prob ------- #
    Posterior <- Posterior / (apply(X = Posterior, MARGIN = 1, sum) %*% ones(1, K))
    Tauik <- Posterior

    # -------- step 7 ------------ #
    # test if the partition is stable (K not changing)
    nit <- 60
    if ((iter >= nit) && (stored_K[iter - (nit - 1)] - K) == 0) {
      Beta <- 0
    }
    # -----------step 8 (26) and (28) ---------------- #
    #############################
    #                           #
    #       M-Step              #
    #                           #
    #############################
    for (k in 1:K) {
      tauik <- Tauik[, k]

      # temp <- repmat(tauik, 1, m)
      # Wk = matrix(data = t(temp), ncol = 1)
      Wk <- matrix(rep(tauik, table(id)), ncol = 1)
      wYk <- sqrt(Wk) * Ytild
      wXk <- sqrt(Wk %*% ones(1, dimBeta)) * Xstack

      betak <- Betak[, k]

      # ----------- (26) ---------------- #
      # update the variance
      sigmak2 <- sum((wYk - wXk %*% betak) ^ 2) / sum(Wk)
      # -----------(28) ---------------- #
      sigmak2 <- (1 - gama) * sigmak2 + gama * Q
      Sigmak2[k] <- sigmak2
    }

    # -----------step 9 (4) ---------------- #
    #############################
    #                           #
    #       E-Step              #
    #                           #
    #############################
    for (k in 1:K) {
      alphak <- Alphak[k]
      betak <- Betak[, k]
      sigmak2 <- Sigmak2[k]

      # fik = normpdf(X,muk,sigmak); %Gaussian density
      z <- ((Ytild - Xstack %*% betak) ^ 2) / sigmak2
      log_fk_xij[, k] <- -0.5 * (log(2 * pi) + log(sigmak2)) - 0.5 * z # [nxm x 1]
      # log-lik for the n_k curves of cluster k
      log_fk_Xi <- unlist(lapply(split(log_fk_xij[,k], id), sum))
      log_Pik_fk_Xi[, k] <- log(alphak) + log_fk_Xi # [nxK]

      log_Pik_Fik[, k] <- log_Pik_fk_Xi[, k]
      # PikFik(:,k) = pik * exp(log_fk_Xi);
    }
    # PikFik = exp(log_Pik_Fik);
    # Posterior = PikFik./(sum(PikFik,2)*ones(1,K));
    # Posterior = exp(log_normalize(log_Pik_fk_Xi));
    log_Posterior <- lognormalize(log_Pik_fk_Xi)
    Posterior <- exp(lognormalize(log_Posterior))
    Tauik <- Posterior

    ##########
    # compute the value of the optimized criterion J (12) #
    #pen_loglik = sum(log(sum(PikFik,2)),1 ) + Beta*n*sum(Alphak.*log(Alphak));
    #pen_loglik = sum(logsumexp(log_Pik_Fik,2),1) + Beta*n*sum(Alphak.*log(Alphak));
    #pen_loglik = (sum(log_Prosterior(:) .* Posterior(:)) - sum(log_Prosterior(:) .* log_Pik_Fik(:)))+ Beta*n*sum(Alphak.*log(Alphak));
    stored_J <- c(stored_J, pen_loglik)
    #     fprintf(1,'EM Iteration : #d  | number of clusters K : #d | penalized loglikelihood: #f \n',iter, K, pen_loglik);
    #########

    # -----------step 10 (25) ---------------- #
    #############################
    #                           #
    #       M-Step              #
    #                           #
    #############################
    BetakOld <- Betak
    for (k in 1:K) {
      tauik <- Tauik[, k]
      # pik(k) = sum(tauik)/n;

      # update of the regression coefficients
      # temp <- repmat(tauik, 1, m)
      # Wk <- matrix(data = t(temp), ncol = 1) # cluster_weights
      Wk <- matrix(rep(tauik, table(id)), ncol = 1)
      wYk <- sqrt(Wk) * Ytild
      wXk <- sqrt(Wk %*% ones(1, dimBeta)) * Xstack
      # maximization w.r.t betak: Weighted least squares
      # betak  = inv(phik'*phik + 0.0001*eye(dimBeta))*phik'*Yk;
      betak <- solve(crossprod(wXk, wXk), crossprod(wXk, wYk))

      Betak[, k] <- betak
    }

    # -----------step 11 ---------------- #
    # test of convergence

    distBetak <- sqrt(apply((Betak - BetakOld) ^ 2, MARGIN = 1, FUN = sum))

    converged <- (max(distBetak) < Epsilon || abs((pen_loglik - pen_loglik_old) / pen_loglik_old) < Epsilon)
    if (is.na(converged)) {
      converged <- FALSE
    }

    pen_loglik_old <- pen_loglik

    iter <- iter + 1

  } # end of the Robust EM loop

  klas <- apply(Posterior, MARGIN = 1, which.max)

  gmm_density <- apply(PikFik, MARGIN = 1, sum)

  params <- list()
  params$label <- cbind(unique(id), klas)
  params$Pik <- Pik
  params$Alphak <- Alphak
  params$Betak <- Betak
  params$Muk <- Xstack %*% Betak
  params$Sigmak2 <- Sigmak2

  return(params)

}
