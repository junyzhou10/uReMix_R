#' @export
emPRM <- function(X, Y, p, nknots, verbose = TRUE) {
  n <- nrow(Y)
  m <- ncol(Y)

  # Construction of the desing matrix
  X <- designmatrix(X, p)$XBeta
  # n regularly sampled curves
  Xstack <- repmat(X, n, 1);#desing matrix [(n*m) x (p+1)]

  # Threshold for testing the convergence
  Epsilon <- 1e-6

  # ------ Step 1 initialization ----- #
  Beta <- 1
  K <- n
  # -----------(28) ---------------- #
  gama <- 1e-6
  Ystack <- t(matrix(data = t(as.numeric(repmat(t(Y), n, 1))), nrow = m)) # [nxn x m];
  dij <- apply(X = (Ystack - repmat(Y, n, 1)) ^ 2, 1, sum)
  dmin <- min(dij[dij > 0])
  Q <- dmin

  Ytild <- matrix(data = t(Y), ncol = 1)

  # Initialize the mixing proportions
  Alphak <- 1 / K * ones(K, 1)
  Pik <- 1 / K * ones(K, 1)

   # Initialize the regression parameters and the variances
  Betak <-  zeros(p+1, n)
  Sigmak2 <- zeros(K, 1)

  for (k in 1:K) {

    # ------- step 2  (27)  ------- #
    # Betak  = inv(Phi'*Phi + 1e-4*eye(dimBeta))*Phi'*Y_in(k,:)';
    betak <- solve(crossprod(X, X), t(X) %*% Y[k, ]) # Inversion problem for spline of order 1 (polynomial degree=1)
    Betak[, k] <- betak
    muk <- X %*% betak
    # Dk = sum((reshape(X,n,m) - reshape(muk',n,m)).^2, 2);
    Dk <- apply(X = (Y - ones(n, 1) %*% t(muk)) ^ 2, MARGIN = 1, sum)
    Dk <- sort(Dk)

    # Sigmak2(k)=  1/m*max(Dk);#Dk(ceil(sqrt(K)));#Dk(ceil(sqrt(K)));#.001;##;#
    Sigmak2[k] <- Dk[ceiling(sqrt(K))] # Dk(ceil(sqrt(K)));#.001;##;#
    # Sigmak2(k) = sum(Y_in(k,:)' - muk);

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
  log_fk_xij <- zeros(n * m, K)
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
    log_fk_Xi <- apply(X = matrix(data = log_fk_xij[, k], m, n), MARGIN = 2, sum) # [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
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
  MaxIter <- 1000
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
      temp <- repmat(tauik, 1, m) # [m x n]
      Wk <- matrix(data = t(temp), ncol = 1) # cluster_weights(:)% [mn x 1]
      # meme chose
      # temp =  repmat(tauik,1,m)';% [m x n]

      # cluster_weights = cluster_weights(:);
      wYk <- sqrt(Wk) * Ytild # fuzzy cluster k
      wXk <- sqrt(Wk %*% ones(1, p+1)) * Xstack # [(n*m)*(M+nknots)]x
      # maximization w.r.t betak: Weighted least squares
      # betak  = inv(phik'*phik + 1e-4*eye(p+1))*phik'*Yk;
      betak <- solve(crossprod(wXk, wXk) + 1e-4 * diag(p+1), crossprod(wXk, wYk))
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
    eta <- min(1, 0.5 ^ floor((m / 2) - 1))
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

      temp <- repmat(tauik, 1, m)
      Wk = matrix(data = t(temp), ncol = 1)
      wYk <- sqrt(Wk) * Ytild
      wXk <- sqrt(Wk %*% ones(1, p+1)) * Xstack

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
      log_fk_Xi <- apply(X = matrix(data = log_fk_xij[, k], nrow = m, ncol = n), MARGIN = 2, sum) # [n x m]:  sum over j=1,...,m: fk_Xi = prod_j sum_k pi_{jk} N(x_{ij},mu_{k},s_{k))
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
      temp <- repmat(tauik, 1, m)
      Wk <- matrix(data = t(temp), ncol = 1) # cluster_weights
      wYk <- sqrt(Wk) * Ytild
      wXk <- sqrt(Wk %*% ones(1, 1+p)) * Xstack
      # maximization w.r.t betak: Weighted least squares
      # betak  = inv(phik'*phik + 0.0001*eye(dimBeta))*phik'*Yk;
      betak <- solve(crossprod(wXk, wXk) + 0.01*diag(p+1), crossprod(wXk, wYk))

      Betak[, k] <- betak
    }

    # -----------step 11 ---------------- #
    # test of convergence

    distBetak <- sqrt(apply(X = (Betak - BetakOld) ^ 2, MARGIN = 1, FUN = sum))

    converged <- (max(distBetak) < Epsilon || abs((pen_loglik - pen_loglik_old) / pen_loglik_old) < Epsilon)
    if (is.na(converged)) {
      converged <- FALSE
    }

    pen_loglik_old <- pen_loglik

    iter <- iter + 1

  } # end of the Robust EM loop

  klas <- apply(X = Posterior, MARGIN = 1, which.max)

  gmm_density <- apply(X = PikFik, MARGIN = 1, sum)

  params <- list()
  params$Pik <- Pik
  params$Alphak <- Alphak
  params$Betak <- Betak
  params$Muk <- X %*% Betak
  params$Sigmak2 <- Sigmak2

  return(params)

}
