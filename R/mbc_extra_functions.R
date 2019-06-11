
draw_future_transactions = function (cal.cbs, draws, T.star = cal.cbs$T.star, sample_size = NULL, ncores = 1, trace = 10000) {
  if (is.null(sample_size)) {
    # nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
    nr_of_draws = sum(sapply(draws$level_2, nrow))
  } else {
    stopifnot(is.numeric(sample_size))
    nr_of_draws <- as.integer(sample_size)
  }
  stopifnot(nr_of_draws >= 1)
  nr_of_cust <- length(draws$level_1)
  # parameters <- varnames(draws$level_1[[1]])
  parameters = c("k", "lambda", "mu", "tau", "z")
  if (nr_of_cust != nrow(cal.cbs))
    stop("mismatch between number of customers in parameters 'cal.cbs' and 'draws'")
  if (is.null(T.star))
    stop("T.star is missing")



  if (length(T.star) == 1)
    T.star <- rep(T.star, nr_of_cust)

  draw_left_truncated_gamma <- function(lower, k, lambda) {
    rand <- runif(1, pgamma(lower, k, k * lambda), 1)
    qgamma(rand, k, k * lambda)
  }

  n_row_cal_cbs = nrow(cal.cbs)

  # x.stars <- array(NA_real_, dim = c(nr_of_draws, nr_of_cust))
  if(ncores == 1) {
    message("Running in 1 core.")
  } else{
    message("Running in parallel in ", ncores, " cores")
  }

  # Ver si se puede mejorar el problema de memoria: https://github.com/tdhock/mclapply-memory
  x.stars = parallel::mclapply(1:n_row_cal_cbs, function(cust){

    x_stars_cust = rep(0.0, nr_of_draws)

    print_flag = cust %% trace == 1
    init_time = Sys.time()
    if(print_flag) cat("Iter", scales::comma(cust), "of", scales::comma(n_row_cal_cbs), "\n")
    Tcal <- cal.cbs$T.cal[cust]
    Tstar <- T.star[cust]
    tx <- cal.cbs$t.x[cust]
    taus <- drop(draws$level_1[[cust]][, "tau"])
    lambdas <- drop(draws$level_1[[cust]][, "lambda"])
    ks <- drop(draws$level_1[[cust]][, "k"])

    stopifnot(length(taus) == length(ks) && length(taus) == length(lambdas))

    if (!is.null(sample_size)) {
      idx <- sample(length(taus), size = sample_size, replace = TRUE)
      taus <- taus[idx]
      ks <- ks[idx]
      lambdas <- lambdas[idx]
    }

    alive <- (taus > Tcal)
    for (draw in which(alive)) {
      itts <- draw_left_truncated_gamma(Tcal - tx, ks[draw], lambdas[draw])
      # Tenía pmax y pmin para escalares. Usar min o max es 1000 veces más rápido.
      minT <- min(Tcal + Tstar - tx, taus[draw] - tx)
      nr_of_itt_draws <- max(10, round(minT * lambdas[draw]))

      itts <- c(itts, rgamma(nr_of_itt_draws * 2,
                             shape = ks[draw],
                             rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rgamma(nr_of_itt_draws * 4,
                               shape = ks[draw],
                               rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rgamma(nr_of_itt_draws * 800,
                               shape = ks[draw],
                               rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT) {
        message("not enough inter-transaction times sampled! cust:",
                cust, " draw:", draw, " ", sum(itts), " < ", minT)
        x_stars_cust[draw] = NA
      } else{
        x_stars_cust[draw] <- sum(cumsum(itts) < minT)
      }

    }

    end_time = Sys.time()
    time_diff = as.numeric(difftime(end_time, init_time, units = "mins"))
    if(print_flag) cat("\tTime difference of", time_diff, "minutes\n\n\n")

    return(x_stars_cust)
  },
  mc.cores = ncores)

  message("Converting to array")
  x.stars = simplify2array(x.stars)

  return(x.stars)
}






pggg_mcmc_DrawParameters = function(cal.cbs, mcmc = 2500, burnin = 500, thin = 50, chains = 2, mc.cores = NULL,
                                    param_init = NULL, trace = 100, log_filename = NULL) {

  if(is.null(log_filename)){
    date_time = gsub(pattern = "[: -]", x = as.character(Sys.time()), replacement = "_", fixed = F)
    log_filename = paste0("log_pggg_draws_", date_time, ".log")
  }
  file.create(log_filename)

  # ** methods to sample heterogeneity parameters {r, alpha, s, beta, t, gamma} **

  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    if (type == "lambda") {
      x <- level_1["lambda", ]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    } else if (type == "mu") {
      x <- level_1["mu", ]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- unlist(hyper_prior[c("s_1", "s_2", "beta_1", "beta_2")])
    } else if (type == "k") {
      x <- level_1["k", ]
      cur_params <- c(level_2["t"], level_2["gamma"])
      hyper <- unlist(hyper_prior[c("t_1", "t_2", "gamma_1", "gamma_2")])
    }
    slice_sample_gamma_parameters(x, cur_params, hyper, steps = 200, w = 0.1)
  }

  # ** methods to sample individual-level parameters **

  draw_k <- function(data, level_1, level_2) {
    pggg_slice_sample("k",
                      x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt,
                      k = level_1["k", ], lambda = level_1["lambda", ],
                      mu = level_1["mu", ], tau = level_1["tau", ],
                      t = level_2["t"], gamma = level_2["gamma"],
                      r = level_2["r"], alpha = level_2["alpha"],
                      s = level_2["s"], beta = level_2["beta"])
  }

  draw_lambda <- function(data, level_1, level_2) {
    pggg_slice_sample("lambda",
                      x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt,
                      k = level_1["k", ], lambda = level_1["lambda", ],
                      mu = level_1["mu", ], tau = level_1["tau", ],
                      t = level_2["t"], gamma = level_2["gamma"],
                      r = level_2["r"], alpha = level_2["alpha"],
                      s = level_2["s"], beta = level_2["beta"])
  }

  draw_mu <- function(data, level_1, level_2) {
    N <- nrow(data)
    tau <- level_1["tau", ]
    s <- level_2["s"]
    beta <- level_2["beta"]

    mu <- rgamma(n = N, shape = s + 1, rate = beta + tau)
    mu[mu == 0 | log(mu) < -30] <- exp(-30)  # avoid numeric overflow
    return(mu)
  }

  draw_tau <- function(data, level_1, level_2) {
    N <- nrow(data)
    x <- data$x
    tx <- data$t.x
    Tcal <- data$T.cal
    lambda <- level_1["lambda", ]
    k <- level_1["k", ]
    mu <- level_1["mu", ]

    # sample z
    p_alive <- pggg_palive(x, tx, Tcal, k, lambda, mu)
    alive <- p_alive > runif(n = N)

    # sample tau
    tau <- numeric(N)

    # Case: still alive - left truncated exponential distribution -> [Tcal, Inf]
    if (any(alive)) {
      tau[alive] <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }

    # Case: churned - distribution of tau truncated to [tx, pmin(tx+1, Tcal)]
    if (any(!alive)) {
      tau[!alive] <- pggg_slice_sample("tau", x = data$x[!alive],
                                       tx = data$t.x[!alive], Tcal = data$T.cal[!alive],
                                       litt = data$litt[!alive],
                                       k = level_1["k", !alive],
                                       lambda = level_1["lambda", !alive],
                                       mu = level_1["mu", !alive],
                                       tau = level_1["tau", !alive],
                                       t = level_2["t"],
                                       gamma = level_2["gamma"],
                                       r = level_2["r"],
                                       alpha = level_2["alpha"],
                                       s = level_2["s"],
                                       beta = level_2["beta"])
    }

    return(tau)
  }



  run_single_chain <- function(chain_id, data, hyper_prior) {

    ## initialize arrays for storing draws ##

    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc - 1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 6))
    dimnames(level_2_draws)[[2]] <- c("t", "gamma", "r", "alpha", "s", "beta")

    cat("Initializing array in chain", chain_id, "...\n")
    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 5, nr_of_cust))
    cat("\tdone initializing array in chain", chain_id, "\n")

    dimnames(level_1_draws)[[2]] <- c("k", "lambda", "mu", "tau", "z")

    ## initialize parameters ##

    level_2 <- level_2_draws[1, ]
    level_2["t"] <- param_init$t
    level_2["gamma"] <- param_init$gamma
    level_2["r"] <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"] <- param_init$s
    level_2["beta"] <- param_init$beta

    level_1 <- level_1_draws[1, , ] # nolint
    level_1["k", ] <- 1
    level_1["lambda", ] <- mean(data$x) / mean(ifelse(data$t.x == 0, data$T.cal, data$t.x))
    level_1["tau", ] <- data$t.x + 0.5 / level_1["lambda", ]
    level_1["z", ] <- as.numeric(level_1["tau", ] > data$T.cal)
    level_1["mu", ] <- 1 / level_1["tau", ]

    ## run MCMC chain ##

    for (step in 1:(burnin + mcmc)) {

      # print
      if (step %% trace == 0){
        time_1 = Sys.time()
        date_time = as.character(Sys.time())
        cat("chain: ", chain_id, " step: ", step, " of ", (burnin + mcmc), " . Datetime: ", date_time, ".....\n", sep = "")
      }

      # Print in log file
      if(step %% 10 == 0){
        cat("chain: ", chain_id, " step: ", step, " of ", (burnin + mcmc), " . Datetime: ", date_time, "\n", sep = "",
            file = log_filename, append = T)
      }


      # store
      if ( (step - burnin) > 0 & (step - 1 - burnin) %% thin == 0) {
        idx <- (step - 1 - burnin) %/% thin + 1
        level_1_draws[idx, , ] <- level_1 # nolint
        level_2_draws[idx, ] <- level_2
      }

      # draw individual-level parameters
      level_1["k", ] <- draw_k(data, level_1, level_2)
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ] <- draw_mu(data, level_1, level_2)
      level_1["tau", ] <- draw_tau(data, level_1, level_2)
      level_1["z", ] <- as.numeric(level_1["tau", ] > data$T.cal)


      # draw heterogeneity parameters
      level_2[c("t", "gamma")] <- draw_gamma_params("k", level_1, level_2, hyper_prior)
      level_2[c("r", "alpha")] <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
      level_2[c("s", "beta")] <- draw_gamma_params("mu", level_1, level_2, hyper_prior)

      time_2 = Sys.time()
      dif_time_mins = as.numeric(difftime(time_2, time_1, units = "mins"))
      unit = ifelse(dif_time_mins < 1, "secs", "mins")
      dif_time = ifelse(dif_time_mins < 1, dif_time_mins*60, dif_time_mins)
      if (step %% trace == 0)
        cat("\tElapsed time:", dif_time, unit, "\n")


    }

    # convert MCMC draws into coda::mcmc objects
    cat("\n\nMCMC finished in chain", chain_id, "\n")
    cat("\n\nMCMC finished in chain", chain_id, "\n", file = log_filename, append = T)
    out = list(
      "level_1" = level_1_draws,
      "level_2" = level_2_draws
    )
    cat("\n\nFinal list created in chain", chain_id, "\n")
    cat("\n\nFinal list created in chain", chain_id, "\n", file = log_filename, append = T)

    return(out)

  }

  # set hyper priors
  hyper_prior <- list(r_1 = 0.001, r_2 = 0.001,
                      alpha_1 = 0.001, alpha_2 = 0.001,
                      s_1 = 0.001, s_2 = 0.001,
                      beta_1 = 0.001, beta_2 = 0.001,
                      t_1 = 0.001, t_2 = 0.001,
                      gamma_1 = 0.001, gamma_2 = 0.001)

  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    try({
      df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)), ]
      param_init <- c(1, 1, BTYD::pnbd.EstimateParameters(df))
      names(param_init) <- c("t", "gamma", "r", "alpha", "s", "beta")
      param_init <- as.list(param_init)
    },
    silent = TRUE)
    if (is.null(param_init))
      param_init <- list(t = 1, gamma = 1, r = 1, alpha = 1, s = 1, beta = 1)
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse = ", "), "\n")
  }

  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal", "litt") %in% names(cal.cbs)))
  stopifnot(all(is.finite(cal.cbs$litt)))

  # run multiple chains - executed in parallel on Unix
  ncores <- ifelse(!is.null(mc.cores), min(chains, mc.cores),
                   ifelse(.Platform$OS.type == "windows", 1, min(chains, detectCores())))

  if (ncores > 1){
    cat("running in parallel on", ncores, "cores\n")
  } else{
    cat("running on 1 core\n")
  }


  draws <- parallel::mclapply(1:chains, function(i) run_single_chain(i, cal.cbs, hyper_prior),
                              mc.cores = ncores,
                              mc.preschedule = F)


  # # merge chains into code::mcmc.list objects
  # cat("\nMerging all chains (", as.character(Sys.time()), ") ...", sep = "")
  # out <- list(level_1 = lapply(1:nrow(cal.cbs), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
  #             level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2)))
  # cat("done (", as.character(Sys.time()), ")\n\n", sep = "")
  # if ("cust" %in% names(cal.cbs))
  #   names(out$level_1) <- cal.cbs$cust
  # return(out)

  cat("\n\nMCMC list created. Returning list as final object.\n\n")
  cat("\n\nMCMC list created. Returning list as final object.\n\n", file = log_filename, append = T)

  return(draws)
}






restructure_samples = function(draws, cust_ids = NULL){
  # Pasa de una lista en la que cada elemento de la lista es una cadena de MCMC a
  # una lista con level_1 y level_2 como en el paquete original

  n_chains = length(draws)
  n_cust = dim(draws[[1]]$level_1)[3]

  message("Doing level 2")
  level_2 = lapply(1:n_chains, function(chain){
    out_tbl_l2 = draws[[chain]]$level_2 %>%
      as_tibble() %>%
      mutate(chain = chain)
    return(out_tbl_l2)
  })
  message("Level 2 completed\n\n")


  message("Doing level 1")
  level_1 = lapply(1:n_cust, function(cust){
    time_1 = Sys.time()

    sims_user = lapply(1:n_chains, function(chain){
      out_l2 = draws[[chain]]$level_1[,,cust]
      # Paste chain number
      out_l2 = cbind(out_l2, rep(chain, nrow(out_l2)))
      colnames(out_l2) = c(colnames(out_l2)[1:5], "chain")
      return(out_l2)
    })
    # Coerce to a single matrix
    sims_user = do.call(rbind, sims_user)

    time_2 = Sys.time()
    if(cust %% 10000 == 1) {
      message("\tCustomer ", cust, " of ", n_cust)
      message("\t", difftime(time_2, time_1, units = "secs"), " secs\n\n")
    }
    return(sims_user)
  })
  message("Level 1 completed\n\n")

  if(!is.null(cust_ids)){
    names(level_1) = cust_ids
  }

  return(list(
    level_1 = level_1,
    level_2 = level_2
  ))

}




get_median_future_trans = function(future_trans){
  future_trans_median = matrixStats::colMedians(future_trans, na.rm = T)
  return(future_trans_median)
}



get_mean_future_trans = function(future_trans){
  future_trans_mean = colMeans(future_trans, na.rm = T)
  return(future_trans_mean)
}





