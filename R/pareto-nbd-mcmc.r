
# TODO: normalize time so that T.cal~1


#' Hierarchical Bayes variant of Pareto/NBD
#'
#' \code{pnbd.mcmc.DrawParameters} samples parameters via MCMC for a given CBS
#' matrix
#' 
#' method 1) If \code{use_data_augmentation==TRUE} then implementation follows chapter
#' 3.2 of Sandeep Conoor's dissertation
#' http://gradworks.umi.com/34/02/3402149.html, i.e. parameter space is expanded
#' with unobserved lifetime tau_i. Note, however, that  we follow the notation
#' of original SMC paper with alpha and beta being the 'rate' parameters of the
#' Gamma distributions, instead of 'scale' parameter.
#' 
#' method 2) If \code{use_data_augmentation==FALSE} then implementation follows
#' Shao-Hui Ma & Jin-Lan Liu paper 
#' http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404, i.e. no data
#' augmentation and draws on individual level need to be done via slice
#' sampling. As such it is 3x slower than method 1)
#' 
#' Estimating parameters via Pareto/NBD MCMC can be 10x slower than Pareto/NBD
#' MLE, which itself can be 10x slower than BG/NBD. Both methods exhibit highly
#' autocorrelated draws of {r, alpha, s, beta} and hence need to be run long, to
#' generate 'enough' draws
#'
#' @param data data.frame with columns 'x', 't.x', 'T.cal'
#' @param mcmc number of MCMC steps
#' @param burnin number of initial MCMC steps which are discarded
#' @param thin only every thin-th MCMC step will be returned
#' @param chains number of MCMC chains to be run
#' @param use_data_augmentation determines MCMC method to be used
#' @param param_init list of 2nd-level parameter start values
#' @param hyper_prior list of hyper parameters for 2nd-level parameters
#' @return list
#' @import coda Rcpp parallel
#' @export
#' @examples
#' #params <- list(r=1.4, alpha=1.3, s=0.7, beta=7)
#' #cbs <- pcnbd.GenerateData(1000, 10, 5, params)$cbs
#' #draws <- pnbd.mcmc.DrawParameters(cbs, mcmc=10000, burnin=10000, thin=10, chains=2)
#' #plot(draws$level_2)
#' #rbind("actual"=unlist(params), "estimated"=summary(draws$level_2, quantiles=0.5)$quantiles)
#' @seealso pcnbd.GenerateData
#' @references http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404
#' @references http://gradworks.umi.com/34/02/3402149.html
pnbd.mcmc.DrawParameters <-
  function(data,
           mcmc = 10000, burnin = 0, thin = 1, chains = 1,
           use_data_augmentation = TRUE,
           param_init = list(r=1, alpha=1, s=1, beta=1),
           hyper_prior = list(r_1=1/1000, r_2=1/1000,
                              alpha_1=1/1000, alpha_2=1/1000,
                              s_1=1/1000, s_2=1/1000,
                              beta_1=1/1000, beta_2=1/1000)) {
  ## Estimated model parameters for given data by running MCMC chain
  ## 
  ##
  ## Returns:
  ##   2-element list
  ##     level_1:  list of coda::mcmc.list objects; one for each customer, containing individual-level draws
  ##     level_2:  coda::mcmc.list object containing draws of heterogeneity parameters
  

  ## methods to sample heterogeneity parameters {r, alpha, s, beta} ##
    
  draw_shape <- function(x, shape, rate, prior1, prior2, steps=20) {

    calc_prior <- function(r) (prior1 - 1) * log(r) - (r * prior2)
    calc_likel <- function(r) (r - 1) * sum(log(x)) + length(x) * (r * log(rate) - lgamma(r))

    cur_shape <- shape
    cur_prior <- calc_prior(cur_shape)
    cur_likel <- calc_likel(cur_shape)
    
    scale <- mean(x)
    
    for (i in 1:steps) {
      
      # generate new proposal
      new_shape <- cur_shape * exp(max(-100, min(100, scale * rt(1, df=3))))
      new_prior <- calc_prior(new_shape)
      new_likel <- calc_likel(new_shape)

      # accept/reject new proposal
      mhratio <- exp(new_prior+new_likel-cur_prior-cur_likel)
      if (mhratio > runif(n=1))
        cur_shape <- new_shape

    }
    return(cur_shape)
  }

  draw_rate <- function(x, shape, rate, prior1, prior2) {
    rgamma(n     = 1, 
           shape = prior1 + length(x) * shape,
           rate  = prior2 + sum(x))
  }

  draw_r <- function(level_1, level_2, hyper_prior) {
    draw_shape(x      = level_1["lambda",], 
               shape  = level_2["r"],
               rate   = level_2["alpha"],
               prior1 = hyper_prior$r_1,
               prior2 = hyper_prior$r_2)
  }
  
  draw_alpha <- function(level_1, level_2, hyper_prior) {
    draw_rate(x      = level_1["lambda",],
              shape  = level_2["r"],
              rate   = level_2["alpha"],
              prior1 = hyper_prior$alpha_1,
              prior2 = hyper_prior$alpha_2)
  }
  
  draw_s <- function(level_1, level_2, hyper_prior) {
    draw_shape(x      = level_1["mu",], 
               shape  = level_2["s"],
               rate   = level_2["beta"],
               prior1 = hyper_prior$s_1,
               prior2 = hyper_prior$s_2)
  }

  draw_beta <- function(level_1, level_2, hyper_prior) {
    draw_rate(x      = level_1["mu",],
              shape  = level_2["s"],
              rate   = level_2["beta"],
              prior1 = hyper_prior$beta_1,
              prior2 = hyper_prior$beta_2)
  }
  
  ## methods to sample individual-level parameters (with data augmentation) ##  
  
  draw_lambda_conoor <- function(data, level_1, level_2) {
    N      <- nrow(data)
    x      <- data[, "x"]
    T.cal  <- data[, "T.cal"]
    tau    <- level_1["tau", ]
    r      <- level_2["r"]
    alpha  <- level_2["alpha"]
    
    lambda <- rgamma(n     = N,
                     shape = r + x,
                     rate  = alpha + pmin(tau, T.cal))
    lambda[lambda==0 | log(lambda) < -70] <- exp(-70) # avoid numeric overflow
    return(lambda)
  }
  
  draw_mu_conoor <- function(data, level_1, level_2) {
    N      <- nrow(data)
    tau    <- level_1["tau", ]  
    s      <- level_2["s"]
    beta   <- level_2["beta"]
    
    mu <- rgamma(n     = N, 
                 shape = s + 1, 
                 rate  = beta + tau)
    mu[mu==0 | log(mu) < -70] <- exp(-70) # avoid numeric overflow
    return(mu)
  }
  
  draw_tau <- function(data, level_1) {
    N      <- nrow(data)
    tx     <- data[, "t.x"]
    T.cal  <- data[, "T.cal"]
    lambda <- level_1["lambda", ]
    mu     <- level_1["mu", ]
  
    mu_lam <- mu + lambda
    t_diff <- T.cal - tx
    
    # sample z
    p_alive <- 1 / (1+(mu/mu_lam)*(exp(mu_lam*t_diff)-1))
    alive   <- p_alive > runif(n=N)
    
    tau <- numeric(N)
    
    # Case: still alive - left truncated exponential distribution -> [T.cal, Inf]
    if (any(alive)) {
      tau[alive]  <- T.cal[alive] + rexp(sum(alive), mu[alive])
    }
    
    # Case: churned     - double truncated exponential distribution -> [tx, T.cal]
    if (any(!alive)) {
      mu_lam_tx   <- pmin(700, mu_lam[!alive] * tx[!alive])
      mu_lam_Tcal <- pmin(700, mu_lam[!alive] * T.cal[!alive])
      # sample with http://en.wikipedia.org/wiki/Inverse_transform_sampling
      rand        <- runif(n=sum(!alive))
      tau[!alive] <- -log( (1-rand)*exp(-mu_lam_tx) + rand*exp(-mu_lam_Tcal)) / mu_lam[!alive]
    }
    
    return(tau)
  }
  
  ## methods to sample individual-level parameters (without data augmentation) ##  
  
  draw_lambda_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("lambda", 
                        x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], 
                        lambda = level_1["lambda",], mu = level_1["mu",], 
                        r = level_2["r"], alpha = level_2["alpha"], 
                        s = level_2["s"], beta = level_2["beta"])
  }
  
  draw_mu_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("mu", 
                        x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], 
                        lambda = level_1["lambda",], mu = level_1["mu",], 
                        r = level_2["r"], alpha = level_2["alpha"], 
                        s = level_2["s"], beta = level_2["beta"])
  }
  
  run_single_chain <- function(chain_id=1) {
    
    ## initialize arrays for storing draws ##
    
    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc-1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim=c(nr_of_draws, 4))
    dimnames(level_2_draws)[[2]] <- c("r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim=c(nr_of_draws, 3, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("lambda", "mu", "tau")
    
    ## initialize parameters ##
    
    level_2          <- level_2_draws[1,]
    level_2["r"]     <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"]     <- param_init$s
    level_2["beta"]  <- param_init$beta
    
    level_1            <- level_1_draws[1,,]
    level_1["lambda",] <- mean(data$x) / mean(ifelse(data$t.x==0, data$T.cal, data$t.x))
    level_1["tau",]    <- data$t.x + 0.5/level_1["lambda",]
    level_1["mu",]     <- 1/level_1["tau",]  
    
    ## run MCMC chain ##
    
    for (step in 1:(burnin+mcmc)) {
      if (step%%100==0) cat("chain:", chain_id, "step:", step, "of", (burnin+mcmc), "\n")
  
      # store
      if ((step-burnin)>0 & (step-1-burnin)%%thin==0) {
        idx <- (step-1-burnin)%/%thin + 1
        level_1_draws[idx,,] <- level_1
        level_2_draws[idx,]  <- level_2
      }
  
      # draw individual-level parameters
      draw_lambda <- if (use_data_augmentation) draw_lambda_conoor else draw_lambda_ma_liu
      draw_mu     <- if (use_data_augmentation) draw_mu_conoor else draw_mu_ma_liu
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ]     <- draw_mu(data, level_1, level_2)
      level_1["tau", ]    <- draw_tau(data, level_1)
      
      # draw heterogeneity parameters
      level_2["alpha"]    <- draw_alpha(level_1, level_2, hyper_prior)
      level_2["r"]        <- draw_r(level_1, level_2, hyper_prior)
      level_2["beta"]     <- draw_beta(level_1, level_2, hyper_prior)
      level_2["s"]        <- draw_s(level_1, level_2, hyper_prior)
    }

    # convert MCMC draws into coda::mcmc objects
    return(list(
      level_1 = lapply(1:nr_of_cust, function(i) mcmc(level_1_draws[,,i], start=burnin, thin=thin)),
      level_2 = mcmc(level_2_draws, start=burnin, thin=thin)))
  }

  # run multiple chains
  draws <- lapply(1:chains, function(i) run_single_chain(i))
  
  # merge chains into code::mcmc.list objects
  return(list(
    level_1 = lapply(1:nrow(data), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
    level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2))))
}
