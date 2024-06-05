


make_data_rstan <- function(df, R2D2_sd=TRUE) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(N=nrow(df),
                    Y=df$licePerFish_rtrt,
                    N_Ye0=sum(zeros),
                    indexes_Ye0=which(zeros),
                    N_Yg0=sum(!zeros),
                    indexes_Yg0=which(!zeros),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(all_of(sim_names))),
                    Xc=as.matrix(df |> select(all_of(paste0("c_", sim_names)))),
                    N_groups=max(df$sepaSiteNum),
                    J_group=df$sepaSiteNum,
                    sim_names=sim_names,
                    R2D2_mean_R2=0.5,
                    R2D2_prec_R2=2,
                    R2D2_cons_D2=rep(0.1, length(sim_names)*ifelse(R2D2_sd, 2, 1)),
                    R2D2_mean_R2_grp=0.5,
                    R2D2_prec_R2_grp=2,
                    R2D2_cons_D2_grp=rep(0.1, length(sim_names)),
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    prior_only=0)
  return(dat_rstan)
}




make_predictions_ensMix <- function(out, newdata, iter=2000, seed=NULL, mode="epred", re=FALSE) {
  library(tidyverse); library(rstan)
  # hurdle component is fitted with centered IP
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]
  
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan(newdata)
  ensIP <- matrix(0, nrow=iter, ncol=nrow(newdata))
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  if(re) {
    r_grp <- rstan::extract(out, pars="r_grp")[[1]]
    r_grp_hu <- rstan::extract(out, pars="r_grp_hu")[[1]]
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(r_grp[iters, dat_ls$J_group[i],,drop=T])
      hu_RE[,i] <- cbind(1, dat_ls$Xc[i,,drop=F]) %*% t(r_grp_hu[iters, dat_ls$J_group[i],])
    }
  } else {
    b_p <- rstan::extract(out, pars="b_p")[[1]]
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(b_p[iters,,drop=F])
    }
  }
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * c(ensIP[,i])
    hu <- Intercept_hu[iters] + c(dat_ls$Xc[i,,drop=F] %*% t(b_hu[iters,,drop=F])) + hu_RE[,i]
    # hu = pr(0)
    if(mode=="epred") {
      preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
    } else if(mode=="predict") {
      preds[,i] <- (1 - rbinom(iter, 1, brms::inv_logit_scaled(hu))) * 
        rlnorm(iter, mu, sigma[iters]) 
    }
  }
  return(preds)
}

make_predictions_candidate <- function(out, newdata, sim, iter=2000, seed=NULL, re=FALSE) {
  library(tidyverse); library(rstan)
  
  newdata <- newdata |> 
    select(sepaSite, sepaSiteNum, date, licePerFish_rtrt, matches(sim))
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]
  
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan(newdata)
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  if(re) {
    r_grp_hu <- rstan::extract(out, pars="r_grp_hu")[[1]]
    for(i in 1:ncol(preds)) {
      hu_RE[,i] <- cbind(1, dat_ls$Xc[i,,drop=F]) %*% t(r_grp_hu[iters, dat_ls$J_group[i],])
    }
  } 
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * dat_ls$X[i,]
    hu <- Intercept_hu[iters] + b_hu[iters,] * dat_ls$Xc[i,] + hu_RE[,i]
    # hu = pr(0)
    preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
  }
  return(preds)
}

