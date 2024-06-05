

library(tidyverse)
library(glue)
library(rstan)
rstan_options(auto_write=T)

source("code/00_fn.R")

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv")

sim_i <- dir("out/biotracker", "sim_i.csv", recursive=T, full.names=T) |>
  map_dfr(~read_csv(.x) |> select(i, fixDepth) |> mutate(path=.x)) |>
  filter((grepl("full", path) & i %in% paste0("0", 1:6)) |
           (grepl("eggT", path) & i %in% str_pad(7:12, 2, "left", "0")) | 
           (grepl("swim", path) & i %in% 13:20)) |>
  arrange(i) |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  filter(sim %in% names(ensFull_df)) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  select(sim, lab_short, lab)
sim_i <- sim_i[17:20,]


# individual simulation models --------------------------------------------

fit_df <- ensFull_df |> 
  filter(year < 2023) |>
  select(rowNum, sepaSite, sepaSiteNum, productionCycleNumber, date, CoGP_period,
         licePerFish_rtrt, liceTreat, any_of(sim_i$sim), any_of(paste0("c_", sim_i$sim)))

for(i in 1:nrow(sim_i)) {
  sim <- sim_i$sim[i]
  dat_rstan <- make_data_rstan(fit_df |> select(-all_of(sim_i$sim[-i]), -all_of(paste0("c_", sim_i$sim[-i]))))
  out_sim <- stan(file="code/stan/candidate_model.stan",
                  model_name=sim,
                  data=dat_rstan,
                  chains=6, cores=6,
                  iter=3000, warmup=2500,
                  pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma"))
  saveRDS(out_sim, glue("out/{sim}_stanfit.rds"))
}
for(i in 1:nrow(sim_i)) {
  sim <- sim_i$sim[i]
  dat_rstan <- make_data_rstan(fit_df |> select(-all_of(sim_i$sim[-i]), -all_of(paste0("c_", sim_i$sim[-i]))))
  out_sim <- stan(file="code/stan/candidate_model_ranef.stan",
                  model_name=sim,
                  data=dat_rstan,
                  chains=6, cores=6,
                  iter=3000, warmup=2500,
                  pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "r_grp_hu", "sd_grp_hu"))
  saveRDS(out_sim, glue("out/{sim}_RE_stanfit.rds"))
}




# unweighted mean models --------------------------------------------------

fit_df <- ensFull_df |>
  filter(year < 2023) |>
  select(rowNum, sepaSite, sepaSiteNum, productionCycleNumber, date, CoGP_period,
         licePerFish_rtrt, liceTreat, starts_with("sim_avg"), starts_with("c_sim_avg"))

sim_avgs <- c("sim_avg2D", "sim_avg3D")
for(i in 1:length(sim_avgs)) {
  sim <- sim_avgs[i]
  dat_rstan <- make_data_rstan(fit_df |> select(-all_of(sim_avgs[-i]), -all_of(paste0("c_", sim_avgs[-i]))))
  out_sim <- stan(file="code/stan/candidate_model.stan",
                  model_name=sim,
                  data=dat_rstan,
                  chains=6, cores=6,
                  iter=3000, warmup=2500,
                  pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma"))
  saveRDS(out_sim, glue("out/{sim}_stanfit.rds"))
}
for(i in 1:length(sim_avgs)) {
  sim <- sim_avgs[i]
  dat_rstan <- make_data_rstan(fit_df |> select(-all_of(sim_avgs[-i]), -all_of(paste0("c_", sim_avgs[-i]))))
  out_sim <- stan(file="code/stan/candidate_model_ranef.stan",
                  model_name=sim,
                  data=dat_rstan,
                  chains=6, cores=6,
                  iter=3000, warmup=2500,
                  pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "r_grp_hu", "sd_grp_hu"))
  saveRDS(out_sim, glue("out/{sim}_RE_stanfit.rds"))
}
