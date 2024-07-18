

library(tidyverse)
library(glue)
library(rstan)
rstan_options(auto_write=T)

source("code/00_fn.R")


# ensembling --------------------------------------------------------------

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv")
years <- unique(ensFull_df$year)

sim_i <- read_csv("out/sim_2019-2023/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predF1", "predF5", "predM", "predMRE", "sim_avg2D", "sim_avg3D", "null"),
           lab_short=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"),
           lab=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
                             paste0("2D.", 1:20), paste0("3D.", 1:20), "Null")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
                                   "2D", "3D", "Null")))




# cross-validation --------------------------------------------------------

mods <- expand_grid(mod=c("ranef"),
                    d=c("3D"))

CV_ensMix <- vector("list", length(years))

for(yr in seq_along(years)) {
  CV_yr <- vector("list", nrow(mods))
  for(i in 1:nrow(mods)) {
    mod <- mods$mod[i]
    train_df <- ensFull_df |>
      filter(year != years[yr]) |>
      select(rowNum, sepaSite, sepaSiteNum, date, licePerFish_rtrt,
             any_of(filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim), 
             any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim)))
    test_df <- ensFull_df |>
      filter(year == years[yr]) |>
      select(rowNum, sepaSite, sepaSiteNum, date, licePerFish_rtrt, 
             any_of(filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim), 
             any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim)))
    dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
    pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              "b_p_uc", "r_grp", "sd_grp")
    if(grepl("pRE_huRE", mods$mod[i])) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    }
    out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{mods$mod[i]}.stan"),
                       model_name=mods$mod[i], data=dat_rstan,
                       chains=3, cores=3, iter=2000, warmup=1500,
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       pars=pars)
    saveRDS(out_ensMix, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_stanfit.rds"))
    saveRDS(dat_rstan, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_standata.rds"))
    CV_yr[[i]] <- make_predictions_ensMix(out_ensMix, test_df, re=T, re_hu=F) |>
      colMeans() |>
      as_tibble() |>
      set_names(paste0("IP_", mods$mod[i], "_", mods$d[i]))
  }
  CV_ensMix[[yr]] <- bind_cols(test_df |> select(rowNum), reduce(CV_yr, bind_cols))
}
reduce(CV_ensMix, bind_rows) |>
  write_csv("out/ensembles/CV_ensMix_predictions.csv")




# full dataset ------------------------------------------------------------

mods <- expand_grid(mod=c("ranef"),
                    d=c("3D"))

for(i in 1:nrow(mods)) {
  mod <- mods$mod[i]
  train_df <- ensFull_df |>
    select(rowNum, sepaSite, sepaSiteNum, date, licePerFish_rtrt,
           any_of(filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim), 
           any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(mods$d[i]=="all", "2D", "")))$sim)))
  dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
  pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
            "b_p_uc", "r_grp", "sd_grp")
  if(grepl("pRE_huRE", mods$mod[i])) {
    pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
  }
  out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{mods$mod[i]}.stan"),
                     model_name=mods$mod[i], data=dat_rstan,
                     chains=3, cores=3, iter=4000, warmup=3000,
                     control=list(adapt_delta=0.95, max_treedepth=20),
                     pars=pars)
  saveRDS(out_ensMix, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_FULL_stanfit.rds"))
  saveRDS(dat_rstan, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_FULL_standata.rds"))
}
