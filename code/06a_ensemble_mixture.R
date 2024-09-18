

library(tidyverse)
library(glue)
library(rstan)
library(recipes)
rstan_options(auto_write=T)

source("code/00_fn.R")


# ensembling --------------------------------------------------------------

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}"))
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
sims_3D <- filter(sim_i, lab_short=="3D")$sim
sims_2D <- filter(sim_i, lab_short=="2D")$sim

site_i <- read_csv("data/farm_sites.csv")
ensFull_LatLon <- read_csv("out/valid_df.csv") |>
  select(rowNum, date, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)


ensMix_splineD4_rec_all <- recipe(licePerFish_rtrt ~ ., data=ensFull_LatLon) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=4) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD4_rec_all, "out/ensembles/recipe_sLonLatD4_all.rds")

ensMix_splineD4_rec_3D <- recipe(licePerFish_rtrt ~ ., 
                                 data=ensFull_LatLon |>
                                   select(-any_of(sims_2D)) |>
                                   select(-any_of(paste0("c_", sims_2D)))) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=4) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD4_rec_all, "out/ensembles/recipe_sLonLatD4_3D.rds")


ensMix_splineD6_rec_all <- recipe(licePerFish_rtrt ~ ., data=ensFull_LatLon) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=6) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD6_rec_all, "out/ensembles/recipe_sLonLatD6_all.rds")

ensMix_splineD6_rec_3D <- recipe(licePerFish_rtrt ~ ., 
                                 data=ensFull_LatLon |>
                                   select(-any_of(sims_2D)) |>
                                   select(-any_of(paste0("c_", sims_2D)))) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=6) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD6_rec_all, "out/ensembles/recipe_sLonLatD6_3D.rds")

ensMix_splineD8_rec_all <- recipe(licePerFish_rtrt ~ ., data=ensFull_LatLon) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=8) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD8_rec_all, "out/ensembles/recipe_sLonLatD8_all.rds")

ensMix_splineD8_rec_3D <- recipe(licePerFish_rtrt ~ ., 
                               data=ensFull_LatLon |>
                                 select(-any_of(sims_2D)) |>
                                 select(-any_of(paste0("c_", sims_2D)))) |>
  step_interact(~easting:northing) |>
  step_bs(easting, northing, easting_x_northing, deg_free=8) |>
  update_role(c(rowNum, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
  update_role_requirements("id", bake=F) |>
  prep()
saveRDS(ensMix_splineD8_rec_all, "out/ensembles/recipe_sLonLatD8_3D.rds")






# cross-validation --------------------------------------------------------

mods <- expand_grid(mod=c("ranef", "sLonLatD4", "sLonLatD6", "sLonLatD8"),
                    d=c("3D", "all"))

CV_ensMix <- vector("list", length(years))

for(yr in rev(seq_along(years))) {
  CV_yr <- vector("list", nrow(mods))
  for(i in 1:nrow(mods)) {
    mod <- mods$mod[i]
    if(grepl("sLonLat", mods$mod[i])) {
      if(mods$d[i]=="3D") {
        if(grepl("D4", mods$mod[i])) {
          train_df <- ensMix_splineD4_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D))))
          test_df <- ensMix_splineD4_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D))))
        } else if(grepl("D6", mods$mod[i])) {
          train_df <- ensMix_splineD6_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D))))
          test_df <- ensMix_splineD6_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D))))
        } else if(grepl("D8", mods$mod[i])) {
          train_df <- ensMix_splineD8_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D))))
          test_df <- ensMix_splineD8_rec_3D |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(sims_3D), 
                          any_of(paste0("c_", sims_3D)))) 
        }
      } else {
        if(grepl("D4", mods$mod[i])) {
          train_df <- ensMix_splineD4_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
          test_df <- ensMix_splineD4_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
        } else if(grepl("D6", mods$mod[i])) {
          train_df <- ensMix_splineD6_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
          test_df <- ensMix_splineD6_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
        } else if(grepl("D8", mods$mod[i])) {
          train_df <- ensMix_splineD8_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) != years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
          test_df <- ensMix_splineD8_rec_all |>
            bake(ensFull_LatLon |> 
                   filter(year(date) == years[yr]) |>
                   select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                          any_of(c(sims_3D, sims_2D)), 
                          any_of(paste0("c_", c(sims_3D, sims_2D)))))
        }
      }
      dat_rstan <- make_data_rstan_sLonLat(train_df)
      pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                paste0("b_s_", c("easting", "northing", "easting_x_northing")))
    } else {
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
    }
    if(grepl("pRE_huRE", mods$mod[i])) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    }
    if(file.exists(glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_stanfit.rds"))) {
      cat("File exists:", mods$d[i], mods$mod[i], years[yr], "\n")
      out_ensMix <- readRDS(glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_stanfit.rds"))
    } else {
      stanMod <- str_remove(mods$mod[i], "D[0-9]")
      out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                         model_name=mods$mod[i], data=dat_rstan,
                         chains=30, cores=30, iter=2250, warmup=2000,
                         # chains=1, cores=1, iter=22, warmup=20,
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         pars=pars)
      saveRDS(out_ensMix, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_stanfit.rds"))
      saveRDS(dat_rstan, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_CV-{years[yr]}_standata.rds")) 
    }
    if(grepl("sLonLat", mods$mod[i])) {
      CV_yr[[i]] <- make_predictions_ensMix_sLonLat(out_ensMix, test_df, mode="epred", iter=10) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mods$mod[i], "_", mods$d[i]))
    } else {
      CV_yr[[i]] <- make_predictions_ensMix(out_ensMix, test_df, re=T, re_hu=grepl("huRE", mods$mod[i])) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mods$mod[i], "_", mods$d[i])) 
    }
  }
  CV_ensMix[[yr]] <- bind_cols(test_df |> select(rowNum), reduce(CV_yr, bind_cols))
}
reduce(CV_ensMix, bind_rows) |>
  write_csv("out/ensembles/CV_ensMix_predictions.csv")




# full dataset ------------------------------------------------------------

mods <- expand_grid(mod=c("ranef", "sLonLatD4", "sLonLatD6", "sLonLatD8"),
                    d=c("3D", "all")) |>
  arrange(desc(mod), d)

for(i in 1:nrow(mods)) {
  mod <- mods$mod[i]
  if(grepl("sLonLat", mods$mod[i])) {
    if(mods$d[i]=="3D") {
      if(grepl("D4", mods$mod[i])) {
        train_df <- ensMix_splineD4_rec_3D |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(sims_3D), 
                        any_of(paste0("c_", sims_3D))))
      } else if(grepl("D6", mods$mod[i])) {
        train_df <- ensMix_splineD6_rec_3D |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(sims_3D), 
                        any_of(paste0("c_", sims_3D))))
      } else if(grepl("D8", mods$mod[i])){
        train_df <- ensMix_splineD8_rec_3D |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(sims_3D), 
                        any_of(paste0("c_", sims_3D))))
      }
    } else {
      if(grepl("D4", mods$mod[i])) {
        train_df <- ensMix_splineD4_rec_all |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(c(sims_3D, sims_2D)), 
                        any_of(paste0("c_", c(sims_3D, sims_2D)))))
      } else if(grepl("D6", mods$mod[i])) {
        train_df <- ensMix_splineD6_rec_all |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(c(sims_3D, sims_2D)), 
                        any_of(paste0("c_", c(sims_3D, sims_2D)))))
      } else if(grepl("D8", mods$mod[i])) {
        train_df <- ensMix_splineD8_rec_all |>
          bake(ensFull_LatLon |> 
                 select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                        any_of(c(sims_3D, sims_2D)), 
                        any_of(paste0("c_", c(sims_3D, sims_2D)))))
      }
    }
    dat_rstan <- make_data_rstan_sLonLat(train_df)
    pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              paste0("b_s_", c("easting", "northing", "easting_x_northing")))
  } else {
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
  }
  if(!file.exists(glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_FULL_stanfit.rds"))) {
    stanMod <- str_remove(mods$mod[i], "D[0-9]")
    out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                       model_name=mods$mod[i], data=dat_rstan,
                       chains=30, cores=15, iter=2250, warmup=2000,
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       pars=pars)
    saveRDS(out_ensMix, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_FULL_stanfit.rds"))
    saveRDS(dat_rstan, glue("out/ensembles/ensMix_{mods$d[i]}_{mods$mod[i]}_FULL_standata.rds"))
  } else {
    cat("File exists:", mods$d[i], mods$mod[i], "\n")
    Sys.sleep(3)
  }
}
