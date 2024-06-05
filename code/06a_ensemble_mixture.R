

library(tidyverse)
library(glue)
library(rstan)
rstan_options(auto_write=T)

source("code/00_fn.R")


# ensembling --------------------------------------------------------------

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

fit_df <- ensFull_df |> 
  filter(year < 2023) |>
  select(rowNum, sepaSite, sepaSiteNum, productionCycleNumber, date, CoGP_period,
         licePerFish_rtrt, liceTreat, 
         any_of(filter(sim_i, lab_short %in% c("3D"))$sim), 
         any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D"))$sim)))


# stan model --------------------------------------------------------------

# dat_rstan <- make_data_rstan(fit_df)
# out_ensMix <- stan(file="code/stan/ensemble_mixture_model.stan",
#                    model_name="ensMix",
#                    data=dat_rstan,
#                    chains=6, cores=6, iter=3000, warmup=2500,
#                    control=list(max_treedepth=20),
#                    pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu"))
# saveRDS(out_ensMix, "out/ensMix_3D_stanfit.rds")
# bayesplot::mcmc_intervals(out_ensMix, regex_pars="b_p") + ggtitle("noRE") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

mod <- "pRE_huRE"
dat_rstan <- make_data_rstan(fit_df, R2D2_sd=F)
out_ensMixRE <- stan(file=glue("code/stan/ensemble_mixture_model_{mod}.stan"),
                          model_name=mod,
                          data=dat_rstan,
                          chains=6, cores=6, iter=3000, warmup=2500,
                          control=list(adapt_delta=0.9, max_treedepth=20),
                          pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                                 "b_p_uc", "r_grp", "sd_grp", "r_grp_hu", "sd_grp_hu"))
saveRDS(out_ensMixRE, glue("out/ensMix_3D_{mod}_stanfit.rds"))
saveRDS(dat_rstan, glue("out/ensMix_3D_{mod}_standata.rds"))
bayesplot::mcmc_intervals(out_ensMixRE, regex_pars="b_p\\[") + ggtitle(mod) +
  scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

dat_rstan <- make_data_rstan(fit_df, R2D2_sd=F)
out_ensMixRE_R2D2 <- stan(file="code/stan/ensemble_mixture_model_pRE_huRE-R2D2.stan",
                      model_name="ensMix_pRE_huRE-R2D2",
                      data=dat_rstan,
                      chains=6, cores=6, iter=3000, warmup=2500,
                      control=list(adapt_delta=0.9, max_treedepth=20),
                      pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", 
                             "b_p_uc", "r_grp", "sd_grp", "r_grp_hu", "sd_grp_hu"))
saveRDS(out_ensMixRE_R2D2, "out/ensMix_3D_pRE_huRE-R2D2_stanfit.rds")
bayesplot::mcmc_intervals(out_ensMixRE_R2D2, regex_pars="b_p\\[") + ggtitle("ensMix_pRE_huRE-R2D2") +
  scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

# dat_rstan <- make_data_rstan(fit_df)
# out_ensMixRE <- stan(file="code/stan/ensemble_mixture_model_ranef.stan",
#                    model_name="ensMix_ranef",
#                    data=dat_rstan,
#                    chains=6, cores=6, iter=3000, warmup=2500,
#                    control=list(adapt_delta=0.9),
#                    pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p_uc", "r_grp", "sd_grp"))
# saveRDS(out_ensMixRE, "out/ensMixRE_3D_stanfit.rds")
# bayesplot::mcmc_intervals(out_ensMixRE, regex_pars="b_p\\[") + ggtitle("RE") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

# dat_rstan <- make_data_rstan(fit_df)
# out_ensMixRE2 <- stan(file="code/stan/ensemble_mixture_model_ranef_2.stan",
#                      model_name="ensMix_ranef_2",
#                      data=dat_rstan,
#                      chains=6, cores=6, iter=3000, warmup=2500,
#                      control=list(adapt_delta=0.9),
#                      pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p_uc", "r_grp", "sd_grp"))
# saveRDS(out_ensMixRE2, "out/ensMixRE2_3D_stanfit.rds")
# sd_grp needs to be limited -- values up to 200!
# bayesplot::mcmc_intervals(out_ensMixRE2, regex_pars="b_p\\[") + ggtitle("RE2") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

# dat_rstan <- make_data_rstan(fit_df)
# out_ensMixRE3 <- stan(file="code/stan/ensemble_mixture_model_ranef_3.stan",
#                       model_name="ensMix_ranef_3",
#                       data=dat_rstan,
#                       chains=6, cores=6, iter=3000, warmup=2500,
#                       control=list(adapt_delta=0.9),
#                       pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p_uc", "r_grp", "sd_grp"))
# saveRDS(out_ensMixRE3, "out/ensMixRE3_3D_stanfit.rds")
# bayesplot::mcmc_intervals(out_ensMixRE3, regex_pars="b_p\\[") + ggtitle("RE3") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

# dat_rstan <- make_data_rstan(fit_df, R2D2_sd=F)
# out_ensMixRE4 <- stan(file="code/stan/ensemble_mixture_model_ranef_4.stan",
#                       model_name="ensMix_ranef_4",
#                       data=dat_rstan,
#                       chains=6, cores=6, iter=3000, warmup=2500,
#                       control=list(adapt_delta=0.9),
#                       pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p_uc", "r_grp", "sd_grp"))
# saveRDS(out_ensMixRE4, "out/ensMixRE4_3D_stanfit.rds")
# bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="b_p\\[") + ggtitle("RE4") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)

# dat_rstan <- make_data_rstan(fit_df, R2D2_sd=F)
# out_ensMixRE5 <- stan(file="code/stan/ensemble_mixture_model_ranef_5.stan",
#                       model_name="ensMix_ranef_5",
#                       data=dat_rstan,
#                       chains=6, cores=6, iter=3000, warmup=2500,
#                       control=list(adapt_delta=0.9, max_treedepth=20),
#                       pars=c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p_uc", "r_grp", "sd_grp", "r_grp_hu", "sd_grp_hu"))
# saveRDS(out_ensMixRE5, "out/ensMixRE5_stanfit.rds")
# bayesplot::mcmc_intervals(out_ensMixRE5, regex_pars="b_p\\[") + ggtitle("RE5") +
#   scale_y_discrete(labels=dat_rstan$sim_names) + xlim(0,1)



bayesplot::mcmc_intervals(out_ensMixRE, regex_pars="sd_grp\\[") + ggtitle("RE") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE2, regex_pars="sd_grp") + ggtitle("RE2") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE3, regex_pars="sd_grp") + ggtitle("RE3") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="sd_grp") + ggtitle("RE4") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE5, regex_pars="sd_grp\\[") + ggtitle("RE5") +
  scale_y_discrete(labels=dat_rstan$sim_names)



bayesplot::mcmc_intervals(out_ensMix, regex_pars="b_hu\\[") + ggtitle("no RE") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE, regex_pars="b_hu\\[") + ggtitle("RE") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE2, regex_pars="b_hu\\[") + ggtitle("RE2") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE3, regex_pars="b_hu\\[") + ggtitle("RE3") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="b_hu\\[") + ggtitle("RE4") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE5, regex_pars="b_hu\\[") + ggtitle("RE5") +
  scale_y_discrete(labels=dat_rstan$sim_names)


# out_ensMix <- readRDS("out/ensMix_stanfit.rds")
bayesplot::mcmc_intervals(out_ensMix, regex_pars="b_p") +
  scale_y_discrete(labels=grep("sim", names(ensFull_df), value=T))

bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="b_p\\[") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="b_p_uc") +
  scale_y_discrete(labels=dat_rstan$sim_names)
bayesplot::mcmc_intervals(out_ensMixRE4, regex_pars="sd_grp") +
  scale_y_discrete(labels=dat_rstan$sim_names)

bayesplot::mcmc_intervals(out_ensMixRE2, regex_pars="b_hu\\[") +
  scale_y_discrete(labels=dat_rstan$sim_names)

library(ggdist); library(scico)
b_p_post <- rstan::extract(out_ensMixRE2, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_rstan$sim_names, sim_i$sim)]) |>
  pivot_longer(everything(), names_to="Simulation", values_to="p") |>
  mutate(lab_short=str_sub(Simulation, 1, 2),
         Simulation=factor(Simulation, levels=levels(sim_i$lab)))
ggplot(b_p_post, aes(p, Simulation, colour=lab_short, fill=lab_short)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, 
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1) +
  scale_colour_manual(values=scico(2, begin=0.2, end=0.7, palette="broc", direction=-1)) +
  scale_fill_manual(values=scico(2, begin=0.2, end=0.7, palette="broc", direction=-1)) +
  scale_y_discrete(limits=rev(sort(unique(b_p_post$Simulation)))) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  labs(x=expression(paste("Mixing weight (", italic(pi[k]), ")"))) +
  theme(legend.position="none")


rgrp_post <- rstan::extract(out_ensMixRE, "r_grp")[[1]]
rgrp_df <- expand_grid( 
  sim=1:dim(rgrp_post)[3],
  sepaSiteNum=1:dim(rgrp_post)[2],
  iter=1:dim(rgrp_post)[1]
) |> 
  mutate(rgrp=c(rgrp_post))
rgrp_df |>
  ggplot(aes(rgrp)) + 
  geom_density() + 
  facet_wrap(~sim, scales="free_y")
rgrp_df |>
  group_by(sim, sepaSiteNum) |>
  summarise(mn=mean(rgrp)) |>
  ggplot(aes(mn, sim)) + geom_point(alpha=0.25)
rgrp_df |>
  group_by(sim, sepaSiteNum) |>
  summarise(mn=mean(rgrp)) |>
  ggplot(aes(mn, sepaSiteNum)) + geom_point() + facet_grid(.~sim)
# 
# 
# 
# ensTest_df <- ensFull_df |> 
#   mutate(sepaSiteNum=as.numeric(factor(sepaSite))) |>
#   filter(year==2023) |>
#   rowwise() |>
#   mutate(IP_avg3D=mean(c_across(any_of(filter(sim_i, lab_short=="3D")$sim))),
#          IP_avg2D=mean(c_across(any_of(filter(sim_i, lab_short=="2D")$sim)))) |>
#   ungroup() 
# 
# ensTest_df <- ensTest_df |>
#   mutate(IP_ensMix=make_predictions_rstan(out_ensMix, 
#                                           newdata=ensTest_df |>
#                                             mutate(sepaSiteNum=as.numeric(factor(sepaSite))) |>
#                                             rename_with(.fn=~str_remove(.x, "IP_"), .cols=starts_with("IP_sim")) |> select(all_of(names(fit_df)))) |>
#            colMeans())
# 
# ensTest_df |>
#   rename_with(.fn=~glue("IP_{.x}"), .cols=starts_with("sim")) |>
#   pivot_longer(starts_with("IP_"), names_to="sim") |>
#   mutate(sim=str_remove(sim, "IP_")) |>
#   group_by(sim) |>
#   summarise(rsq=yardstick::rsq_vec(value, truth=licePerFish_rtrt),
#             r=cor(value, licePerFish_rtrt, use="pairwise"),
#             N=n(),
#             prop_treat=mean(liceTreat=="TRUE"),
#             prop_0=mean(licePerFish_rtrt==0)) |>
#   ungroup() |>
#   mutate(dates="all") |> 
#   arrange(desc(rsq))
