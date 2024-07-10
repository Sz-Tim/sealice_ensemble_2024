# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk




# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(rstan)
library(yardstick)
library(ggpubr)
library(cowplot)
library(scico)
library(ggdist)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R")

cmr <- readRDS("../../00_misc/cmr_cmaps.RDS")

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(liceTreat=factor(liceTreat),
         lice_g05=factor(licePerFish_rtrt^4 > 0.5))

site_i <- read_csv("data/farm_sites.csv") 
sim_i <- read_csv("out/sim_2019-2023/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predF1", "predF5", "predMix", "sim_avg2D", "sim_avg3D", "null"),
           lab_short=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"),
           lab=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "Mean2D", "Mean3D", 
                             paste0("2D.", 1:20), paste0("3D.", 1:20), "Null")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "Mean2D", "Mean3D", 
                                   "2D", "3D", "Null")))



# ensemble results --------------------------------------------------------

ensTest_df <- ensFull_df |>
  filter(year==2023) |>
  group_by(sepaSite) |>
  mutate(nObs=sum(licePerFish_rtrt > 0)) |>
  ungroup() |>
  filter(nObs > 5) |>
  mutate(week=floor(week(date)/2)) |>
  left_join(ensFull_df |>
              filter(year < 2023) |>
              mutate(week=floor(week(date)/2)) |>
              group_by(week) |>
              summarise(IP_null=mean(licePerFish_rtrt)) |>
              ungroup()) |>
  select(-week)

# Mixing ensemble
out_ensMix_3D <- readRDS("out/ensembles/ensMix_3D_ranef_stanfit.rds")

# Forecasting ensemble validation predictions 
out_ensFcst <- full_join(
  dirf("out/ensembles", "^lice.*_best_preds_1wk") |>
    map(~readRDS(.x) |> 
          select(.row, any_of(c(".pred", ".pred_TRUE"))) |>
          group_by(.row) |> slice_tail(n=1) |> ungroup()) |>
    reduce(full_join) |>
    rename(IP_predF1=.pred, pr_predF1=.pred_TRUE),
  dirf("out/ensembles", "^lice.*_best_preds_5wk") |>
    map(~readRDS(.x) |> 
          select(.row, any_of(c(".pred", ".pred_TRUE"))) |>
          group_by(.row) |> slice_tail(n=1) |> ungroup()) |>
    reduce(full_join) |>
    rename(IP_predF5=.pred, pr_predF5=.pred_TRUE)) |>
  rename(rowNum=.row)
  
# Merge predictions
ensTest_df <- ensTest_df |>
  mutate(
    IP_predMix=make_predictions_ensMix(
      out_ensMix_3D,
      newdata=ensTest_df |>
        select(-contains("sim_avg"), -any_of(filter(sim_i, lab_short=="2D")$sim)),
      re=T, re_hu=F) |>
        colMeans(),
    ) |>
  left_join(out_ensFcst)

pred_sims <- map_dfc(grep("^sim_[0-9]", sim_i$sim, value=T), 
                     ~glue("out/candidates/{.x}_RE_stanfit.rds") |>
                       readRDS() |>
                       make_predictions_candidate(ensTest_df, .x, re=T) |>
                       colMeans() |>
                       as_tibble() |>
                       set_names(paste0("IP_", .x))) 
pred_avgs <- map_dfc(grep("^sim_avg", sim_i$sim, value=T), 
                     ~glue("out/ensembles/{.x}_RE_stanfit.rds") |>
                       readRDS() |>
                       make_predictions_candidate(ensTest_df, .x, re=T) |>
                       colMeans() |>
                       as_tibble() |>
                       set_names(paste0("IP_", .x))) 
ensTest_df <- ensTest_df |>
  select(-matches("sim")) |>
  bind_cols(pred_sims) |>
  bind_cols(pred_avgs)
write_csv(ensTest_df, "out/ensemble_oos.csv")



# performance plot --------------------------------------------------------

# Evaluation metrics:
# - Global, all values: R2, r, RMSE, ROC-AUC, PR-AUC
# - Median(within-site values): R2, r, RMSE, ROC-AUC, PR-AUC
# - Values (among site mean): R2, r, RMSE

ensTest_df <- read_csv("out/ensemble_oos.csv") 
ensTest_df <- ensTest_df |>
  mutate(liceTreat=factor(liceTreat),
         lice_g05=factor(licePerFish_rtrt^4 > 0.5)) 

# Global
metrics_global <- full_join(
  ensTest_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sim) |>
    summarise(rsq=rsq_vec(value, truth=licePerFish_rtrt),
              r=cor(value, licePerFish_rtrt, use="pairwise"),
              rmse=rmse_vec(value, truth=licePerFish_rtrt),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup(),
  ensTest_df |>
    select(-starts_with("IP_predF")) |> rename_with(~str_replace(.x, "pr_pred", "IP_pred")) |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    filter(!is.na(value)) |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sim) |>
    summarise(ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
              PR_AUC=pr_auc_vec(value, truth=lice_g05, event_level="second")) |>
    ungroup()
)


# Median within-site
metrics_within <- full_join(
  ensTest_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sepaSite, sim) |>
    summarise(rsq=rsq_vec(value, truth=licePerFish_rtrt),
              r=cor(value, licePerFish_rtrt, use="pairwise"),
              rmse=rmse_vec(value, truth=licePerFish_rtrt),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    group_by(sim) |>
    summarise(rsq=mean(rsq, na.rm=T),
              r=mean(r, na.rm=T),
              rmse=mean(rmse, na.rm=T),
              N=mean(N, na.rm=T),
              prop_g05=mean(prop_g05),
              prop_0=mean(prop_0, na.rm=T)) |>
    ungroup(),
  ensTest_df |>
    select(-starts_with("IP_predF")) |> rename_with(~str_replace(.x, "pr_pred", "IP_pred")) |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    filter(!is.na(value)) |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sepaSite, sim) |>
    summarise(ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
              PR_AUC=pr_auc_vec(value, truth=lice_g05, event_level="second")) |>
    group_by(sim) |>
    summarise(ROC_AUC=mean(ROC_AUC, na.rm=T),
              PR_AUC=mean(PR_AUC, na.rm=T)) |>
    ungroup()
)

# Site means
metrics_among <- ensTest_df |>
  pivot_longer(starts_with("IP_"), names_to="sim") |>
  mutate(sim=str_remove(sim, "IP_")) |>
  group_by(sepaSite, sim) |>
  summarise(value=mean(value), 
            licePerFish_rtrt=mean(licePerFish_rtrt),
            prop_g05=mean(lice_g05=="TRUE"),
            prop_0=mean(licePerFish_rtrt==0)) |>
  group_by(sim) |>
  summarise(rsq=rsq_vec(value, truth=licePerFish_rtrt),
            r=cor(value, licePerFish_rtrt, use="pairwise"),
            rmse=rmse_vec(value, truth=licePerFish_rtrt),
            N=n(),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0)) |>
  ungroup()

all_metrics_df <- bind_rows(
  metrics_global |> mutate(type="global"),
  metrics_within |> mutate(type="within"),
  metrics_among |> mutate(type="among")
) |>
  pivot_longer(any_of(c("rmse", "rsq", "r", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  filter(metric %in% c("rmse", "r", "ROC_AUC", "PR_AUC")) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rmse", "rsq", "r"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "RMSE", "R^2", "r"))) |>
  left_join(sim_i) |>
  arrange(lab) |>
  mutate(type=factor(type, levels=c("global", "within", "among"))) |>
  drop_na()

ggplot(all_metrics_df, aes(type, value, colour=lab_short, shape=lab_short, size=lab_short)) + 
  geom_point() + 
  facet_wrap(~metric, nrow=1, scales="free") + 
  scale_colour_manual(values=c(1, 2, 3, 4, 5, 4, 5, "grey")) + 
  scale_shape_manual(values=c(19, 19, 19, 1, 1, 1, 1, 4)) + 
  scale_size_manual(values=c(3, 3, 3, 2, 2, 1, 1, 2))


all_metrics_labs <- all_metrics_df |>
  filter(metric=="r",
         type=="within",
         grepl("(Null|Mean|Ens)", lab_short)) |>
  arrange(lab) |>
  mutate(label=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "'2D'", "'3D'", "Null"),
         value=seq(1.01, 0.81, length.out=6)) %>%
  bind_rows(., 
            . |> 
              filter(grepl("Mean", lab_short)) |>
              mutate(lab=c("2D.1", "3D.1"),
                     lab_short=c("2D", "3D"),
                     label=c(NA, NA)))

all_metrics_df |>
  arrange(desc(lab)) |>
  ggplot() + 
  geom_point(aes(type, value, colour=lab_short, shape=lab_short, size=lab_short, alpha=lab_short)) + 
  # geom_text(data=all_metrics_df |> slice_head(n=1),
  #           aes(label=paste0("n: ", N, "\n",
  #                            round(prop_treat*100), "% > threshold\n",
  #                            round(prop_0*100), "% 0's")),
  #           x=1.35, y=0.975, size=2.5, colour="grey30") +
  # annotate("rect", xmin=5.03, xmax=6, ymin=0.775, ymax=1.1, colour="grey80", fill="white", linewidth=0.2) +
  geom_text(data=all_metrics_labs, aes(type, value, label=label, colour=lab_short),
            hjust=0, nudge_x=0.6, vjust=0.5, size=2.5, parse=T) +
  geom_point(data=all_metrics_labs, position=position_nudge(x=0.2),
             aes(type, value, colour=lab_short, shape=lab_short, size=lab_short), alpha=1) +
  # scale_x_discrete(labels=scales::parse_format(), expand=expansion(mult=c(0.1, 0.25))) +
  scale_y_continuous("Out-of-sample score (2023)", breaks=seq(0, 1, by=0.5),
                     minor_breaks=seq(0, 1, by=0.1), limits=c(0, 1), oob=scales::oob_keep) +
  facet_grid(.~metric, scales="free_x", space="free_x", labeller=label_parsed) + 
  scale_colour_manual(values=c("black", "black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=-1),
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=-1),
                               "grey50")) +
  scale_linetype_manual(values=c(1, 3, 1, 2, 2, 1, 1, 3)) +
  scale_linewidth_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.25, 0.25, 0.5)) +
  scale_shape_manual(values=c(19, 1, 19, 4, 4, 1, 1, 3)) +
  scale_size_manual(values=c(rep(1.5, 5), rep(0.5, 2), 1)) +
  scale_alpha_manual(values=c(1, 1, 1, 1, 1, 0.5, 0.5, 1)) +
  theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=0.5),
        legend.position="none")
ggsave("figs/pub/validation_metrics_2023.png", width=6, height=3.5)








# mixing proportions ------------------------------------------------------

# Mixing ensemble
out_ensMixRE <- readRDS("out/ensembles/ensMix_3D_ranef_stanfit.rds")
dat_ensMixRE <- readRDS("out/ensembles/ensMix_3D_ranef_standata.rds")

b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensMixRE$sim_names, sim_i$sim)]) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="Simulation", values_to="p") |>
  mutate(lab_short=str_sub(Simulation, 1, 2),
         Simulation=factor(Simulation, levels=levels(sim_i$lab)))
p <- ggplot(b_p_post, aes(p, Simulation, colour=lab_short, fill=lab_short)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4", 
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_y_discrete(limits=rev(sort(unique(b_p_post$Simulation)))) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Ensemble mixing weight (", italic(pi[k]), ")")),
       y="Simulation") +
  theme(legend.position="none")


sim_key <- read_csv("out/sim_2019-2023/sim_i.csv") |> 
  mutate(sim=paste0("sim_", i)) |>
  select(sim, salinityMort, eggTemp, swimSpeed) |> 
  inner_join(sim_i) |>
  rename(Simulation=lab)
swim_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, swimSpeed)) |>
  group_by(swimSpeed, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(swimSpeed) |>
  mutate(swimSpeed=factor(swimSpeed, levels=c(0.0001, 0.0004, 0.0007, 0.001), 
                          labels=paste(c("0.01", "0.04", "0.07", "0.10"), "cm/s")))
eggTemp_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, eggTemp)) |>
  group_by(eggTemp, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(eggTemp) |>
  mutate(eggTemp=factor(eggTemp, levels=c(TRUE, FALSE), 
                        labels=c("f(Temp.)", "Constant")))
salMort_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, salinityMort)) |>
  group_by(salinityMort, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(salinityMort) |>
  mutate(salinityMort=factor(salinityMort, levels=c(TRUE, FALSE), 
                             labels=c("f(Sal.)", "Constant")))

p_ls <- vector("list", 3)
p_ls[[1]] <- swim_post |>
  ggplot(aes(p, swimSpeed)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Mixing weight (", italic(pi[k]), ")")),
       y="Vertical swim speed") +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
p_ls[[2]] <- eggTemp_post |>
  ggplot(aes(p, eggTemp)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Mixing weight (", italic(pi[k]), ")")),
       y="Egg production") +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
p_ls[[3]] <- salMort_post |>
  ggplot(aes(p, salinityMort)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Ensemble mixing weight (", italic(pi[k]), ")")),
       y="Larval mortality") +
  theme(legend.position="none")

p_full <- cowplot::plot_grid(
  p,
  cowplot::plot_grid(plotlist=p_ls, ncol=1, align="v", axis="tblr", labels=c("b", "c", "d")),
  ncol=2, labels=c("a", ""))
ggsave("figs/pub/ens_mix_p.png", p_full, width=7, height=6)


# Posterior summaries
swim_post |> group_by(swimSpeed) |> sevcheck::get_intervals(p)
eggTemp_post |> group_by(eggTemp) |> sevcheck::get_intervals(p)
salMort_post |> group_by(salinityMort) |> sevcheck::get_intervals(p)

# Comparisons
swim_post |> group_by(iter) |> 
  summarise(dConst_f=first(p)-last(p)) |> 
  ungroup() |> 
  summarise(pG0=mean(dConst_f>0),
            mnDiff=mean(dConst_f),
            CI_025=quantile(dConst_f, probs=0.025),
            CI_05=quantile(dConst_f, probs=0.05),
            CI_10=quantile(dConst_f, probs=0.1),
            CI_90=quantile(dConst_f, probs=0.9),
            CI_95=quantile(dConst_f, probs=0.95),
            CI_975=quantile(dConst_f, probs=0.975))
eggTemp_post |> group_by(iter) |> 
  summarise(dConst_f=first(p)-last(p)) |> 
  ungroup() |> 
  summarise(pG0=mean(dConst_f>0),
            mnDiff=mean(dConst_f),
            CI_025=quantile(dConst_f, probs=0.025),
            CI_05=quantile(dConst_f, probs=0.05),
            CI_10=quantile(dConst_f, probs=0.1),
            CI_90=quantile(dConst_f, probs=0.9),
            CI_95=quantile(dConst_f, probs=0.95),
            CI_975=quantile(dConst_f, probs=0.975))
salMort_post |> group_by(iter) |> 
  summarise(dConst_f=first(p)-last(p)) |> 
  ungroup() |> 
  summarise(pG0=mean(dConst_f>0),
            mnDiff=mean(dConst_f),
            CI_025=quantile(dConst_f, probs=0.025),
            CI_05=quantile(dConst_f, probs=0.05),
            CI_10=quantile(dConst_f, probs=0.1),
            CI_90=quantile(dConst_f, probs=0.9),
            CI_95=quantile(dConst_f, probs=0.95),
            CI_975=quantile(dConst_f, probs=0.975))



# scatterplots ------------------------------------------------------------

ensTest_df <- read_csv("out/ensemble_oos.csv")

preds_df <- ensTest_df |>
  select(rowNum, licePerFish_rtrt, liceTreat, IP_null, IP_sim_avg2D, IP_sim_avg3D,
         IP_predMix, IP_predF1, IP_predF5) |>
  pivot_longer(starts_with("IP")) |>
  mutate(name=factor(name, 
                     levels=paste0("IP_", c("null", "sim_avg2D", "sim_avg3D", 
                                            "predF1", "predF5", "predMix")),
                     labels=c("Null", "Mean['2D']", "Mean['3D']",
                              "Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']"))) 
preds_df |>
  filter(licePerFish_rtrt > 0) |>
  ggplot(aes(value, licePerFish_rtrt)) + 
  stat_smooth(data=preds_df, linetype=2, linewidth=0.5, alpha=0.5,
              colour="cadetblue", fill="cadetblue") +
  geom_point(size=0.75, shape=1, alpha=0.1) +
  geom_jitter(data=preds_df |> filter(licePerFish_rtrt==0), 
              size=0.75, shape=1, alpha=0.1, width=0, height=0.025) + 
  facet_wrap(~name, scales="free_x", labeller=label_parsed) +
  labs(x="Prediction",
       y=expression(paste("(Mean"~~italic("L. salmonis")~~"per fish)"^0.25)))
ggsave("figs/pub/predictions_2023_scatterplot.png", width=7, height=5)




# infection pressure through time -----------------------------------------

set.seed(1003)
mod <- "3D_ranef"
out_ensMixRE <- readRDS(glue("out/ensembles/ensMix_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(dat_ensMixRE$sim_names) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="sim", values_to="p") |>
  nest(p=c(iter, p))

influx_df <- readRDS("out/sim_2019-2023/processed/connectivity_day.rds") |>
  select(sepaSite, date, sim, influx_m2) |>
  mutate(sim=paste0("sim_", sim),
         influx_m2=replace_na(influx_m2, 0))
date_seq <- sort(unique(influx_df$date))
ens_ls <- vector("list", length(date_seq))
for(i in seq_along(date_seq)) {
  ens_ls[[i]] <- influx_df |>
    filter(date==date_seq[i]) |>
    inner_join(b_p_post |> mutate(p=map(p, ~.x[iter_sub,])), by=join_by("sim")) |>
    unnest(p) |>
    mutate(wtIP=influx_m2 * p) |>
    group_by(sepaSite, date, iter) |>
    summarise(ens_IP=sum(wtIP), .groups="keep") |>
    group_by(sepaSite, date) |>
    summarise(lice_mn=mean(ens_IP),
              lice_q005=quantile(ens_IP, probs=0.005),
              lice_q025=quantile(ens_IP, probs=0.025),
              lice_q975=quantile(ens_IP, probs=0.975),
              lice_q995=quantile(ens_IP, probs=0.995),
              .groups="keep") |>
    ungroup()
  cat("Finished", as.character(date_seq[i]), "\n")
}
ens_ls |>
  reduce(bind_rows) |>
  saveRDS("out/sim_2019-2023/processed/influx_ens.rds")

influx_ens <- readRDS("out/sim_2019-2023/processed/influx_ens.rds")

thresholds <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)
thresh_cols <- c("white", viridis::turbo(length(thresholds)+1))

fig_influx <- influx_ens |>
  group_by(date) |>
  rename(lice=lice_mn) |>
  summarise(lt_t1=mean(lice==0),
            lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
            lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
            lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
            lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
            lt_t6=mean(between(lice, thresholds[5], thresholds[6])),
            lt_t7=mean(lice > thresholds[6])) |>
  ungroup() |>
  pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
  mutate(threshold=factor(threshold, 
                          labels=c("0", 
                                   paste(thresholds[1:5], "-", thresholds[2:6]),
                                   paste(">", thresholds[6]))),
         threshold_num=as.numeric(threshold)) |>
  filter(threshold_num != 1) |>
  ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
  geom_area(colour="grey30", linewidth=0.05, outline.type="both") +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous("Proportion of active farms", limits=c(0,1)) +
  scale_fill_viridis_b(expression(paste("Copepodids" %.% "m"^-2 %.% "h"^-1)),
                       option="turbo", begin=0.05,
                       breaks=c(2.5, 3.5, 4.5, 5.5, 6.5),
                       labels=c("0.0001", "0.001", "0.01", "0.1", "1")) +
  theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
        panel.grid.minor.x=element_line(colour="grey", linewidth=0.2),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=11),
        legend.position="bottom", 
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm"),
        strip.text=element_text(size=11))
ggsave("figs/pub/ens_influx_daily.png", fig_influx, width=7, height=4, dpi=400)

influx_ens |>
  group_by(date) |>
  summarise(lice_CI95width=median(lice_q975-lice_q025)) |>
  ggplot(aes(date, lice_CI95width)) + geom_point(alpha=0.25)

influx_ens |>
  # filter(sepaSite=="AAC3") |>
  mutate(week=round_date(date, "month")) |>
  group_by(sepaSite, week) |>
  summarise(across(where(is.numeric), mean)) |>
  ggplot(aes(week, lice_mn^0.25, ymin=lice_q025^0.25, ymax=lice_q975^0.25)) + 
  geom_ribbon(alpha=0.5, colour=NA) + 
  geom_line() + 
  facet_wrap(~sepaSite)

influx_ens |>
  # filter(sepaSite=="AAC3") |>
  mutate(week=round_date(date, "month")) |>
  group_by(sepaSite, week) |>
  summarise(across(where(is.numeric), mean)) |>
  mutate(year=year(week),
         week_std=ymd(paste(2023, month(week), day(week), sep="-"))) |>
  ggplot(aes(week, lice_mn^0.25, group=week)) + 
  stat_slabinterval() +
  theme_bw() 








# ensemble density calculations -------------------------------------------

mod <- "ranef"
out_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]]


# particle densities ------------------------------------------------------

library(future); library(furrr)
f <- dirf("out/sim_2019-2023/processed/weekly", "Mature")
ens_ls <- vector("list", length(f))
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_CL005=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CL995=c(0,0),
                  ens_CI95width=c(0,0),
                  ens_CI99width=c(0,0),
                  sim_sd=c(0,0))
plan(multisession, workers=30)
for(i in 1:length(f)) {
  timestep <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5)))
  ps_i <- readRDS(f[i]) |>
    select(i, all_of(dat_ensMixRE$sim_names))
  ps_i_mx <- ps_i |> select(-i) |> as.matrix()
  ps_i_mx[is.na(ps_i_mx)] <- 0
  ens_ls[[i]] <- ps_i |>
    mutate(row=row_number()) |>
    mutate(ens=future_map(row, ~apply(b_p_post[iter_sub,], 1, function(x) sum(x * ps_i_mx[.x,]))),
           ens_mn=map_dbl(ens, mean),
           ens_CL005=map_dbl(ens, ~quantile(.x, probs=0.005)),
           ens_CL025=map_dbl(ens, ~quantile(.x, probs=0.025)),
           ens_CL975=map_dbl(ens, ~quantile(.x, probs=0.975)),
           ens_CL995=map_dbl(ens, ~quantile(.x, probs=0.995)),
           ens_CI95width=ens_CL975-ens_CL025,
           ens_CI99width=ens_CL995-ens_CL005,
           sim_sd=apply(ps_i_mx, 1, sd)) |>
    select(i, starts_with("ens_"), sim_sd)
  ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ens_ls[[i]]$ens_mn)))
  ps_lims$ens_CL005 <- range(c(ps_lims$ens_CL005, range(ens_ls[[i]]$ens_CL005)))
  ps_lims$ens_CL025 <- range(c(ps_lims$ens_CL025, range(ens_ls[[i]]$ens_CL025)))
  ps_lims$ens_CL975 <- range(c(ps_lims$ens_CL975, range(ens_ls[[i]]$ens_CL975)))
  ps_lims$ens_CL995 <- range(c(ps_lims$ens_CL995, range(ens_ls[[i]]$ens_CL995)))
  ps_lims$ens_CI95width <- range(c(ps_lims$ens_CI95width, range(ens_ls[[i]]$ens_CI95width))) 
  ps_lims$sim_sd <- range(c(ps_lims$sim_sd, range(ens_ls[[i]]$sim_sd)))
  cat("Finished", as.character(timestep), "\n")
  gc()
}
plan(sequential)
saveRDS(ps_lims, "out/sim_2019-2023/processed/ps_lims.rds")

timesteps <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f, "_t_", 2)[,2], 1, -5)))
ens_df <- map2_dfr(ens_ls, timesteps, ~.x |> mutate(date=.y))
# TODO: Replace NAs with 0s
saveRDS(ens_df, "out/sim_2019-2023/processed/ens_weekly.rds")

ens_avg <- ens_df |>
  group_by(i) |>
  summarise(across(where(is.numeric), mean)) |>
  ungroup()
saveRDS(ens_avg, "out/sim_2019-2023/processed/ens_avg.rds")




# maps --------------------------------------------------------------------

# Left side: Ensemble mean(copepodid density)
# Right side: Ensemble mean(weekly CI width)
ens_avg <- readRDS("out/sim_2019-2023/processed/ens_avg.rds")

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) +
  theme(legend.position=c(0.285, 0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.4, "cm"),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        axis.title=element_blank())
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) +
  theme(legend.position="none",
        axis.title=element_blank())
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) +
  theme(legend.position="none",
        axis.title=element_blank())

mn_lims <- c(0, 1)
# mn_lims <- c(0, quantile(ens_avg$ens_mn, probs=0.995))
mn_breaks <- c(0, 0.01, 0.25, 1)

# ci_lims <- c(0, 0.06)
# ci_breaks <- c(0, 5e-7, 5e-6)
ci_lims <- c(0, 0.1)
ci_breaks <- c(0, 1e-6, 1e-4)
ci_labs <- c("0", "1e-6", "1e-4")

ens_avg <- ens_avg |>
  mutate(ens_mn=pmin(ens_mn, mn_lims[2]),
         ens_CI95width=pmin(ens_CI95width, ci_lims[2]))

ens_map <- vector("list", 6)
# Mean densities
ens_map[[1]] <- westcoms_panel + 
  geom_sf(data=ens_avg |> right_join(mesh_sf, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c("",
                       option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) +
  annotate("text", x=79000, y=547000, label="Ensemble mean", size=3) +
  annotate("text", x=79000, y=525000, 
           label=expression("cop." %.% "m"^"-2" %.% "h"^"-1"), parse=T, size=3)
ens_map[[2]] <- linnhe_panel + 
  geom_sf(data=ens_avg |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[3]] <- skye_panel + 
  geom_sf(data=ens_avg |> right_join(skye_mesh, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[4]] <- westcoms_panel + 
  geom_sf(data=ens_avg |> right_join(mesh_sf, y=_),
          aes(fill=ens_CI95width), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs) +
  annotate("text", x=79000, y=547000, label="95% CI width", size=3) +
  annotate("text", x=79000, y=525000, 
           label=expression("cop." %.% "m"^-2 %.% "h"^-1), parse=T, size=3)
ens_map[[5]] <- linnhe_panel + 
  geom_sf(data=ens_avg |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_CI95width), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs)
ens_map[[6]] <- skye_panel + 
  geom_sf(data=ens_avg |> right_join(skye_mesh, y=_),
          aes(fill=ens_CI95width), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs)

plot_grid(plotlist=ens_map, ncol=2, nrow=3, labels="auto", byrow=FALSE,
          rel_heights=c(2.1, 0.98, 1.23), rel_widths=c(1, 1)) |>
  ggsave("figs/pub/ens_map.png", plot=_, width=4.75, height=9.55, dpi=300)





# maps of r ---------------------------------------------------------------

mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
ensTest_df <- read_csv("out/ensemble_oos.csv")

r_info <- tibble(breaks=seq(-1, 1, by=0.25)) |>
  mutate(break_labs=as.character(round(breaks, 1)),
         break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
         letter=letters[row_number()],
         mdpt=(breaks + (lead(breaks)-breaks)/2)) 

p_ls <- vector("list", 3)
mods <- c("IP_predF1", "IP_predF5", "IP_predMix")
for(i in seq_along(mods)) {
  mod_lab <- as.character(filter(sim_i, grepl(str_sub(mods[i], 4, -1), sim))$lab) |>
    str_remove("Ens\\['") |> str_remove("']")
  farm_r.df <- ensTest_df |> 
    rename_with(~"predColumn", .cols=matches(mods[i])) |>
    group_by(sepaSite) |> 
    summarise(r=cor(licePerFish_rtrt, predColumn)) |> 
    filter(!is.na(r)) |> 
    mutate(letter=cut(r, breaks=r_breaks, labels=letters[1:(length(r_breaks)-1)])) |>
    left_join(site_i)
  farm_r_count <- farm_r.df |>
    count(letter) |>
    mutate(scaled=n/max(n)) |>
    full_join(r_info |> select(letter, mdpt) |> drop_na()) |>
    mutate(n=replace_na(n, 0),
           scaled=replace_na(scaled, 0)) |>
    arrange(letter)
  low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
                        y=rep(c(0, 54800/nrow(farm_r_count)), each=2) + 652000)
  x_rng <- diff(range(low_polygon$x))
  y_rng <- diff(range(low_polygon$y))
  farm_r_bar.df <- map_dfr(1:nrow(farm_r_count), 
                           ~low_polygon |> 
                             mutate(mdpt=farm_r_count$mdpt[.x],
                                    x=if_else(x==max(x), 
                                              x, 
                                              max(x)-x_rng*farm_r_count$scaled[.x]),
                                    y=y + (y_rng*(.x-1)))
  )
  farm_r_count_labs <- farm_r_bar.df |>
    group_by(mdpt) |>
    summarise(x=min(x), y=mean(y)) |>
    ungroup() |>
    left_join(farm_r_count) |>
    mutate(prop=paste0(round(n/sum(n)*100), "%"))
  
  col_lab <- expr(r:~Ens[!!mod_lab])
  p_ls[[i]] <- farm_r.df |> 
    ggplot() + 
    geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
    geom_point(aes(easting, northing, fill=r), shape=21, size=2, 
               position=position_jitter(width=2e3, height=2e3, seed=2)) +
    geom_polygon(data=farm_r_bar.df, aes(x, y, fill=mdpt, group=mdpt), colour="grey10", linewidth=0.15) +
    geom_text(data=farm_r_count_labs, aes(x, y, label=prop), 
              size=2.25, hjust=1, vjust=0.5, nudge_x=-1000) +
    colorspace::scale_fill_binned_diverging(
      name=col_lab, palette="Blue-Red 3", rev=T, 
      limits=c(-1,1), breaks=r_info$breaks, labels=r_info$break_labs,
      l1=20, l2=90, p2=2) +
    scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                       breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
    scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                       limits=c(75000, 235000), oob=scales::oob_keep) +
    theme(legend.position="inside",
          legend.position.inside=c(ifelse(i==3, 0.195, 0.185), 0.203),
          legend.background=element_blank(),
          legend.key.height=unit(0.395, "cm"),
          legend.key.width=unit(0.0, "cm"),
          legend.text=element_text(size=6),
          legend.title=element_text(size=8, vjust=1, hjust=1),
          legend.ticks=element_line(colour="grey10", linewidth=0.25),
          legend.ticks.length=unit(0.04, "cm"),
          axis.title=element_blank()) 
}
ggarrange(p_ls[[1]], p_ls[[2]], p_ls[[3]], nrow=1, common.legend=F, labels="auto") |> 
  ggsave("figs/pub/ens_farm-r_map.png", plot=_ , width=9, height=5.5, dpi=300)




ensTest_df |>
  filter(sepaSite %in% filter(farm_r.df, r < 0)$sepaSite) |>
  ggplot(aes(IP_predMix, licePerFish_rtrt)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~sepaSite)
ensTest_df |>
  ggplot(aes(IP_predMix, licePerFish_rtrt)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~sepaSite, nrow=8, scales="free")



# maps of R2 --------------------------------------------------------------

# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# 
# # R2_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
# R2_breaks <- seq(0, 1, by=0.2)
# farm_R2.df <- ensTest_df |>
#   group_by(sepaSite) |> 
#   summarise(R2=yardstick::rsq_vec(licePerFish_rtrt, IP_predMix),
#             nObs=mean(nObs)) |> 
#   filter(!is.na(R2)) |> 
#   mutate(R2_range=cut(R2, breaks=R2_breaks, labels=letters[1:(length(R2_breaks)-1)])) |>
#   left_join(site_i)
# farm_R2_count <- farm_R2.df |>
#   count(R2_range) |>
#   mutate(scaled=n/max(n),
#          R2=(R2_breaks + (lead(R2_breaks)-R2_breaks)/2)[-length(R2_breaks)])
# low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
#                       y=rep(c(0, 54800/nrow(farm_R2_count)), each=2) + 652000)
# x_rng <- diff(range(low_polygon$x))
# y_rng <- diff(range(low_polygon$y))
# farm_R2_bar.df <- map_dfr(1:nrow(farm_R2_count), 
#                           ~low_polygon |> 
#                             mutate(R2=farm_R2_count$R2[.x],
#                                    x=if_else(x==max(x), 
#                                              x, 
#                                              max(x)-x_rng*farm_R2_count$scaled[.x]),
#                                    y=y + (y_rng*(.x-1)))
# )
# farm_R2_count_labs <- farm_R2_bar.df |>
#   group_by(R2) |>
#   summarise(x=min(x), y=mean(y)) |>
#   ungroup() |>
#   left_join(farm_R2_count) |>
#   mutate(prop=paste0(round(n/sum(n)*100), "%"))
# 
# set.seed(2)
# p <- farm_R2.df |> 
#   ggplot() + 
#   geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
#   geom_jitter(aes(easting, northing, fill=R2), shape=21, size=2, width=2e3, height=2e3) +
#   geom_polygon(data=farm_R2_bar.df, aes(x, y, fill=R2, group=R2), colour="grey10", linewidth=0.15) +
#   geom_text(data=farm_R2_count_labs, aes(x, y, label=prop), 
#             size=2.25, hjust=1, vjust=0.5, nudge_x=-1000) +
#   colorspace::scale_fill_binned_sequential(name=expression(R^2~"(2023)"), palette="Oslo", rev=F,
#                                            limits=c(0,1), breaks=R2_breaks) +
#   scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
#                      breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
#   scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
#                      limits=c(75000, 235000), oob=scales::oob_keep) +
#   theme(legend.position="inside",
#         legend.position.inside=c(0.175, 0.207),
#         legend.background=element_blank(),
#         legend.key.height=unit(0.394, "cm"),
#         legend.key.width=unit(0.0, "cm"),
#         legend.text=element_text(size=6),
#         legend.title=element_text(size=8, vjust=1, hjust=1),
#         legend.ticks=element_line(colour="grey10", linewidth=0.25),
#         legend.ticks.length=unit(0.04, "cm"),
#         axis.title=element_blank())
# ggsave("figs/pub/ens_farm-R2_map_02.png", p, width=3, height=5.5, dpi=300)

# 
# ensTest_df |>
#   filter(sepaSite %in% filter(farm_R2.df, R2 > 0.5)$sepaSite) |>
#   ggplot(aes(IP_predMix, licePerFish_rtrt)) + 
#   geom_point(alpha=0.5) + 
# #   facet_wrap(~sepaSite)
# ensTest_df |>
#   ggplot(aes(IP_predMix, licePerFish_rtrt)) +
#   geom_point(alpha=0.2, size=0.5) +
#   geom_line(aes(group=sepaSite),stat="smooth", method="lm", se=F, formula=y~x, 
#             col="cadetblue", alpha=0.5) +
#   geom_line(stat="smooth", method="lm", se=F, formula=y~x, col="cadetblue", linewidth=1)


# variability among simulations -------------------------------------------

c_df <- readRDS("out/sim_2019-2023/processed/connectivity_wk.rds")
p <- vector("list", 3)
p[[1]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_rng=diff(range(influx_m2, na.rm=T))) |> 
  group_by(sepaSite) |> 
  summarise(md_rng=median(influx_rng, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=md_rng), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
p[[2]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_CV=sd(influx_m2, na.rm=T)/mean(influx_m2, na.rm=T)) |> 
  group_by(sepaSite) |> 
  summarise(mn_CV=mean(influx_CV, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=mn_CV), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
p[[3]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2, na.rm=T)) |> 
  group_by(sepaSite) |> 
  summarise(mn_sd=mean(influx_sd^0.25, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=mn_sd), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
ggpubr::ggarrange(plotlist=p, nrow=1, common.legend=F)

c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2^0.25, na.rm=T)) |>
  ggplot(aes(week, influx_sd)) + geom_line() + facet_wrap(~sepaSite)
c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2^0.25, na.rm=T)) |>
  ggplot(aes(week, influx_sd, group=sepaSite)) + geom_line(alpha=0.25)
c_df |> 
  mutate(influx_m2=influx_m2^0.25) |>
  ggplot(aes(week, influx_m2, colour=sim)) + geom_line() + facet_wrap(~sepaSite)




# vertical distributions --------------------------------------------------

set.seed(1003)
mod <- "ranef"
out_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=500)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensMixRE$sim_names, sim_i$sim)]) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="lab", values_to="p") |>
  nest(p=c(iter, p))

z_iter <- readRDS("out/temp_summary_daily_z.rds") |>
  ungroup() |>
  left_join(sim_i |> select(sim, lab)) |>
  select(lab, region, day, z, N) |>
  left_join(b_p_post |> mutate(p=map(p, ~.x[iter_sub,]))) |>
  unnest(p) |>
  mutate(wtN=N * p) |>
  group_by(region, day, z, iter) |>
  summarise(ens_N=sum(wtN)) |>
  ungroup() 
gc()
z_ens <- z_iter |>
  group_by(region, day, z) |>
  summarise(ens_mn=mean(ens_N)) |>
  group_by(region, day) |>
  mutate(across(starts_with("ens"), ~.x/sum(.x), .names="prop_{.col}"),
         region=factor(region, levels=c("WeStCOMS", "Linnhe", "Skye"),
                       labels=c("Full domain", "Loch Linnhe", "Skye")))
gc()

z_ens |>
  ggplot(aes(day, prop_ens_mn, fill=z, colour=z, group=z)) +
  geom_area(outline.type="upper", linewidth=0.2) +
  scale_y_continuous("Ensemble proportion of copepodids (daily)") +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_fill_viridis_b("Depth bin (m)", direction=-1,
                       breaks=c(1, 2, 5, 10, 20, 30, 50)-0.01,
                       labels=c(1, 2, 5, 10, 20, 30, 50)) +
  scale_colour_viridis_b("Depth bin (m)", direction=-1,
                         breaks=c(1, 2, 5, 10, 20, 30, 50)-0.01,
                         labels=c(1, 2, 5, 10, 20, 30, 50)) +
  facet_grid(region~.) +
  theme(axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.2, "cm"), 
        legend.key.width=unit(1.5, "cm"))
ggsave("figs/pub/ens_z_distribution.png", width=4.5, height=8)










# hourly animation --------------------------------------------------------

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))
site_i <- read_csv("data/farm_sites_2023.csv") |> 
  st_as_sf(coords=c("easting", "northing"), crs=27700)

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))


f <- dirf("D:/sealice_ensembling/out/sim_2023-MarMay/processed/hourly", "Mature")
lims <- readRDS("out/sim_2023-AprMay/processed/hourly_Mature_pslims.rds")
lims <- tibble(ens_mn=c(0, 1),
               ens_CI95width=c(0, 0.5))
lims_N <- tibble(ens_mn=c(0, 1),
                 ens_CI95width=c(0, 0.02))

for(i in seq_along(f)) {
  timestep <- ymd_hms("2023-04-01 00:00:00") + 
    dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5))-1)
  
  if(file.exists(glue("D:/sealice_ensembling/figs/temp/westcoms_{format(timestep, '%F_%H')}.png"))) {
    next
  }
  
  ps_i <- readRDS(f[i]) |> filter(ens_mn > 0) |>
    mutate(ens_mn=pmin(ens_mn, lims$ens_mn[2]))
  if(nrow(ps_i)==0) next
  
  # WeStCOMS
  fig_a <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_mn), colour=NA) + 
    geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
    scale_fill_viridis_c("Ensemble mean cop./m2 (4th-rt)", option="turbo", limits=lims$ens_mn) + 
    ggtitle(format(timestep, "%b-%d %H:%M")) 
  fig_b <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_CI95width), colour=NA) + 
    geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
    scale_fill_viridis_c("Ensemble 95% CI width (4th-rt)", option="turbo", limits=lims$ens_CI95width) + 
    ggtitle(format(timestep, "%b-%d %H:%M"))
  ggpubr::ggarrange(fig_a, fig_b, nrow=1, common.legend=FALSE) |>
    ggsave(filename=glue("D:/sealice_ensembling/figs/temp/westcoms_{format(timestep, '%F_%H')}.png"), 
           plot=_, width=6.5, height=8)
  fig_a <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_mn^4), colour=NA) + 
    scale_fill_viridis_c("Ensemble mean cop./m2", option="turbo", limits=lims_N$ens_mn) + 
    ggtitle(format(timestep, "%b-%d %H:%M")) 
  fig_b <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_CI95width^4), colour=NA) + 
    scale_fill_viridis_c("Ensemble 95% CI width", option="turbo", limits=lims_N$ens_CI95width) + 
    ggtitle(format(timestep, "%b-%d %H:%M"))
  ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
    ggsave(filename=glue("D:/sealice_ensembling/figs/temp/westcoms-N_{format(timestep, '%F_%H')}.png"), 
           plot=_, width=6.5, height=8)
  
  # Linnhe
  ps_linnhe <- ps_i |> filter(i %in% linnhe_mesh$i)
  if(nrow(ps_linnhe) > 0) {
    fig_a <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_mn), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2 (4th-rt)", option="turbo", limits=lims$ens_mn) + 
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_CI95width), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 95% CI width (4th-rt)", option="turbo", limits=lims$ens_CI95width) + 
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("D:/sealice_ensembling/figs/temp/linnhe_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=7, height=4.5)
    fig_a <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_mn^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2", option="turbo", limits=lims_N$ens_mn) + 
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_CI95width^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 95% CI width", option="turbo", limits=lims_N$ens_CI95width) + 
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("D:/sealice_ensembling/figs/temp/linnhe-N_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=7, height=4.5)
  }
  
  # Skye
  ps_skye <- ps_i |> filter(i %in% skye_mesh$i)
  if(nrow(ps_skye) > 0) {
    fig_a <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_mn), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2 (4th-rt)", option="turbo", limits=lims$ens_mn) + 
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_CI95width), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 95% CI width (4th-rt)", option="turbo", limits=lims$ens_CI95width) + 
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("D:/sealice_ensembling/figs/temp/skye_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=8, height=6)
    fig_a <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_mn^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2", option="turbo", limits=lims_N$ens_mn) + 
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_CI95width^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 95% CI width", option="turbo", limits=lims_N$ens_CI95width) + 
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("D:/sealice_ensembling/figs/temp/skye-N_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=8, height=6)
  }
  gc() 
}

library(av)
sets <- c("westcoms_", "linnhe_", "skye_", 
          "westcoms-N_", "linnhe-N_", "skye-N_")
for(i in sets) {
  dirf("D:/sealice_ensembling/figs/temp", glue("{i}.*png")) |>
    av_encode_video(glue("figs/hourly_anim_{i}2023-AprMay.mp4"),
                      framerate=12)   
}
