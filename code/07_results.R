# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk




# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(ggpubr)
library(cowplot)
library(scico)
library(colorspace)
library(ggdist)
library(rstan)
library(tidymodels)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R")

cmr <- readRDS("../../00_misc/cmr_cmaps.RDS")

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(liceTreat=factor(liceTreat))

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
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predF1", "predF5", "predMRE", "sim_avg2D", "sim_avg3D", "null"),
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
  mutate(week=floor(week(date)/2)) |>
  left_join(ensFull_df |>
              filter(year < 2023) |>
              mutate(week=floor(week(date)/2)) |>
              group_by(week) |>
              summarise(IP_null=mean(licePerFish_rtrt)) |>
              ungroup()) |>
  select(-week)

# Mixing ensemble
# out_ensMix <- readRDS("out/ensMix_stanfit.rds")
# out_ensMix3D <- readRDS("out/ensMix_3D_stanfit.rds")
# out_ensMixRE <- readRDS("out/ensMix_pRE_huRE_stanfit.rds")
out_ensMixRE3D <- readRDS("out/ensMix_3D_pRE_huRE_stanfit.rds")
# out_ensMixRER2D2 <- readRDS("out/ensMix_pRE_huRE-R2D2_stanfit.rds")
# out_ensMixRER2D23D <- readRDS("out/ensMix_3D_pRE_huRE-R2D2_stanfit.rds")

# Forecasting ensemble validation predictions 
out_ensFcst <- full_join(
  dirf("out", "^_lice.*_best_preds_1wk") |>
    map(~readRDS(.x) |> 
          select(.row, any_of(c(".pred", ".pred_TRUE"))) |>
          group_by(.row) |> slice_tail(n=1) |> ungroup()) |>
    reduce(full_join) |>
    rename(IP_predF1=.pred, pr_predF1=.pred_TRUE),
  dirf("out", "^_lice.*_best_preds_5wk") |>
    map(~readRDS(.x) |> 
          select(.row, any_of(c(".pred", ".pred_TRUE"))) |>
          group_by(.row) |> slice_tail(n=1) |> ungroup()) |>
    reduce(full_join) |>
    rename(IP_predF5=.pred, pr_predF5=.pred_TRUE)) |>
  rename(rowNum=.row)
  
# Merge predictions
ensTest_df <- ensTest_df |>
  # mutate(IP_predM=make_predictions_ensMix(out_ensMix, 
  #                                         newdata=ensTest_df |> 
  #                                           select(-contains("sim_avg")), 
  #                                         re=FALSE) |> colMeans()) |>
  # mutate(IP_predM3D=make_predictions_ensMix(out_ensMix3D,
  #                                           newdata=ensTest_df |>
  #                                             select(-any_of(filter(sim_i, lab_short=="2D")$sim),
  #                                                    -contains("sim_avg")),
  #                                           re=FALSE) |> colMeans()) |>
  # mutate(IP_predMRE=make_predictions_ensMix(out_ensMixRE, 
  #                                           newdata=ensTest_df |>
  #                                             select(-contains("sim_avg")),
  #                                           re=TRUE) |> colMeans()) |>
  mutate(IP_predMRE=make_predictions_ensMix(out_ensMixRE3D, 
                                          newdata=ensTest_df |> 
                                            select(-any_of(filter(sim_i, lab_short=="2D")$sim),
                                                   -contains("sim_avg")),
                                              re=TRUE) |> colMeans()) |>
  # mutate(IP_predMRER2D2=make_predictions_ensMix(out_ensMixRER2D2, 
  #                                                 newdata=ensTest_df |> 
  #                                                   select(-contains("sim_avg")),
  #                                                 re=TRUE) |> colMeans()) |>
  # mutate(IP_predMRER2D23D=make_predictions_ensMix(out_ensMixRER2D23D, 
  #                                                 newdata=ensTest_df |> 
  #                                                   select(-any_of(filter(sim_i, lab_short=="2D")$sim),
  #                                                          -contains("sim_avg")),
  #                                                 re=TRUE) |> colMeans()) |>
  left_join(out_ensFcst)

ensTest_df |>
  summarise(across(starts_with("IP_pred"), ~cor(licePerFish_rtrt, .x, use="pairwise"))) |>
  pivot_longer(everything()) |>
  ggplot(aes(name, value)) + geom_point()

pred_sims <- map_dfc(grep("^sim", sim_i$sim, value=T), 
                     ~make_predictions_candidate(readRDS(glue("out/{.x}_RE_stanfit.rds")),
                                                 ensTest_df, .x, re=T) |>
                       colMeans() %>%
                       as_tibble() |>
                       set_names(paste0("IP_", .x))) 
ensTest_df <- ensTest_df |>
  select(-matches("sim")) |>
  bind_cols(pred_sims)
write_csv(ensTest_df, "out/_ensemble_oos.csv")



# performance plot --------------------------------------------------------

licePerFish_metrics <- ensTest_df |>
  pivot_longer(starts_with("IP_"), names_to="sim") |>
  mutate(sim=str_remove(sim, "IP_")) |>
  left_join(sim_i) |>
  group_by(lab, lab_short) |>
  summarise(rsq=rsq_vec(value, truth=licePerFish_rtrt),
            r=cor(value, licePerFish_rtrt, use="pairwise"),
            rmse=mae_vec(value, truth=licePerFish_rtrt),
            N=n(),
            prop_treat=mean(liceTreat=="TRUE"),
            prop_0=mean(licePerFish_rtrt==0)) |>
  ungroup()
liceTreat_metrics <- ensTest_df |>
  select(-starts_with("IP_predF")) |> rename_with(~str_replace(.x, "pr_pred", "IP_pred")) |>
  mutate(across(starts_with("IP_sim"), ~.x/max(.x))) |>
  pivot_longer(starts_with("IP_"), names_to="sim") |>
  filter(!is.na(value)) |>
  mutate(sim=str_remove(sim, "IP_")) |>
  left_join(sim_i) |>
  group_by(lab, lab_short) |>
  summarise(ROC_AUC=roc_auc_vec(value, truth=liceTreat, event_level="second"),
            PR_AUC=pr_auc_vec(value, truth=liceTreat, event_level="second")) |>
  ungroup()


all_metrics_df <- full_join(licePerFish_metrics, liceTreat_metrics) |>
  pivot_longer(any_of(c("rmse", "rsq", "r", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  filter(metric %in% c("rsq", "r", "ROC_AUC", "PR_AUC")) |>
  mutate(metric=factor(metric, levels=c("PR_AUC", "ROC_AUC", "rmse", "rsq", "r"),
                       labels=c("'PR-AUC'", "'ROC-AUC'", "rmse", "R^2", "r"))) 
all_metrics_labs <- all_metrics_df |>
  filter(metric=="r",
         grepl("(Null|Mean|Ens)", lab_short)) |>
  arrange(lab) |>
  mutate(label=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['Mix']", "'2D Mean'", "'3D'", "Null"))
all_metrics_df |>
  ggplot() + 
  geom_line(aes(metric, value, colour=lab_short, group=lab,
                linetype=lab_short, linewidth=lab_short, alpha=lab_short)) +
  geom_point(aes(metric, value, colour=lab_short, shape=lab_short, size=lab_short)) + 
  geom_text(data=all_metrics_df |>
              slice_head(n=1),
            aes(label=paste0("n: ", N, "\n",
                             round(prop_treat*100), "% > threshold\n",
                             round(prop_0*100), "% 0's")),
            x=1.35, y=0.925, size=2.5) +
  geom_text(data=all_metrics_labs, aes(metric, value, label=label, colour=lab_short),
            hjust=0, nudge_x=0.1, vjust=0, size=3, parse=T) +
  scale_x_discrete(labels=scales::parse_format(), expand=expansion(mult=c(0.1, 0.25))) +
  scale_y_continuous("Out-of-sample score (2023)", breaks=seq(0, 1, by=0.5),
                     minor_breaks=seq(0, 1, by=0.1), limits=c(0, 1)) +
  scale_colour_manual(values=c("black", "black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=-1),
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=-1),
                               "grey50")) +
  scale_linetype_manual(values=c(1, 3, 1, 2, 2, 1, 1, 3)) +
  scale_linewidth_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.25, 0.25, 0.5)) +
  scale_shape_manual(values=c(19, 1, 19, 4, 4, 1, 1, 3)) +
  scale_size_manual(values=c(rep(2.5, 5), rep(0.75, 3))) +
  scale_alpha_manual(values=c(1, 1, 1, 1, 1, 0.5, 0.5, 1)) +
  theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=0.5),
        legend.position="none")# + 
  # facet_wrap(~metric, scales="free", nrow=1)
ggsave("figs/pub/validation_metrics_2023.png", width=3.5, height=3.5)




# mixing proportions ------------------------------------------------------

# Mixing ensemble
mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensMixRE$sim_names, sim_i$sim)]) |>
  pivot_longer(everything(), names_to="Simulation", values_to="p") |>
  mutate(lab_short=str_sub(Simulation, 1, 2),
         Simulation=factor(Simulation, levels=levels(sim_i$lab)))
p <- ggplot(b_p_post, aes(p, Simulation, colour=lab_short, fill=lab_short)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, 
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1) +
  scale_colour_manual(values=scico(2, begin=0.2, end=0.7, palette="broc", direction=-1)[2]) +
  scale_fill_manual(values=scico(2, begin=0.2, end=0.7, palette="broc", direction=-1)[2]) +
  scale_y_discrete(limits=rev(sort(unique(b_p_post$Simulation)))) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  labs(x=expression(paste("Mixing weight (", italic(pi[k]), ")"))) +
  theme(legend.position="none")
p
ggsave("figs/pub/ens_mix_p.png", p, width=5, height=4)








# scatterplots ------------------------------------------------------------

ensTest_df <- read_csv("out/_ensemble_oos.csv")

preds_df <- ensTest_df |>
  select(rowNum, licePerFish_rtrt, liceTreat, IP_null, IP_sim_avg2D, IP_sim_avg3D,
         IP_predMRE3D, IP_predF1, IP_predF5) |>
  pivot_longer(starts_with("IP")) |>
  mutate(name=factor(name, 
                     levels=paste0("IP_", c("null", "sim_avg2D", "sim_avg3D", 
                                            "predF1", "predF5", "predMRE3D")),
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
mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=10)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensMixRE$sim_names, sim_i$sim)]) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="lab", values_to="p") |>
  nest(p=c(iter, p))
c_iter <- dirrf("out/biotracker", "connectivity_day.rds") |>
  map_dfr(~readRDS(.x) |> select(sepaSite, date, sim, influx_m2) |>
            mutate(path=.x)) |>
  mutate(sim=paste0("sim_", sim),
         influx_m2=replace_na(influx_m2, 0)) |>
  select(-path) |>
  left_join(sim_i |> select(sim, lab)) |>
  select(-sim) |>
  inner_join(b_p_post |> mutate(p=map(p, ~.x[iter_sub,]))) |>
  unnest(p) |>
  mutate(wtIP=influx_m2 * p) |>
  group_by(sepaSite, date, iter) |>
  summarise(ens_IP=sum(wtIP)) |>
  ungroup() 
gc()
c_ens <- c_iter |>
  group_by(sepaSite, date) |>
  summarise(lice=mean(ens_IP),
            lice_q05=quantile(ens_IP, probs=0.05),
            lice_q95=quantile(ens_IP, probs=0.95)) |>
  ungroup()

thresholds <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)
thresh_cols <- c("white", viridis::turbo(length(thresholds)+1))

fig_connect <- c_ens |>
  group_by(date) |>
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
ggsave("figs/pub/ens_influx_daily_.png", fig_connect, width=7, height=4, dpi=400)

rm(c_iter); rm(c_ens); gc()







# ensemble density calculations -------------------------------------------

mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]]


# particle densities ------------------------------------------------------

library(future); library(furrr)
f <- dirf("out/sim_2019-2023/processed/weekly", "Mature")
ens_ls <- vector("list", length(f))
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CI95width=c(0,0))
plan(multisession, workers=40)
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
           ens_CL025=map_dbl(ens, ~quantile(.x, probs=0.025)),
           ens_CL975=map_dbl(ens, ~quantile(.x, probs=0.975)),
           ens_CI95width=ens_CL975-ens_CL025) |>
    select(i, starts_with("ens_"))
  ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ens_ls[[i]]$ens_mn)))
  ps_lims$ens_CL025 <- range(c(ps_lims$ens_CL025, range(ens_ls[[i]]$ens_CL025)))
  ps_lims$ens_CL975 <- range(c(ps_lims$ens_CL975, range(ens_ls[[i]]$ens_CL975)))
  ps_lims$ens_CI95width <- range(c(ps_lims$ens_CI95width, range(ens_ls[[i]]$ens_CI95width))) 
  cat("Finished", as.character(timestep), "\n")
  gc()
}
plan(sequential)



# maps --------------------------------------------------------------------

# Left side: Ensemble mean(copepodid density)
# Right side: Ensemble mean(weekly CI width)
set.seed(1003)
mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=5)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensMixRE$sim_names, sim_i$sim)]) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="lab", values_to="p") |>
  nest(p=c(iter, p))
ps_iter <- dirrf("out/biotracker", "psteps_avg_rtrt") |>
  map_dfr(~readRDS(.x) |> st_drop_geometry() |> as_tibble() |> 
            select(i, sim, mean_N) |>
            mutate(path=.x)) |>
  mutate(mean_N=replace_na(mean_N, 0)) |>
  select(-path) |>
  left_join(sim_i) |>
  select(-sim) |>
  inner_join(b_p_post |> mutate(p=map(p, ~.x[iter_sub,]))) |>
  unnest(p) |>
  mutate(wtN=mean_N * p) |>
  group_by(i, iter) |>
  summarise(ens_N=sum(wtN)) |>
  ungroup() 
gc()
ps_ens <- ps_iter |>
  mutate(ens_N_nat=ens_N^4) |>
  group_by(i) |>
  summarise(ens_mn=mean(ens_N),
            ens_sd=sd(ens_N),
            ens_q025=quantile(ens_N, probs=0.025),
            ens_q975=quantile(ens_N, probs=0.975),
            ensN_mn=mean(ens_N_nat),
            ensN_sd=sd(ens_N_nat),
            ensN_q025=quantile(ens_N_nat, probs=0.025),
            ensN_q975=quantile(ens_N_nat, probs=0.975)) |>
  ungroup() |>
  mutate(ens_95CIwidth=ens_q975-ens_q025,
         ens_CV=ens_sd/ens_mn,
         ensN_95CIwidth=ensN_q975-ensN_q025,
         ensN_CV=ensN_sd/ensN_mn)
rm(ps_iter); gc()

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))
fyne_mesh <- mesh_sf |> 
  st_crop(c(xmin=182000, xmax=245000, ymin=630000, ymax=714000))

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) +
  theme(legend.position=c(0.285, 0.13),
        legend.background=element_blank(),
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.42, "cm"),
        legend.title=element_text(size=8, hjust=0.5),
        legend.text=element_text(size=7))
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) +
  theme(legend.position="none")
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) +
  theme(legend.position="none")

mn_lims <- c(0, max(ps_ens$ens_mn))
mn_breaks <- c(0, 0.01, 0.25, 1)

ci_lims <- c(0, 0.05)
ci_breaks <- c(0, 1e-7, 2e-6)

ps_ens <- ps_ens |>
  mutate(ens_95CIwidth=pmin(ens_95CIwidth, ci_lims[2]))

ens_map <- vector("list", 6)
# Mean densities
ens_map[[1]] <- westcoms_panel + 
  geom_sf(data=ps_ens |> right_join(mesh_sf, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(expression(paste(atop("Ensemble mean", "cop." %.% "m"^-2 %.% "h"^-1))),
                       option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[2]] <- linnhe_panel + 
  geom_sf(data=ps_ens |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[3]] <- skye_panel + 
  geom_sf(data=ps_ens |> right_join(skye_mesh, y=_),
          aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[4]] <- westcoms_panel + 
  geom_sf(data=ps_ens |> right_join(mesh_sf, y=_),
          aes(fill=ens_95CIwidth), colour=NA) + 
  scale_fill_scico(expression(paste(atop("95% CI width", "cop." %.% "m"^-2 %.% "h"^-1))),
                   palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_breaks)
ens_map[[5]] <- linnhe_panel + 
  geom_sf(data=ps_ens |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_95CIwidth), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_breaks)
ens_map[[6]] <- skye_panel + 
  geom_sf(data=ps_ens |> right_join(skye_mesh, y=_),
          aes(fill=ens_95CIwidth), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_breaks)

plot_grid(plotlist=ens_map, ncol=2, nrow=3, labels="auto", byrow=FALSE,
          rel_heights=c(2.1, 0.98, 1.23), rel_widths=c(1, 1)) |>
  ggsave("figs/pub/ens_map_.png", plot=_, width=4.75, height=9.55, dpi=300)




skye_panel + 
  geom_sf(data=ps_ens |> right_join(skye_mesh, y=_),
          aes(fill=ensN_CV), colour=NA) + 
  scale_fill_viridis_c(option="turbo")




# vertical distributions --------------------------------------------------

set.seed(1003)
mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

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
fyne_mesh <- mesh_sf |> 
  st_crop(c(xmin=182000, xmax=245000, ymin=630000, ymax=714000))

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) +
  theme(legend.position=c(0.285, 0.13),
        legend.background=element_blank(),
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.42, "cm"),
        legend.title=element_text(size=8, hjust=0.5),
        legend.text=element_text(size=7))
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) +
  theme(legend.position="none")
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) +
  theme(legend.position="none")


f <- dirf("out/sim_2023-AprMay/processed/hourly", "Mature")
lims <- readRDS("out/sim_2023-AprMay/processed/hourly/ps_lims.rds")

for(i in seq_along(f)) {
  timestep <- ymd_hms("2023-04-01 00:00:00") + 
    dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5))-1)
  ps_i <- readRDS(f[i]) |>
    filter(ens_mn > 0)
  fig_a <- mesh_sf |> 
    inner_join(ps_i) |>
    ggplot() + 
    geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
    geom_sf(aes(fill=ens_mn), colour=NA) + scale_fill_viridis_c(option="turbo") + 
    ggtitle(format(timestep, "%b-%d %H:%M"))
  fig_b <- mesh_sf |> 
    inner_join(ps_i) |>
    ggplot() + 
    geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
    geom_sf(aes(fill=ens_CI95width), colour=NA) + scale_fill_viridis_c(option="turbo") + 
    ggtitle(format(timestep, "%b-%d %H:%M"))
}
