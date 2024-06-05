# Forecasting ensembles
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk

# This script tunes and validates forecasting ensembles, selecting the best 
# performing model. Optimization is performed separately for predicting the
# mean number of lice per fish and the probability that the mean lice per fish
# exceeds the Code of Good Practice treatment threshold. Optimization is also
# performed separately for predictions 1 and 5 weeks in the future.


# setup -------------------------------------------------------------------
library(tidyverse); library(glue); library(scico)
library(tidymodels); library(DALEXtra); library(butcher)
library(baguette); # bag_mars
# library(discrim); # naive_Bayes
# library(bonsai); # lightgbm
library(future)
theme_set(theme_bw() + theme(panel.grid=element_blank()))

gridSize <- 100
cores <- 15







# compile dataset ---------------------------------------------------------

ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(liceTreat=factor(liceTreat))


advance <- c(1, 5)
preds <- vector("list", length(advance))
for(i in seq_along(advance)) {
  # expanding window validation: simulate weekly forecasts starting with 2023
  # - fit model each 'week' (2019-04-01 to 2023-week) and predict the on-farm
  #   lice 1 week and 5 weeks into the future
  folds <- sliding_period(ensFull_df, index=date, period="week", skip=52*3, lookback=Inf,
                          assess_start=advance[i], assess_stop=advance[i]+1)
  
  # Ensemble formulas -------------------------------------------------------
  mod_terms <- list(
    sim=paste(grep("sim", names(ensFull_df), value=T), collapse=" + "),
    sim_date=paste(grep("(sim|yday)", names(ensFull_df), value=T), collapse=" + ")
  )
  # 1. licePerFish_rtrt ~ IP + yday
  licePerFish_form <- map(mod_terms, ~formula(paste("licePerFish_rtrt ~", 
                                                    str_remove(.x, " \\+ lag_liceTreat"))))
  # 1. liceTreat ~ IP + yday
  liceTreat_form <- map(mod_terms, ~formula(paste("liceTreat ~", 
                                                  str_remove(.x, " \\+ lag_licePerFish_rtrt"))))
  
  # Candidate ensemble models -----------------------------------------------
  # Ridge, Elastic Net, Lasso, Random Forest, XGBoost, Neural Network, CART, KNN
  licePerFish_models <- list(
    enet=linear_reg(penalty=tune(), mixture=tune()) |>
      set_engine("glmnet") |> set_mode("regression"),
    rf=rand_forest(trees=tune(), min_n=tune()) |>
      set_engine("randomForest") |> set_mode("regression"),
    # xgb=boost_tree(tree_depth=tune(), trees=tune(), learn_rate=tune(), mtry=tune(),
    #                min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("xgboost") |> set_mode("regression"),
    # lxgb=boost_tree(tree_depth=tune(), trees=tune(), learn_rate=tune(), mtry=tune(),
    #                 min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("lightgbm") |> set_mode("regression"),
    # bart=bart(trees=tune(), prior_terminal_node_coef=tune(),
    #           prior_terminal_node_expo=tune(), prior_outcome_range=tune()) |>
    #   set_engine("dbarts") |> set_mode("regression"),
    mars=bag_mars(num_terms=tune(), prod_degree=tune(), prune_method=tune()) |>
      set_engine("earth") |> set_mode("regression"),
    # svmL=svm_linear(cost=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("regression"),
    # svmP=svm_poly(cost=tune(), degree=tune(), scale_factor=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("regression"),
    svmR=svm_rbf(cost=tune(), rbf_sigma=tune(), margin=tune()) |>
      set_engine("kernlab") |> set_mode("regression")#,
    # mlp=mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
    #   set_engine("nnet") |> set_mode("regression"),
    # knn=nearest_neighbor(neighbors=tune(), weight_func=tune(), dist_power=tune()) |>
    #   set_engine("kknn") |> set_mode("regression")
  )
  
  liceTreat_models <- list(
    enet=logistic_reg(penalty=tune(), mixture=tune()) |>
      set_engine("glmnet") |> set_mode("classification"),
    rf=rand_forest(trees=tune(), min_n=tune()) |>
      set_engine("randomForest") |> set_mode("classification"),
    # xgb=boost_tree(tree_depth=tune(), trees=tune(), learn_rate=tune(), mtry=tune(),
    #                min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("xgboost") |> set_mode("classification"),
    # lxgb=boost_tree(tree_depth=tune(), trees=tune(), learn_rate=tune(), mtry=tune(),
    #                 min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("lightgbm") |> set_mode("classification"),
    # bart=bart(trees=tune(), prior_terminal_node_coef=tune(),
    #           prior_terminal_node_expo=tune(), prior_outcome_range=tune()) |>
    #   set_engine("dbarts") |> set_mode("classification"),
    mars=bag_mars(num_terms=tune(), prod_degree=tune(), prune_method=tune()) |>
      set_engine("earth") |> set_mode("classification"),
    # nbayes=naive_Bayes(smoothness=tune(), Laplace=tune()) |>
    #   set_engine("klaR") |> set_mode("classification"),
    # svmL=svm_linear(cost=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("classification"),
    # svmP=svm_poly(cost=tune(), degree=tune(), scale_factor=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("classification"),
    svmR=svm_rbf(cost=tune(), rbf_sigma=tune(), margin=tune()) |>
      set_engine("kernlab") |> set_mode("classification")#,
    # mlp=mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
    #   set_engine("nnet") |> set_mode("classification"),
    # knn=nearest_neighbor(neighbors=tune(), weight_func=tune(), dist_power=tune()) |>
    #   set_engine("kknn") |> set_mode("classification")
  )
  
  
  # Tune ensembles ----------------------------------------------------------
  # licePerFish
  plan(multisession, workers=cores)
  licePerFish_wfs <- workflow_set(
    preproc=licePerFish_form,
    models=licePerFish_models
  ) |>
    workflow_map("tune_grid", 
                 resamples=folds, 
                 grid=gridSize,
                 metrics=metric_set(rsq),
                 control=control_grid(save_pred=T, 
                                      save_workflow=T),
                 verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished licePerFish tuning, advance:", advance[i], "\n")
  autoplot(licePerFish_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/_licePerFish_ranks_{advance[i]}wk.png"), width=15, height=5)
  
  map(c("rsq"), 
      ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
        filter(.metric==.x) |>
        select(rank, .metric, mean, model, wflow_id, .config))
  map(c("rsq"), 
      ~rank_results(licePerFish_wfs[-1,], rank_metric=.x, select_best=TRUE) |>
        filter(.metric==.x) |>
        select(rank, .metric, mean, model, wflow_id, .config)) |>
    saveRDS(glue("out/_licePerFish_ranks_{advance[i]}wk.rds"))
  gc()
  
  ## Best fits
  licePerFish_best_mod <- rank_results(licePerFish_wfs, rank_metric="rsq", select_best=TRUE)
  licePerFish_best_wf <- licePerFish_wfs |> 
    extract_workflow(licePerFish_best_mod$wflow_id[1])
  licePerFish_best_results <- licePerFish_wfs |> 
    extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
    select_best(metric="rsq")
  licePerFish_final_fit <- licePerFish_best_wf |>
    finalize_workflow(licePerFish_best_results) |>
    fit(data=ensFull_df)
  licePerFish_best_preds <- licePerFish_wfs |> 
    extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
    collect_predictions() |>
    filter(.config==licePerFish_best_mod$.config[1])
  
  saveRDS(licePerFish_best_wf, glue("out/_licePerFish_best_wf_{advance[i]}wk.rds"))
  saveRDS(licePerFish_best_results, glue("out/_licePerFish_best_results_{advance[i]}wk.rds"))
  saveRDS(licePerFish_final_fit, glue("out/_licePerFish_best_fitted_{advance[i]}wk.rds"))
  saveRDS(licePerFish_best_preds, glue("out/_licePerFish_best_preds_{advance[i]}wk.rds"))
  
  
  
  # liceTreat
  plan(multisession, workers=cores)
  liceTreat_wfs <- workflow_set(
    preproc=liceTreat_form,
    models=liceTreat_models
  ) |>
    workflow_map("tune_grid", 
                 resamples=folds, 
                 grid=gridSize,
                 metrics=metric_set(roc_auc, pr_auc),
                 control=control_grid(save_pred=T, 
                                      save_workflow=T, 
                                      event_level="second"),
                 verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished liceTreat tuning, advance:", advance[i], "\n")
  autoplot(liceTreat_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/_liceTreat_ranks_{advance[i]}wk.png"), width=15, height=5)
  map(c("roc_auc", "pr_auc"), 
      ~rank_results(liceTreat_wfs, rank_metric=.x, select_best=TRUE) |>
        filter(.metric==.x) |>
        select(rank, .metric, mean, model, wflow_id, .config))
  map(c("roc_auc", "pr_auc"), 
      ~rank_results(liceTreat_wfs, rank_metric=.x, select_best=TRUE) |>
        filter(.metric==.x) |>
        select(rank, .metric, mean, model, wflow_id, .config)) |>
    saveRDS(glue("out/_liceTreat_ranks_{advance[i]}wk.rds"))
  gc()
  
  
  # Best fits
  liceTreat_best_mod <- rank_results(liceTreat_wfs, rank_metric="roc_auc", select_best=TRUE)
  liceTreat_best_wf <- liceTreat_wfs |> 
    extract_workflow(liceTreat_best_mod$wflow_id[1])
  liceTreat_best_results <- liceTreat_wfs |> 
    extract_workflow_set_result(id=liceTreat_best_mod$wflow_id[1]) |>
    select_best(metric="roc_auc")
  liceTreat_final_fit <- liceTreat_best_wf |>
    finalize_workflow(liceTreat_best_results) |>
    fit(data=ensFull_df)
  liceTreat_best_preds <- liceTreat_wfs |> 
    extract_workflow_set_result(id=liceTreat_best_mod$wflow_id[1]) |>
    collect_predictions() |>
    filter(.config==liceTreat_best_mod$.config[1])
  saveRDS(liceTreat_best_wf, glue("out/_liceTreat_best_wf_{advance[i]}wk.rds"))
  saveRDS(liceTreat_best_results, glue("out/_liceTreat_best_results_{advance[i]}wk.rds"))
  saveRDS(liceTreat_final_fit, glue("out/_liceTreat_best_fitted_{advance[i]}wk.rds"))
  saveRDS(liceTreat_best_preds, glue("out/_liceTreat_best_preds_{advance[i]}wk.rds"))
  
  # preds[[i]] <- list(licePerFish=licePerFish_best_preds |> 
  #                      rename_with(~paste0("IP_predF_", advance[i], "wk"), starts_with(".pred")) |>
  #                      group_by(.row) |> slice_tail(n=1) |> ungroup(),
  #                    liceTreat=liceTreat_best_preds |> 
  #                      select(-.pred_FALSE) |>
  #                      rename_with(~paste0("pr_treat_", advance[i], "wk"), starts_with(".pred_TRUE")) |>
  #                      group_by(.row) |> slice_tail(n=1) |> ungroup()
  # )
  gc()
}











# plots -------------------------------------------------------------------


# 
# 
# pred_df <- imap(preds, 
#                 ~bind_cols(.x$licePerFish |> select(.row, starts_with("IP")),
#                            .x$liceTreat |> select(starts_with("pr")))) |>
#   reduce(full_join)
# 
# ensTest_df <- ensFull_df |>
#   right_join(pred_df, by=join_by(rowNum==.row)) |>
#   rename_with(.fn=~glue("IP_{.x}"), .cols=starts_with("sim")) |>
#   rowwise() |>
#   mutate(IP_avg=mean(c_across(starts_with("IP_sim")))) |>
#   ungroup() |>
#   mutate(month=month(date)) |>
#   left_join(ensFull_df |>
#               filter(year < 2023) |>
#               mutate(month=month(date)) |>
#               group_by(month) |>
#               summarise(IP_null=mean(licePerFish_rtrt)) |>
#               ungroup()) |>
#   select(-month)
# 
# write_csv(ensTest_df, "out/ensemble_oos.csv")
# 
# 
# 
# # licePerFish: Performance metrics
# licePerFish_metrics <- bind_rows(
#   ensTest_df |>
#     filter(between(month(date), 2, 6)) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
#               rsq=rsq_vec(value, truth=licePerFish_rtrt),
#               r=cor(value, licePerFish_rtrt, use="pairwise"),
#               N=n(),
#               prop_treat=mean(liceTreat=="TRUE"),
#               prop_0=mean(licePerFish_rtrt==0)) |>
#     ungroup() |>
#     mutate(dates="Feb-Jun"),
#   ensTest_df |>
#     filter(!between(month(date), 2, 6)) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
#               rsq=rsq_vec(value, truth=licePerFish_rtrt),
#               r=cor(value, licePerFish_rtrt, use="pairwise"),
#               N=n(),
#               prop_treat=mean(liceTreat=="TRUE"),
#               prop_0=mean(licePerFish_rtrt==0)) |>
#     ungroup() |>
#     mutate(dates="Jul-Jan"),
#   ensTest_df |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
#               rsq=rsq_vec(value, truth=licePerFish_rtrt),
#               r=cor(value, licePerFish_rtrt, use="pairwise"),
#               N=n(),
#               prop_treat=mean(liceTreat=="TRUE"),
#               prop_0=mean(licePerFish_rtrt==0)) |>
#     ungroup() |>
#     mutate(dates="all")
# )
# 
# 
# 
# 
# # liceTreat: Performance metrics
# liceTreat_metrics <- bind_rows(
#   ensTest_df |>
#     filter(between(month(date), 2, 6)) |>
#     select(-starts_with("IP_pred")) |> rename_with(~str_replace(.x, "pr_treat", "IP_pred")) |>
#     mutate(across(starts_with("IP_sim"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     filter(!is.na(value)) |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(ROC_AUC=roc_auc_vec(value, truth=liceTreat, event_level="second"),
#               PR_AUC=pr_auc_vec(value, truth=liceTreat, event_level="second")) |>
#     ungroup() |>
#     mutate(dates="Feb-Jun"),
#   ensTest_df |>
#     filter(!between(month(date), 2, 6)) |>
#     select(-starts_with("IP_pred")) |> rename_with(~str_replace(.x, "pr_treat", "IP_pred")) |>
#     mutate(across(starts_with("IP_sim"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     filter(!is.na(value)) |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(ROC_AUC=roc_auc_vec(value, truth=liceTreat, event_level="second"),
#               PR_AUC=pr_auc_vec(value, truth=liceTreat, event_level="second")) |>
#     ungroup() |>
#     mutate(dates="Jul-Jan"),
#   ensTest_df |>
#     select(-starts_with("IP_pred")) |> rename_with(~str_replace(.x, "pr_treat", "IP_pred")) |>
#     mutate(across(starts_with("IP_sim"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     filter(!is.na(value)) |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     summarise(ROC_AUC=roc_auc_vec(value, truth=liceTreat, event_level="second"),
#               PR_AUC=pr_auc_vec(value, truth=liceTreat, event_level="second")) |>
#     ungroup() |>
#     mutate(dates="all")
# )
# 
# 
# 
# all_metrics_df <- full_join(licePerFish_metrics, liceTreat_metrics) |>
#   pivot_longer(all_of(c("rmse", "rsq", "r", "ROC_AUC", "PR_AUC")), names_to="metric") |>
#   filter(metric %in% c("rsq", "r", "ROC_AUC")) |>
#   mutate(metric=factor(metric, levels=c("rsq", "PR_AUC", "r", "ROC_AUC"),
#                        labels=c("R^2", "'PR-AUC'", "r", "'ROC-AUC'")),
#          dates=factor(dates, levels=c("Feb-Jun", "Jul-Jan", "all"),
#                       labels=c("Feb-Jun: >0.5 lice per fish (CoGP treatment threshold)", 
#                                "Jul-Jan: >1  lice per fish (CoGP treatment threshold)",
#                                "All dates")))
# all_metrics_df |>
#   ggplot() + 
#   geom_line(aes(metric, value, colour=lab_short, group=lab, linetype=lab_short, linewidth=lab_short)) +
#   geom_point(aes(metric, value, colour=lab_short, shape=lab_short, size=lab_short)) + 
#   geom_text(data=all_metrics_df |> group_by(dates) |> slice_head(n=1),
#             aes(label=paste0("n: ", N, "\n",
#                              round(prop_treat*100), "% > threshold\n",
#                              round(prop_0*100), "% 0's")),
#             x=1.4, y=0.94, size=2.5) +
#   scale_x_discrete(labels=scales::parse_format()) +
#   scale_y_continuous("Score (out-of-sample)", breaks=seq(0, 1, by=0.5), 
#                      minor_breaks=seq(0, 1, by=0.1), limits=c(0, 1)) + 
#   scale_colour_manual("Simulation type: ", 
#                       values=c("black", "grey30", 
#                                scico(2, begin=0.2, end=0.7, palette="broc", direction=-1),
#                                "grey30")) + 
#   scale_linetype_manual("Simulation type: ", values=c(1, 2, 1, 1, 3)) +
#   scale_linewidth_manual("Simulation type: ", values=c(1, 0.75, 0.1, 0.1, 0.5)) +
#   scale_shape_manual("Simulation type: ", values=c(19, 19, 1, 1, 4)) + 
#   scale_size_manual("Simulation type: ", values=c(2, 1, 1, 1, 1)) + 
#   facet_grid(.~dates, labeller=label_wrap_gen(30)) +
#   theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
#         panel.grid.minor.y=element_line(colour="grey90", linewidth=0.1),
#         axis.title.x=element_blank(),
#         axis.text.x=element_text(vjust=0.5),
#         legend.position="bottom",
#         legend.key.height=unit(0.2, "cm"),
#         legend.key.width=unit(0.2, "cm"))
# 
# ggsave("figs/validation_metrics_ensemble.png", width=7, height=4.5)
# 
# 
# all_metrics_df |> 
#   filter(lab %in% c("Null", "Ensemble")) |> 
#   arrange(dates, metric, lab) |>
#   group_by(dates, metric) |>
#   summarise(valueEns=first(value),
#             delta=first(value)-last(value),
#             delta_pct_2D=(first(value)-last(value))/last(value)*100,
#             skill=delta/(1-last(value))) |>
#   ungroup() |>
#   select(dates, metric, valueEns, skill) |>
# full_join(all_metrics_df |> 
#             filter(lab %in% c("2D.1", "Ensemble")) |> 
#             arrange(dates, metric, lab) |>
#             group_by(dates, metric) |>
#             summarise(valueEns=first(value),
#                       delta=first(value)-last(value),
#                       delta_pct_2D=(first(value)-last(value))/last(value)*100,
#                       skill=delta/(1-last(value))) |>
#             ungroup() |>
#             select(dates, metric, delta_pct_2D)) |>
#   full_join(all_metrics_df |> 
#               filter(lab %in% c("3D.1", "Ensemble")) |> 
#               arrange(dates, metric, lab) |>
#               group_by(dates, metric) |>
#               summarise(valueEns=first(value),
#                         delta=first(value)-last(value),
#                         delta_pct_3D=(first(value)-last(value))/last(value)*100) |>
#               ungroup() |>
#               select(dates, metric, delta_pct_3D))
# 
# 
# 
# 
# 
# bind_rows(
#   ensTest_df |>
#     filter(between(month(date), 2, 6)) |>
#     select(-IP_pred) |> rename(IP_pred=pr_liceTreat) |>
#     mutate(across(starts_with("IP_"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     roc_curve(value, truth=liceTreat, event_level="second") |>
#     ungroup() |>
#     mutate(dates="Feb-Jun"),
#   ensTest_df |>
#     filter(!between(month(date), 2, 6)) |>
#     select(-IP_pred) |> rename(IP_pred=pr_liceTreat) |>
#     mutate(across(starts_with("IP_"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     roc_curve(value, truth=liceTreat, event_level="second") |>
#     ungroup() |>
#     mutate(dates="Jul-Jan"),
#   ensTest_df |>
#     select(-IP_pred) |> rename(IP_pred=pr_liceTreat) |>
#     mutate(across(starts_with("IP_"), ~.x/max(.x))) |>
#     pivot_longer(starts_with("IP_"), names_to="sim") |>
#     mutate(sim=str_remove(sim, "IP_")) |>
#     left_join(sim_i) |>
#     group_by(lab, lab_short) |>
#     roc_curve(value, truth=liceTreat, event_level="second") |>
#     ungroup() |>
#     mutate(dates="all")
# ) |>
#   ggplot(aes(1-specificity, sensitivity, group=lab, 
#              colour=lab_short, linewidth=lab_short, linetype=lab_short)) + 
#   geom_abline(colour="grey90") + geom_line() + 
#   geom_text(data=liceTreat_metrics |> filter(lab=="Ensemble"), 
#             aes(label=paste0("Ens. AUC:\n", round(ROC_AUC, 3))),
#             x=0.1, y=0.9, size=3, show.legend=FALSE) +
#   scale_colour_manual("Simulation type: ", 
#                       values=c("black", "grey30", 
#                                scico(2, begin=0.2, end=0.7, palette="broc", direction=-1))) + 
#   scale_linetype_manual("Simulation type: ", values=c(1, 3, 1, 1)) +
#   scale_linewidth_manual("Simulation type: ", values=c(1, 0.5, 0.1, 0.1)) +
#   coord_equal(xlim=c(0,1), ylim=c(0,1)) +
#   facet_grid(dates~.) +
#   theme(legend.position=c(0.8, 0.1),
#         legend.background=element_blank())
# ggsave("figs/validation_ROC.png", width=4, height=9)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Explainers --------------------------------------------------------------
# 
# 
# licePerFish_explain_fit <- licePerFish_ |> extract_workflow() |> fit(data=ensTest_df)
# licePerFish_explainer <- explain_tidymodels(
#   licePerFish_explain_fit, 
#   data = as.matrix(ensTrain_df |> select(starts_with("sim"), starts_with("yday"))), 
#   y = ensTrain_df$licePerFish_rtrt,
#   label = licePerFish_best_mod
# )
# pdp_sims <- model_profile(licePerFish_explainer, N = 100, 
#                         variables = paste0("sim_0", 1:8))
# 
# as_tibble(pdp_sims$agr_profiles) %>%
#   left_join(sim_i |> select(sim, lab_short, lab2), by=join_by(`_vname_`==sim)) |>
#   ggplot(aes(`_x_`, `_yhat_`, group = lab2, colour=lab_short)) +
#   geom_line(alpha = 0.8) +
#   labs(x = "IP_adult", 
#        y = "licePerFish_hat", 
#        color = NULL)
# pdp_yday <- model_profile(licePerFish_explainer, N = 100, 
#                           variables = c("ydayCos", "ydaySin", "ydayCosXSin"))
# as_tibble(pdp_yday$agr_profiles) %>%
#   ggplot(aes(`_x_`, `_yhat_`, group = `_vname_`, colour= `_vname_`)) +
#   geom_line(alpha = 0.8) +
#   labs(x = "yday", 
#        y = "licePerFish_hat", 
#        color = NULL)
# 
# 
