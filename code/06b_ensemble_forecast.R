# Forecasting ensembles
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk

# This script tunes and validates forecasting ensembles, selecting the best 
# performing model. Optimization is performed separately for predicting the
# mean number of lice per fish and the probability that the mean lice per fish
# exceeds the Code of Good Practice treatment threshold. Optimization is also
# performed separately for predictions 1 and5 weeks in the future.


# setup -------------------------------------------------------------------
library(tidyverse); library(glue); #library(scico)
library(tidymodels); #library(DALEXtra); library(butcher)
# library(baguette); # bag_mars
# library(discrim); # naive_Bayes
# library(bonsai); # lightgbm
library(finetune)
library(future)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
options(tidymodels.dark = TRUE)

gridSize <- 10
cores <- 10







# compile dataset ---------------------------------------------------------

ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5))

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

fit_df <- ensFull_df |>
  select(rowNum, sepaSite, sepaSiteNum, productionCycleNumber, year, date,
         licePerFish_rtrt, lice_g05, contains("yday"),
         any_of(filter(sim_i, lab_short %in% c("3D"))$sim),
         any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D"))$sim)))
rm(ensFull_df); gc()


advance <- c(1, 5)
for(i in seq_along(advance)) {
  # expanding window cross-validation: simulate weekly forecasts within each year
  # - fit model each week using all non-focus years + that year up to the week
  #   and predict the on-farm lice 1 week and 5 weeks into the future
  
  
  # folds <- group_vfold_cv(fit_df, group=year)
  # fold_rowNums <- tibble(.row=fit_df$rowNum, rowNum=fit_df$rowNum)
  years <- unique(fit_df$year)
  folds_ls <- vector("list", length(years))
  for(y in seq_along(years)) {
    fit_df_y <- fit_df |>
      mutate(date_mod=if_else(year==years[y], date+dyears(12), date)) |>
      arrange(date_mod)
    folds_ls[[y]] <- sliding_period(fit_df_y, index=date_mod, period="week", lookback=Inf, complete=F,
                                    skip=max(which(year(unique(fit_df_y$date_mod)) < 2025))-advance[i],
                                    assess_start=advance[i], assess_stop=advance[i])
    rm(fit_df_y)
  }
  folds_merged <- reduce(folds_ls, bind_rows) |>
    filter(map_lgl(splits, ~length(.x$out_id)>0))
  folds <- manual_rset(folds_merged$splits, paste0("Slice", str_pad(1:nrow(folds_merged), 3, "left", "0")))
  fold_rowNums <- folds |>
    mutate(rowNum=map(splits, ~.x$data$rowNum[.x$out_id]),
           .row=map(splits, ~.x$out_id)) |>
    select(id, .row, rowNum) |>
    unnest(c(".row", "rowNum"))
  
  
  # Ensemble formulas -------------------------------------------------------
  mod_terms <- list(
    # sim=paste(grep("^sim", names(fit_df), value=T), collapse=" + "),
    # sim_date=paste(grep("(^sim|yday)", names(fit_df), value=T), collapse=" + "),
    # c_sim=paste(grep("c_sim", names(fit_df), value=T), collapse=" + "),
    c_sim_date=paste(grep("(c_sim|yday)", names(fit_df), value=T), collapse=" + ")#,
    # sim_3D=paste(grep("^sim", names(fit_df_3D), value=T), collapse=" + "),
    # sim_date_3D=paste(grep("(^sim|yday)", names(fit_df_3D), value=T), collapse=" + "),
    # c_sim_3D=paste(grep("c_sim", names(fit_df_3D), value=T), collapse=" + "),
    # c_sim_date_3D=paste(grep("(c_sim|yday)", names(fit_df_3D), value=T), collapse=" + ")
  )
  # 1. licePerFish_rtrt ~ IP + yday
  licePerFish_form <- map(mod_terms, ~formula(paste("licePerFish_rtrt ~", .x)))
  # 1. lice_g05 ~ IP + yday
  liceBinary_form <- map(mod_terms, ~formula(paste("lice_g05 ~", .x)))
  
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
    # lxgb=boost_tree(tree_depth=tune(), trees=1000, stop_iter=tune(), learn_rate=tune(), mtry=tune(),
    #                 min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("lightgbm", num_leaves=tune()) |> set_mode("regression"),
    # bart=bart(trees=tune(), prior_terminal_node_coef=tune(),
    #           prior_terminal_node_expo=tune(), prior_outcome_range=tune()) |>
    #   set_engine("dbarts") |> set_mode("regression"),
    # mars=bag_mars(num_terms=tune(), prod_degree=tune(), prune_method=tune()) |>
    #   set_engine("earth") |> set_mode("regression")#,
    # svmL=svm_linear(cost=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("regression"),
    # svmP=svm_poly(cost=tune(), degree=tune(), scale_factor=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("regression"),
    # svmR=svm_rbf(cost=tune(), rbf_sigma=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("regression"),
    mlp=mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
      set_engine("nnet") |> set_mode("regression")#,
    # knn=nearest_neighbor(neighbors=tune(), weight_func=tune(), dist_power=tune()) |>
    #   set_engine("kknn") |> set_mode("regression")
  )
  
  liceBinary_models <- list(
    enet=logistic_reg(penalty=tune(), mixture=tune()) |>
      set_engine("glmnet") |> set_mode("classification"),
    rf=rand_forest(trees=tune(), min_n=tune()) |>
      set_engine("randomForest") |> set_mode("classification"),
    # xgb=boost_tree(tree_depth=tune(), trees=tune(), learn_rate=tune(), mtry=tune(),
    #                min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("xgboost") |> set_mode("classification"),
    # lxgb=boost_tree(tree_depth=tune(), trees=1000, stop_iter=tune(), learn_rate=tune(), mtry=tune(),
    #                 min_n=tune(), loss_reduction=tune(), sample_size=tune()) |>
    #   set_engine("lightgbm", num_leaves=tune()) |> set_mode("classification"),
    # bart=bart(trees=tune(), prior_terminal_node_coef=tune(),
    #           prior_terminal_node_expo=tune(), prior_outcome_range=tune()) |>
    #   set_engine("dbarts") |> set_mode("classification"),
    # mars=bag_mars(num_terms=tune(), prod_degree=tune(), prune_method=tune()) |>
    #   set_engine("earth") |> set_mode("classification"),
    # nbayes=naive_Bayes(smoothness=tune(), Laplace=tune()) |>
    #   set_engine("klaR") |> set_mode("classification"),
    # svmL=svm_linear(cost=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("classification"),
    # svmP=svm_poly(cost=tune(), degree=tune(), scale_factor=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("classification"),
    # svmR=svm_rbf(cost=tune(), rbf_sigma=tune(), margin=tune()) |>
    #   set_engine("kernlab") |> set_mode("classification"),
    mlp=mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
      set_engine("nnet") |> set_mode("classification")#,
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
    workflow_map(#"tune_grid",
      "tune_race_anova", 
      resamples=folds, 
      grid=gridSize,
      metrics=metric_set(rmse),
      control=control_race(save_pred=T, 
                           save_workflow=T,
                           parallel_over="everything",
                           verbose=F,
                           verbose_elim=T,
                           burn_in=30),
      # control=control_grid(save_pred=T,
      #                      save_workflow=T,
      #                      parallel_over="everything"),
      verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished licePerFish tuning, advance:", advance[i], "\n")
  autoplot(licePerFish_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/licePerFish_ranks_{advance[i]}wk.png"), width=15, height=5)
  
  for(m in c("rmse")) {
    map(m, 
        ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
          filter(.metric==.x) |>
          select(rank, .metric, mean, model, wflow_id, .config))
    map(m,
        ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
          filter(.metric==.x) |>
          select(rank, .metric, mean, model, wflow_id, .config)) |>
      saveRDS(glue("out/ensembles/licePerFish_ranks_{advance[i]}wk_{m}.rds"))
    gc()
    
    ## Best fits
    licePerFish_best_mod <- rank_results(licePerFish_wfs, rank_metric=m, select_best=TRUE)
    licePerFish_best_wf <- licePerFish_wfs |> 
      extract_workflow(licePerFish_best_mod$wflow_id[1])
    licePerFish_best_results <- licePerFish_wfs |> 
      extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
      select_best(metric=m)
    licePerFish_final_fit <- licePerFish_best_wf |>
      finalize_workflow(licePerFish_best_results) |>
      fit(data=fit_df)
    licePerFish_best_preds <- licePerFish_wfs |> 
      extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
      collect_predictions() |>
      filter(.config==licePerFish_best_mod$.config[1]) |>
      left_join(fold_rowNums) 
    
    saveRDS(licePerFish_best_wf, glue("out/ensembles/licePerFish_best_wf_{advance[i]}wk_{m}.rds"))
    saveRDS(licePerFish_best_results, glue("out/ensembles/licePerFish_best_results_{advance[i]}wk_{m}.rds"))
    saveRDS(licePerFish_final_fit, glue("out/ensembles/licePerFish_best_fitted_{advance[i]}wk_{m}.rds"))
    write_csv(licePerFish_best_preds, glue("out/ensembles/CV_ensFc-{advance[i]}_{m}.csv")) 
  }
  
  
  
  # liceBinary
  plan(multisession, workers=cores)
  liceBinary_wfs <- workflow_set(
    preproc=liceBinary_form,
    models=liceBinary_models
  ) |>
    workflow_map(#"tune_grid",
      "tune_race_anova",
      resamples=folds,
      grid=gridSize,
      metrics=metric_set(roc_auc),
      control=control_race(save_pred=T, 
                           save_workflow=T,
                           parallel_over="everything",
                           verbose=F,
                           verbose_elim=T,
                           burn_in=30,
                           event_level="second"),
      # control=control_grid(save_pred=T,
      #                      save_workflow=T,
      #                      parallel_over="everything",
      #                      event_level="second"),
      verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished liceBinary tuning, advance:", advance[i], "\n")
  autoplot(liceBinary_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/liceBinary_ranks_{advance[i]}wk.png"), width=15, height=5)
  
  for(m in c("roc_auc")) {
    
    # Best fits
    liceBinary_best_mod <- collect_metrics(liceBinary_wfs) |>
      filter(.metric==m) |>
      arrange(desc(mean)) |>
      slice_max(mean, n=1, with_ties=F) |>
      select(.metric, mean, model, wflow_id, .config)
    liceBinary_best_wf <- liceBinary_wfs |>
      extract_workflow(liceBinary_best_mod$wflow_id[1])
    liceBinary_best_results <- liceBinary_wfs |>
      extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
      select_best(metric=m)
    liceBinary_final_fit <- liceBinary_best_wf |>
      finalize_workflow(liceBinary_best_results) |>
      fit(data=fit_df)
    liceBinary_best_preds <- liceBinary_wfs |>
      extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
      collect_predictions() |>
      filter(.config==liceBinary_best_mod$.config[1]) |>
      left_join(fold_rowNums)
    
    saveRDS(liceBinary_best_mod, glue("out/ensembles/liceBinary_ranks_{advance[i]}wk_{m}.rds"))
    saveRDS(liceBinary_best_wf, glue("out/ensembles/liceBinary_best_wf_{advance[i]}wk_{m}.rds"))
    saveRDS(liceBinary_best_results, glue("out/ensembles/liceBinary_best_results_{advance[i]}wk_{m}.rds"))
    saveRDS(liceBinary_final_fit, glue("out/ensembles/liceBinary_best_fitted_{advance[i]}wk_{m}.rds"))
    write_csv(liceBinary_best_preds, glue("out/ensembles/CV_ensFc-{advance[i]}_{m}.csv"))
  }
  
  gc()
}




