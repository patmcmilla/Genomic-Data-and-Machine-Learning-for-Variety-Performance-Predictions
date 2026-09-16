## Bayesian hyperparameter optimization for the gBART model.
##
## Tunes the number of trees, the regularization strength, the two tree structure priors and the
## residual variance prior for the gBART specification reported in the manuscript: the
## environmental factor together with the SNP markers, with no year term. The search space is the
## one given in supplementary section S3.
##
## Tuning uses trial years 2004 to 2008 only, so the 2009 to 2018 evaluation years of the main
## analysis are never seen here. Candidates are scored on held-out trial years, by the mean over
## the (year, location) cells of the held-out year of the correlation between observed and
## predicted yield, averaged over the folds.

library(tidyverse)
options(java.parameters = c("-Xmx32g", "--add-modules=jdk.incubator.vector")) ## Increase memory available to BART machine.
library(bartMachine)
library(rBayesianOptimization)

## The grouping the model is tuned at. The main analysis runs both 'location' and 'area'.
Grouping <- "location"

## Load in required datasets ####
Yield.dat <- read_csv(file = "Data/Winter Wheat Trial Data.csv") %>%
  filter(category == "ww") %>%
  select(-category)

Genetic.dat <- read_csv(file = "Data/occc_gs_updated_8k.csv") %>%
  select(-"X1", -"plate_id", -"y") %>%
  rename_with(~gsub("_.*", "", .x), -variety) %>%
  rename_with(~gsub("X", "V", .x), -variety)

Merged.dat <- Yield.dat %>%
  filter(area != 5, # Too few observations
         location != 'mc') %>% # Too few observations
  dplyr::select("yield", "variety", "location", "year", "area") %>%
  right_join(., Genetic.dat, by = "variety") %>%
  drop_na(yield) %>%
  mutate(location = as.factor(location),
         area = as.factor(area),
         variety = as.factor(variety)) %>%
  filter(year >= 2004)

## The SNP columns, which is everything in the genotype file except the variety name.
Marker.cols <- names(Genetic.dat)[-1]

## Tuning window only. The 2009 to 2018 evaluation years never enter this script.
Tuning.dat <- Merged.dat %>%
  filter(year < 2009)

## Tuning folds: validate on each year from 2004 to 2008, training on the tuning years that precede it. 
Folds <- map(2005:2008, function(Val.Year){
  list(Val.Year = Val.Year,
       Train = filter(Tuning.dat, year < Val.Year),
       Test = filter(Tuning.dat, year == Val.Year))
})

Make.X <- function(dat){
  dat %>%
    select(all_of(Grouping), all_of(Marker.cols)) %>%
    as.matrix() %>%
    as.data.frame()
}

Score.Fold <- function(Fold, Cfg){
  Mod <- bartMachine(X = Make.X(Fold$Train),
                     y = Fold$Train$yield,
                     num_trees = Cfg$num_trees,
                     k = Cfg$k,
                     alpha = Cfg$alpha,
                     beta = Cfg$beta,
                     nu = Cfg$nu,
                     q = Cfg$q,
                     mem_cache_for_speed = T, # To remove the speed for memory trade-off. It doesn't work for large P small N problems such as this
                     run_in_sample = F) # Don't compute in sample statistics to help improve running time

  Fold$Test %>%
    mutate(Prediction = predict(Mod, Make.X(Fold$Test))) %>%
    group_by(year, location) %>%
    summarise(Correlation = cor(yield, Prediction), .groups = "drop") %>%
    pull(Correlation) %>%
    mean(na.rm = TRUE)
}

## Search space ####
Grid <- list(k = c(1, 2, 3, 5),
             alpha = c(0.75, 0.95),
             beta = c(1, 2, 4),
             ## The three sigma prior settings named by Chipman et al., kept paired not crossed.
             nu_q = list(c(3, 0.90), c(3, 0.99), c(10, 0.75)))

Snap <- function(x, n) max(1L, min(n, as.integer(round(x))))

gBART.HypOpt.fun <- function(num.trees, i_k, i_alpha, i_beta, i_nuq){
  nq <- Grid$nu_q[[Snap(i_nuq, length(Grid$nu_q))]]
  Cfg <- list(num_trees = round(num.trees),
              k = Grid$k[Snap(i_k, length(Grid$k))],
              alpha = Grid$alpha[Snap(i_alpha, length(Grid$alpha))],
              beta = Grid$beta[Snap(i_beta, length(Grid$beta))],
              nu = nq[1], q = nq[2])

  Scores <- tryCatch(map_dbl(Folds, Score.Fold, Cfg = Cfg),
                     error = function(e){
                       message("  fit failed: ", conditionMessage(e))
                       NA_real_
                     })

  list(Score = if (all(is.na(Scores))) 0 else mean(Scores, na.rm = TRUE),
       Pred = 1)
}

## Hyperparameter optimization ####
Search.bounds <- list(num.trees = c(25L, 200L),
                      i_k = c(0.51, length(Grid$k) + 0.49),
                      i_alpha = c(0.51, length(Grid$alpha) + 0.49),
                      i_beta = c(0.51, length(Grid$beta) + 0.49),
                      i_nuq = c(0.51, length(Grid$nu_q) + 0.49))

Opt.result <- BayesianOptimization(FUN = gBART.HypOpt.fun,
                                   bounds = Search.bounds,
                                   init_points = 5,
                                   n_iter = 50,
                                   acq = "poi")

Best <- Opt.result$Best_Par
Results <- tibble(Parameter = c("num_trees", "k", "alpha", "beta", "nu", "q"),
                  Value = c(round(Best[["num.trees"]]),
                            Grid$k[Snap(Best[["i_k"]], length(Grid$k))],
                            Grid$alpha[Snap(Best[["i_alpha"]], length(Grid$alpha))],
                            Grid$beta[Snap(Best[["i_beta"]], length(Grid$beta))],
                            Grid$nu_q[[Snap(Best[["i_nuq"]], length(Grid$nu_q))]][1],
                            Grid$nu_q[[Snap(Best[["i_nuq"]], length(Grid$nu_q))]][2]))

write.csv(Results, "results/gBART Hyperparameters.csv", row.names = FALSE)
write.csv(Opt.result$History, "results/gBART Hyperparameter History.csv", row.names = FALSE)
