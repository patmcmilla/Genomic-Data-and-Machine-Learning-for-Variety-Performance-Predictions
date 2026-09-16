## Variety prediction analysis for winter wheat multi-environment trials.
##
## Fits the five models reported in the manuscript, BLUP, gBLUP, BART, gBART and RKHS, over a
## rolling temporal split: train on the 1, 3 or 5 years preceding each testing year and predict
## variety yields in that year, for testing years 2009 to 2018, at both the individual trial
## location and the mega-environment scale.

library(tidyverse)
library(AGHmatrix) ## For calculating G matrix
library(BGLR) ## For gBLUP, BLUP, and RKHS Bayesian estimates
## Increase memory available to BART machine.
options(java.parameters = c("-Xmx32g", "--add-modules=jdk.incubator.vector"))
library(bartMachine) ## For BART

Day_run <- format(Sys.Date(), "%y%m%d")
 
## Set which groupings and training set sizes to run over.
Run.Groupings <- c('location', 'area')
Run.TrainYears <- c(1, 3, 5)
Run.TestYears <- 2009:2018

## BGLR sampler settings. saveAt keeps the sampler traces out of the working directory.
N.iter <- 6000
Burn.in <- 1000
BGLR.saveAt <- file.path(tempdir(), "")

## The gBART configuration selected by gBARTHypOpt.R. Every value the search covers is set
## here, including the three that coincide with the bartMachine defaults, so that rerunning the
## search and pasting its result in cannot leave part of the configuration silently unapplied.
gBART.pars <- list(num_trees = 125, k = 5, alpha = 0.95, beta = 2, nu = 3, q = 0.9)

## BART is not tuned by the search. Only k is moved off its default.
BART.k <- 3

## Bandwidths for the RKHS kernel averaging model
Bandwidths <- c(0.2, 1, 5)

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

set.seed(42)

## BLUP & gBLUP - Main effects for variety, location and year, plus all three pairwise interactions.
Fit.BLUP <- function(Training.data, Testing.data, Grouping){
  Censored.dat <- Testing.data %>%
    mutate(yield = NA) %>%
    bind_rows(Training.data, .)

  Censored.dat$variety <- factor(x=Censored.dat$variety, ordered=TRUE)

  Z.mat <- model.matrix(~factor(Censored.dat[[Grouping]])-1)
  ZVAR <- model.matrix(~Censored.dat$variety-1)
  ZYR <- model.matrix(~factor(Censored.dat$year)-1)

  Marker.dat <- Genetic.dat %>%
    filter(variety %in% Censored.dat$variety) %>%
    column_to_rownames("variety") %>%
    as.matrix()

  Var.order <- sub("^Censored.dat\\$variety", "", colnames(ZVAR))
  Marker.dat <- Marker.dat[Var.order,]

  ZZ <- tcrossprod(Z.mat)
  ZYY <- tcrossprod(ZYR)

  ## inclusion of the year term is dependent on more than one year of data being in the training set.
  if (n_distinct(Training.data$year) > 1) {
    KYE <- ZZ*ZYY
    diag(KYE) <- diag(KYE)+1/200 ;KYE <- KYE/mean(diag(KYE))
    L4 <- t(chol(KYE))
  }

  ## Builds the ETA for whichever relationship matrix is passed in
  Make.ETA <- function(G.Mat){
    L <- t(chol(G.Mat))
    ZL <- ZVAR%*%L

    ZAZ <- tcrossprod(ZL)

    K <- ZZ*ZAZ
    diag(K) <- diag(K)+1/200 ;K<-K/mean(diag(K))
    L3 <- t(chol(K))

    ETA <- list(ENV=list(X=Z.mat,model='BRR'),
                PED=list(X=ZL,model='BRR'),
                AxE=list(X=L3,model='BRR'))
    
    ## inclusion of the year term is dependent on more than one year of data being in the training set.
    if (n_distinct(Training.data$year) > 1) {
      KGY <- ZYY*ZAZ
      diag(KGY) <- diag(KGY)+1/200 ;KGY<-KGY/mean(diag(KGY))
      L5 <- t(chol(KGY))

      ETA <- c(ETA, list(YEAR=list(X=ZYR,model='BRR'),
                         YxE=list(X=L4,model='BRR'),
                         GxY=list(X=L5,model='BRR')))
    }
    ETA
  }

  ## BLUP, where the relationship matrix is the identity and all varieties are assumed unrelated
  fm1.BLUP <- BGLR(y=Censored.dat$yield,
                   ETA=Make.ETA(diag(nrow(Marker.dat))),
                   nIter=N.iter,
                   burnIn=Burn.in,
                   saveAt=BGLR.saveAt,
                   verbose = F)

  ## gBLUP, where the relationship matrix is the VanRaden additive genomic relationship matrix
  G.Mat <- Gmatrix(SNPmatrix = Marker.dat,
                   missingValue = NA,
                   maf = 0.05,
                   method = "VanRaden")

  diag(G.Mat)=diag(G.Mat)+1/1e4

  fm1.gBLUP <- BGLR(y=Censored.dat$yield,
                    ETA=Make.ETA(G.Mat),
                    nIter=N.iter,
                    burnIn=Burn.in,
                    saveAt=BGLR.saveAt,
                    verbose = F)

  tibble(BLUP = tail(fm1.BLUP$yHat, nrow(Testing.data)),
         gBLUP = tail(fm1.gBLUP$yHat, nrow(Testing.data)))
}

## RKHS
Fit.RKHS <- function(Training.data, Testing.data, Grouping){
  Censored.dat <- Testing.data %>%
    mutate(yield = NA) %>%
    bind_rows(Training.data, .)

  Censored.dat$variety <- factor(x=Censored.dat$variety, ordered=TRUE)

  Z.mat <- model.matrix(~factor(Censored.dat[[Grouping]])-1)
  ZVAR <- model.matrix(~Censored.dat$variety-1)
  ZYR <- model.matrix(~factor(Censored.dat$year)-1)

  Marker.dat <- Genetic.dat %>%
    filter(variety %in% Censored.dat$variety) %>%
    column_to_rownames("variety") %>%
    as.matrix()

  Var.order <- sub("^Censored.dat\\$variety", "", colnames(ZVAR))
  Marker.dat <- Marker.dat[Var.order,]

  ZZ <- tcrossprod(Z.mat)
  ZYY <- tcrossprod(ZYR)

  ## Squared Euclidean distance between marker profiles, averaged over the markers and
  ## normalized by the median off diagonal distance so that the bandwidth is scale free.
  Xs <- scale(Marker.dat, center = TRUE, scale = TRUE)
  Xs[is.na(Xs)] <- 0 # markers with zero variance in this split
  D <- as.matrix(dist(Xs, method = "euclidean"))^2 / ncol(Xs)
  D <- D / median(D[upper.tri(D)])

  Norm.K <- function(K){ diag(K) <- diag(K)+1e-6 ; K/mean(diag(K)) }

  ## The kernels are built at the variety level and expanded to the plot level with ZVAR, so
  ## repeated measures of the same variety keep their own rows and no BLUEs are formed.
  ## Term structure mirrors Make.ETA in Fit.BLUP: main effects for variety, environment and
  ## year, plus all three pairwise interactions, with the Gaussian kernel standing in for G.
  ETA <- list(ENV=list(X=Z.mat,model='BRR'))

  ## inclusion of the year term is dependent on more than one year of data being in the training set.
  Year.term <- n_distinct(Training.data$year) > 1
  if (Year.term) {
    ETA <- c(ETA, list(YEAR=list(X=ZYR,model='BRR'),
                       YxE=list(K=Norm.K(ZZ*ZYY),model='RKHS')))
  }

  for(h in Bandwidths){
    KG <- Norm.K(ZVAR%*%exp(-h*D)%*%t(ZVAR))

    ETA[[paste0('G_h', h)]] <- list(K=KG, model='RKHS')
    ETA[[paste0('GE_h', h)]] <- list(K=Norm.K(KG*ZZ), model='RKHS')

    if (Year.term) {
      ETA[[paste0('GY_h', h)]] <- list(K=Norm.K(KG*ZYY), model='RKHS')
    }
  }

  fm1.RKHS <- BGLR(y=Censored.dat$yield,
                   ETA=ETA,
                   nIter=N.iter,
                   burnIn=Burn.in,
                   saveAt=BGLR.saveAt,
                   verbose = F)

  tibble(RKHS = tail(fm1.RKHS$yHat, nrow(Testing.data)))
}

## BART
Fit.BART <- function(Training.data, Testing.data, Grouping){
  Train.Y <- Training.data %>%
    select(yield) %>%
    as.vector()

  Test.Y <- Testing.data %>%
    select(yield) %>%
    as.vector()

  ## BART
  nongenetic.Train.X <- Training.data %>%
    select(variety, all_of(Grouping)) %>%
    as.matrix() %>%
    as.data.frame()

  nongenetic.Test.X <- Testing.data %>%
    select(variety, all_of(Grouping)) %>%
    as.matrix() %>%
    as.data.frame()

  BART.mod <- bartMachine(X = nongenetic.Train.X,
                          y = Train.Y[[1]],
                          k = BART.k,
                          q = 0.9,
                          mem_cache_for_speed = T,
                          run_in_sample = F) # Don't compute in sample statistics to help improve running time

  BART.preds <- bart_predict_for_test_data(BART.mod,
                                           Xtest = nongenetic.Test.X,
                                           ytest = Test.Y[[1]])

  ## gBART
  Train.X <- Training.data %>%
    select(all_of(Grouping), all_of(Marker.cols)) %>%
    as.matrix() %>%
    as.data.frame()

  Test.X <- Testing.data %>%
    select(all_of(Grouping), all_of(Marker.cols)) %>%
    as.matrix() %>%
    as.data.frame()

  gBART.mod <- bartMachine(X = Train.X,
                           y = Train.Y[[1]],
                           num_trees = gBART.pars$num_trees,
                           k = gBART.pars$k,
                           alpha = gBART.pars$alpha,
                           beta = gBART.pars$beta,
                           nu = gBART.pars$nu,
                           q = gBART.pars$q,
                           mem_cache_for_speed = T, # To remove the speed for memory trade-off. It doesn't work for large P small N problems such as this
                           run_in_sample = F)

  gBART.preds <- bart_predict_for_test_data(gBART.mod,
                                            Xtest = Test.X,
                                            ytest = Test.Y[[1]])

  tibble(BART = BART.preds$y_hat,
         gBART = gBART.preds$y_hat)
}

## Rolling analysis ####
Run.Model.fun <- function(Model, Grouping, TrainYears, TestYear){
  message("Fitting ", Model, " for ", TestYear, ", ", Grouping, ", ", TrainYears, " training years")

  Restricted.data <- Merged.dat %>%
    mutate(year = as.numeric(as.character(year))) %>%
    filter(year >= TestYear - TrainYears)

  Training.data <- Restricted.data %>%
    filter(year < TestYear)

  Testing.data <- Restricted.data %>%
    filter(year == TestYear)

  Preds <- get(paste0("Fit.", Model))(Training.data, Testing.data, Grouping)

  Testing.data %>%
    select(yield:area) %>%
    bind_cols(Preds) %>%
    pivot_longer(cols = -c(yield, variety, location, year, area),
                 names_to = "Method",
                 values_to = "Prediction") %>%
    mutate(TrainingYears = TrainYears,
           Grouping = if_else(Grouping == "area", "MegaEnvironment", "Location"))
}

## Set which models to fit.
## Run the analysis
Model.runs <- expand_grid(Model = Run.Models,
                          Grouping = Run.Groupings,
                          TrainYears = Run.TrainYears,
                          TestYear = Run.TestYears)

Prediction.Results <- pmap_dfr(Model.runs, Run.Model.fun)

## Summarise the predictions
Correlation.Results <- Prediction.Results %>%
  group_by(Grouping, TrainingYears, Method, year, location) %>%
  summarise(Correlation = cor(yield, Prediction),
            RankCorr = cor(yield, Prediction, method = "kendall"),
            .groups = "drop") %>%
  group_by(Grouping, TrainingYears, Method) %>%
  summarise(Correlation = mean(Correlation, na.rm = T),
            RankCorr = mean(RankCorr, na.rm = T),
            .groups = "drop")

MSPE.Results <- Prediction.Results %>%
  group_by(Grouping, TrainingYears, Method) %>%
  summarise(MSPE = mean((yield - Prediction)^2),
            .groups = "drop")

Joined.Results <- Correlation.Results %>%
  left_join(., MSPE.Results, by = c("Grouping", "TrainingYears", "Method"))

## Save results
write.csv(Prediction.Results, paste0("results/Raw Predictions_", Day_run, ".csv"), row.names = FALSE)
write.csv(Joined.Results, paste0("results/Accuracy Summary_", Day_run, ".csv"), row.names = FALSE)
