######################################################
# Author: Jonas Wolber

# Date of last change: 24.09.22

# Description: This is the script for the PatternChrome algorithm. PatternChrome can be divided into three major steps:
# 1. Feature extraction: patterns predictive of gene expression are added until training accuracy reaches 99.9%.
# 2. Backward elimination: patterns that are redundant as evaluated on the validation set are removed.
# 3. Hyperparameter tuning: XGBoost hyperparameters are optimized on the validation set.
# 4. Final gene expression prediction: A model based on the training set and using the hyperparameter values is used to predict
# gene expression levels of the test set genes.

# Preparation: Before running the code check how many cores  you have available. You can do so by activating the parallel 
# package and run detectCores(). Set the variable workers to max the amount of available cores. The more workers the faster
# the code execution. Make sure that you have all the input data in the correct directory (ChIP-Seq in Binned_sequencing_data,
# RNA_seq in PatternChrome directory) and that you have created a Datasets directory to store the data if desired. 
# Make sure that all packages are pre-installed.
######################################################

#### Paths ####
setwd("~/PatternChrome")
dataset_path <- "~/PatternChrome/Datasets"
main_path <- "~/PatternChrome/Binned_sequencing_data/"

#### Load libraries ####
library(parallel)
library(xgboost)
library(caTools)
library(pROC)
library(pso)
library(dplyr)

#### Parameters ####
# Metaparameters
workers <- 3 # number of cres used
options(warn=-1)
report <- TRUE # whether to report progress during training
report_frequency <- 5 # report frequency during PSO optimization
rounding_decimal <- 5

# Training parameters
num_cps <- 7
objective_lower <- c(1,1,1.0001,0.25,floor(num_cps/2),rep(0,num_cps))
objective_upper <- c(floor(10000/bin_size),floor(10000/bin_size),5.999,0.75,num_cps,rep(1,num_cps))
max_train_accuracy <- 0.999
num_sample_genes <- 3000
bin_size <- 50
min_distance <- (num_cps*3)^2
max_distance <- max_distance <- (floor(10000/bin_size/3))^2

# PSO parameters during training and hyperparameter tuning
initial_swarm_size <- 20
initial_maxit_stagnate <- 3
initial_maxit <- 20
exploitation_values <- c(0.8, 0.4)
c_p <- 2.05
c_g <- 2.05

# XGBoost parameters during training and backward elimination
xgb_early_stopping_rounds <- 10
nrounds <- 50
eta <- 0.2

# Hyperparameter tuning parameters
hyperparameter_lower <- c(300, 0.005, 0, 1, 1, 1, 0, 0.1)
hyperparameter_upper <- c(700, 0.2, 10, 10, 5, 10, 10, 0.7)
hp_swarm_size <- 20
hp_maxit_stagnate <- 3
hp_maxit <- 20

# Number of k-fold crossvalidation checks
k_fold <- 10

#### Objective function ####
Objective <- function(pars){
  if(!between((pars[1]-pars[2])^2, min_distance, max_distance)){return(-Inf)}
  switch (pars[3],
          hm_train <- H3K4me3_train[,pars[1]:pars[2]],
          hm_train <- H3K4me1_train[,pars[1]:pars[2]],
          hm_train <- H3K36me3_train[,pars[1]:pars[2]],
          hm_train <- H3K27me3_train[,pars[1]:pars[2]],
          hm_train <- H3K9me3_train[,pars[1]:pars[2]]
  ) 
  train_genes <- sample(train_genes, num_sample_genes)
  mp <- pars[6:(5+pars[5])]
  checked_positions <- 1:(ncol(hm_train)-pars[5])
  optim_df <- cbind(train_df[train_genes,],parSapply(cl,train_genes, function(g){
    sum(sapply(checked_positions, function(pos){
      cor(mp, hm_train[g,pos:(pos+length(mp)-1)])>pars[4]
    }), na.rm = T)
  }))
  round(auc(response = optim_df[,1], 
    predictor = predict(xgb.train(
    data = xgb.DMatrix(data = as.matrix(optim_df[,-1]), label = optim_df[,1]), 
    nrounds = nrounds, eta = eta, verbose = 0, objective = "binary:logistic", nthread = workers,
    eval.metric = "auc", tree_method = "hist"), 
    xgb.DMatrix(data = as.matrix(optim_df[,-1]), label = optim_df[,1])), quiet = T),rounding_decimal)
}

#### Hyperparametertuning function ####
hyperparameter_tuning <- function(hyperparams){
  round(auc(response = getinfo(xgb_validation, name = "label"), 
    predictor = predict(xgb.train(data = xgb_train, verbose = 0,  nthread = workers, 
    nrounds = hyperparams[1], eta = hyperparams[2], gamma = hyperparams[3], max_depth = round(hyperparams[4]), 
    lambda = hyperparams[5], alpha = hyperparams[6], min_child_weight = hyperparams[7], subsample = hyperparams[8], 
    objective = "binary:logistic", eval.metric = "auc",tree_method = "hist"), xgb_validation), quiet = T),rounding_decimal)
}

#### Main for loop ####
load("~/PatternChrome/RNA_seq.RData")
setwd(main_path)
cell_lines <-  list.dirs(full.names = F, path = main_path)[-1]
stats <- c()
date <- "24_09"
for (cell_line in cell_lines) {
  print(cell_line)
  t1 <- Sys.time()
  #### Load CHIP_Seq data ####
  setwd(paste(main_path, cell_line, sep = ""))
  load(paste("H3K4me3_",bin_size, "bp_bins.RData", sep =""))
  load(paste("H3K4me1_",bin_size, "bp_bins.RData", sep =""))
  load(paste("H3K36me3_",bin_size, "bp_bins.RData", sep =""))
  load(paste("H3K27me3_",bin_size, "bp_bins.RData", sep =""))
  load(paste("H3K9me3_",bin_size, "bp_bins.RData", sep =""))
  
  #### Load RNA_Seq data and binarize and sort it ####
  RNA <- as.numeric(RNA_seq[rownames(H3K27me3),cell_line])
  RNA <- (RNA > median(RNA)) + 0
  names(RNA) <- rownames(H3K27me3)
  
  #### Train, validation and test subsets ####
  set.seed(42)
  train_genes <- sample(names(RNA),6601)
  validation_test_genes <- names(RNA)[!names(RNA) %in% train_genes]
  H3K4me3_train <- H3K4me3[train_genes,]
  H3K4me1_train <- H3K4me1[train_genes,]
  H3K36me3_train <- H3K36me3[train_genes,]
  H3K27me3_train <- H3K27me3[train_genes,]
  H3K9me3_train <- H3K9me3[train_genes,]
  H3K4me3_validation_test <- H3K4me3[validation_test_genes,]
  H3K4me1_validation_test <- H3K4me1[validation_test_genes,]
  H3K36me3_validation_test <- H3K36me3[validation_test_genes,]
  H3K27me3_validation_test <- H3K27me3[validation_test_genes,]
  H3K9me3_validation_test <- H3K9me3[validation_test_genes,]
  train_df <- data.frame(GE = RNA[train_genes])
  validation_test_df <- data.frame(GE = RNA[validation_test_genes])
  rm(H3K4me3,H3K4me1,H3K36me3,H3K27me3,H3K9me3)
  
  #### Initialisation ####
  improvement_threshold <- 0
  no_improvement_round <- 0
  maxit <-initial_maxit
  swarm_size <- initial_swarm_size
  maxit_stagnate <- initial_maxit_stagnate
  improving <- TRUE
  patterns <- data.frame(Cell_line = 0, Start = 0, End = 0, HM = 0, MP_threshold = 0, Width = 0)
  for (p in seq_len(num_cps)) {
    patterns <- cbind.data.frame(patterns, Pattern = 0)
    colnames(patterns)[length(patterns)] <- paste("Point", p, sep = "_")
  }
  cl <- makePSOCKcluster(workers)
  while(improving){
    
    #### Optimization function ####
    res <- psoptim(par = rep(NA, length(objective_upper)), fn = Objective, 
                   lower = objective_lower, upper = objective_upper, 
                   control = list(
                     fnscale = -1,vectorize = T, abstol = -1, 
                     s = swarm_size, w = exploitation_values, maxit = maxit, c.p = c_p, c.g = c_g, 
                     trace = report, REPORT = report_frequency,
                     maxit.stagnate = maxit_stagnate))
    if(-res$value > improvement_threshold){
      patterns <- rbind(patterns, c(cell_line, res$par))
      switch (res$par[3],
              hm_train <- H3K4me3_train[,res$par[1]:res$par[2]],
              hm_train <- H3K4me1_train[,res$par[1]:res$par[2]],
              hm_train <- H3K36me3_train[,res$par[1]:res$par[2]],
              hm_train <- H3K27me3_train[,res$par[1]:res$par[2]],
              hm_train <- H3K9me3_train[,res$par[1]:res$par[2]]
      ) 
      switch (res$par[3],
              hm_validation_test <- H3K4me3_validation_test[,res$par[1]:res$par[2]],
              hm_validation_test <- H3K4me1_validation_test[,res$par[1]:res$par[2]],
              hm_validation_test <- H3K36me3_validation_test[,res$par[1]:res$par[2]],
              hm_validation_test <- H3K27me3_validation_test[,res$par[1]:res$par[2]],
              hm_validation_test <- H3K9me3_validation_test[,res$par[1]:res$par[2]]
      ) 
      MP_threshold <- res$par[4]
      mp <- res$par[6:(5+res$par[5])]
      checked_positions <- 1:(ncol(hm_train)-length(mp))
      clusterExport(cl, c("hm_train", "hm_validation_test", "MP_threshold", "mp", "checked_positions"))
      train_df <- cbind.data.frame(train_df, parSapply(cl,rownames(train_df),function(g){
        sum(sapply(checked_positions, function(pos){cor(mp, hm_train[g,pos:(pos+length(mp)-1)])>MP_threshold}),na.rm = T)
      }))
      validation_test_df <- cbind.data.frame(validation_test_df, parSapply(cl, rownames(validation_test_df),
        function(g){
          sum(sapply(checked_positions, function(pos){cor(mp, hm_validation_test[g,pos:(pos+length(mp)-1)])>MP_threshold}),na.rm = T)
        }))
      if(report){print(paste("Train accuracy: ", -res$value))} 
      if(-res$value >= max_train_accuracy){
        improving <- FALSE
      }  
    }
    else {
      no_improvement_round <- no_improvement_round + 1
      maxit <- maxit + 1
      swarm_size <- swarm_size + 1
      maxit_stagnate <- initial_maxit + 1
      if(no_improvement_round == maxit_stagnate){
        improving <- FALSE
      }
    }  
  }
  stopCluster(cl)
  print(paste("Training finished with ", -res$value, " accuracy. 
    Number of trained patterns: ", ncol(train_df)-1, ". Time: ", Sys.time()-t1))
  colnames(train_df) <- c("GE",paste("Pattern", 1:(ncol(train_df)-1), sep="_"))
  colnames(validation_test_df) <- colnames(train_df)
  original_train_df <- train_df
  auc_scores <- c()
  
  #### K-fold crossvalidation ####
  for(k in 1:k_fold){
    
    #### Split into validation and test df ####
    set.seed(k)
    split <- sample.split(validation_test_df$GE, SplitRatio = 0.5)
    train_df <- original_train_df
    validation_df <- validation_test_df[split,]
    test_df <- validation_test_df[!split,]
    
    #### Backward elimination ####
    validation_accuracy <- round(auc(response = validation_df[,1], 
      predictor = predict(xgb.train(data = xgb.DMatrix(data = as.matrix(train_df[,-1]), label = train_df[,1]),
      nrounds = nrounds, eta = eta, nthread = workers, objective = "binary:logistic", verbose = 0, 
      eval.metric = "auc", tree_method = "hist"), 
      xgb.DMatrix(data = as.matrix(validation_df[,-1]), label = validation_df[,1])), quiet = T), 
      rounding_decimal)
    pattern <- ncol(train_df)
    while(pattern != 2){
      be_train_df <- train_df[,-pattern]
      be_validation_df <- validation_df[,-pattern]
      be_test_df <- test_df[,-pattern]
      be_accuracy <- round(auc(response = be_validation_df[,1], 
        predictor = predict(xgb.train(
          data = xgb.DMatrix(data = as.matrix(be_train_df[,-1]), label = be_train_df[,1]), 
          nthread = workers, nrounds = nrounds, eta = eta, objective = "binary:logistic", verbose = 0, 
          eval.metric = "auc"), 
          xgb.DMatrix(data = as.matrix(be_validation_df[,-1]), label = be_validation_df[,1])), quiet = T), 
          rounding_decimal)
      if(be_accuracy > validation_accuracy){
        train_df <- be_train_df
        validation_df <- be_validation_df
        test_df <- be_test_df
        validation_accuracy <- be_accuracy
        pattern <- ncol(train_df)
      }
      pattern <- pattern - 1
    }

    #### XGB hyperparameter tuning ####
    xgb_train <- xgb.DMatrix(data = as.matrix(train_df[,-1]), label = train_df[,1])
    xgb_validation <- xgb.DMatrix(data = as.matrix(validation_df[,-1]), label = validation_df[,1])
    res <- psoptim(par = rep(NA, length(hyperparameter_lower)), fn = hyperparameter_tuning, 
      lower = hyperparameter_lower, upper = hyperparameter_upper, 
        control = list(fnscale = -1,vectorize = T, abstol = -1, 
          trace = report, REPORT = report_frequency, 
          s = hp_swarm_size, w = exploitation_values, maxit = hp_maxit, c.p = c_p, c.g = c_g, 
          maxit.stagnate = hp_maxit_stagnate))

    #### Final test accuracy of k-fold ####
    xgb_train <- xgb.DMatrix(data = as.matrix(train_df[,-1]), label = train_df[,1])
    xgb_validation <- xgb.DMatrix(data = as.matrix(validation_df[,-1]), label = validation_df[,1])
    test_accuracy <- round(auc(response = test_df[,1], predictor = predict(xgb.train(data = xgb_train, nthread = workers,
      nrounds = res$par[1], eta = res$par[2], gamma = res$par[3], max_depth = round(res$par[4]), 
      lambda = res$par[5], alpha = res$par[6], min_child_weight = res$par[7], subsample = res$par[8], 
      objective = "binary:logistic", tree_method = "exact",
      watchlist=list(train = xgb_train, test = xgb_validation), 
      early_stopping_rounds = xgb_early_stopping_rounds, verbose = 0, eval.metric = "auc"), 
      xgb.DMatrix(data = as.matrix(test_df[,-1]), label = test_df[,1])), quiet = T), rounding_decimal)
    auc_scores <- c(auc_scores, test_accuracy)
  }
  
  #### Get mean and SD of AUC score for cell line ####
  test_accuracy <- mean(auc_scores)
  test_sd <- sd(auc_scores)
  cell_line_stats <- c(cell_line, test_accuracy, test_sd, Sys.time()-t1, ncol(train_df)-1)
  print(paste("Run completed for cell line ", cell_line, ". Duration: ", Sys.time()-t1, " minutes/hours. Test AUC score: ", test_accuracy, "+-", test_sd, ". Number of patterns learned: ", ncol(train_df)-1))
  
  #### Save statistics and dataframes ####
  setwd(dataset_path)
  stats <- rbind(stats, cell_line_stats)
  dir.create(cell_line)
  setwd(paste(dataset_path,"/", cell_line, sep=""))
  patterns <- patterns[-1,]
  patterns$Start <- floor(as.numeric(patterns$Start))
  patterns$End <- floor(as.numeric(patterns$End))
  patterns$HM <- floor(as.numeric(patterns$HM))
  patterns$Width <- floor(as.numeric(patterns$Width))
  save(patterns, file = paste(cell_line, date, "patterns.RData", sep = "_"))
  train_df <- original_train_df
  save(train_df, file = paste(cell_line, date, "train_df.RData", sep = "_"))
  save(validation_test_df, file = paste(cell_line, date, "validation_test_df.RData", sep = "_"))
}
stopCluster(cl)
setwd(dataset_path)
colnames(stats) <- c("Cell line", "Mean AUC score", "AUC score SD", "Run time", "Number of trained patterns")
rownames(stats) <- 1:length(cell_lines)
write.csv(stats, file = paste(date, "stats.csv"))
