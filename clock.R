#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

rm(list = ls())  # Clear the workspace

# Determine if running in RStudio; set the number of cores accordingly
if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set working directory
  n_cores = 4 
} 

#basic clock function
build_clock_leave_one_out <- function(dog_id,
                                      metaData,
                                      perc_meth,
                                      alph = 0.5,
                                      nfolds = 10
) {
  # Read in meta info with known chronological ages/sex and technical variables.
  all_info <- metaData
  
  lid_pids <- colnames(perc_meth)
  meta <- all_info[,c("lid_pid", "Age_at_sample","dog_id", "Cohort", "Breed_size", "Sex")]
  
  n_regions <- nrow(perc_meth)
  
  # Make sure the 'Age_at_sample' column to numeric
  meta$Age_at_sample <- as.numeric(meta$Age_at_sample)
  
  # Read in data - This should be the imputed df you are working with
  epi <- perc_meth
  
  #filter epi for dogs you have meta data for and the metadata for the perc meth you have data for
  epi <- epi[, colnames(epi) %in% meta$lid_pid]
  meta <- meta[meta$lid_pid %in% colnames(epi),]
  
  # Reorder the columns in 'epi' to match the order of 'meta$lid'
  meta <- meta[order(meta$lid_pid),]
  epi <- epi[,order(colnames(epi))]
  
  # Remove test subject(s)
  # SAMP indexes from 1 to N samples
  test_lid_pid <- metaData$lid_pid[metaData$dog_id == dog_id]
  train <- epi[, colnames(epi)!=test_lid_pid]
  test <- epi[, colnames(epi)==test_lid_pid]
  
  # Create a vector of training and test ages for elastic net model construction
  trainage <- meta$Age_at_sample[colnames(epi)!=test_lid_pid]
  testage <- meta$Age_at_sample[colnames(epi)==test_lid_pid]
  test_id <- meta$dog_id[colnames(epi)==test_lid_pid]
  test_lid <- meta$lid_pid[colnames(epi)==test_lid_pid]
  
  #### Elastic-net model building ####
  
  # Using N-fold internal CV, train the elastic net model using the training data
  # Note with larger sample sizes, N-fold internal CV becomes intractable
  
  # safety catch
  for (row_ in 1:nrow(train)){
    if (sum(is.na(train[row_,])) > 0){
      stop(paste("Row", row_, "has NAs. Please impute before running this function."))
    }
  }
  
  ## Transpose the train matrix to be samples x features
  model <- glmnet::cv.glmnet(t(train), trainage, nfolds = nfolds, alpha = alph, standardize = FALSE)
  
  # Predict age using the test sample from parameters that minimized MSE during internal CV
  predicted <- predict(model, newx = t(test), s = "lambda.min")
  
  return(data.frame("lid_pid" = test_lid_pid, "Predicted_age" = predicted[1,1]))
}

perc_meth_imputed <- readRDS("methylImp2_perc_meth_imputed.rds")

metaData <- readRDS("dap_rrbs-metaData.rds")

# make sure the dogs are in the precicsion cohort. "first_rrbs" == "yes" if it is the earlier RRBS sample of that dog.
dogs_oi <- metaData$dog_id[grepl("precision", metaData$Cohort) & metaData$first_rrbs == "yes",]

# make sure the perc_meth and the metaData lid_pids are in each other(they should be)
perc_meth <- perc_meth_imputed[,colnames(perc_meth_imputed) %in% metaData$lid_pid]
metaData <- metaData[metaData$lid_pid %in% colnames(perc_meth),]

#### make clock
if (T){
  dir.create("clock_leaveoneout")
  for (dog_id in (metaData$dog_id)){
    print(dog_id)
    file_path = paste0("clock_leaveoneout/", dog_id, ".table")
    if (file.exists(paste(file_path))){next}
    write.table(dog_id, file_path, col.names = F, row.names = F)
    res_ <- build_clock_leave_one_out(dog_id, 0.5, metaData, perc_meth)
    write.table(res_, file_path, col.names = F, row.names = F)
  }
  
  i=1
  for (file in list.files("clock_leaveoneout/")){
    res_ <- read.table(paste0("clock_leaveoneout/", file))
    if (i == 1){
      res <- res_
      i=2
    } else {
      res <- rbind(res, res_)
    }
  }
  colnames(res) <- c("lid_pid", "predicted")
  saveRDS(res, "clock_result.rds")
}
