#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript
#### Blaise Mariner 
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com

library_list <- c(
  "tidyverse",
  "PQLseq2",
  "parallel",
  "doParallel"
)
lapply(library_list, require, character.only = TRUE)

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

# load in the kinship / genetic relatedness matrix
RelatednessMatrix_ <- readRDS(file = "GRM_dog_id.rds")

# get the baseline precision samples
metaData <- readRDS("dap_rrbs-metaData.rds")
metaData <- metaData[metaData$first_rrbs == 'yes' & grepl("precision", metaData$Cohort),]

# mean center the predicted adult height (breed size)
metaData$predicted_height <- scale(metaData$predicted_height, scale = F)

# load in the coverage and methylation matricies. Promoters here as an example.
coverage_all_chr <- readRDS(file =  paste0("coverage_all_chr-PROMOTERS.rds"))
methylation_all_chr <- readRDS(file =  paste0("methylation_all_chr-PROMOTERS.rds"))

# make sure everything is in each other
metaData <- metaData[metaData$lid_pid %in% colnames(coverage_all_chr),]

RelatednessMatrix_<- RelatednessMatrix_[rownames(RelatednessMatrix_) %in% metaData$dog_id, 
                                        colnames(RelatednessMatrix_) %in% metaData$dog_id]
RelatednessMatrix_ <- RelatednessMatrix_[order(rownames(RelatednessMatrix_)), order(colnames(RelatednessMatrix_))]
metaData <- metaData[metaData$dog_id %in% colnames(RelatednessMatrix_),]
metaData <- metaData[order(metaData$dog_id),]

#change the colnames and rownames of the GRM to the lid_pid information to match the colnames of the methylation and coverage mats
rownames(RelatednessMatrix_) <- colnames(RelatednessMatrix_) <- metaData$lid_pid
RelatednessMatrix_<- RelatednessMatrix_[rownames(RelatednessMatrix_) %in% colnames(coverage_all_chr), colnames(RelatednessMatrix_) %in% colnames(coverage_all_chr)]
RelatednessMatrix_ <- RelatednessMatrix_[order(rownames(RelatednessMatrix_)), order(colnames(RelatednessMatrix_))]

coverage_all_chr <- coverage_all_chr[,colnames(coverage_all_chr) %in% metaData$lid_pid]
coverage_all_chr <- coverage_all_chr[,order(colnames(coverage_all_chr))]
methylation_all_chr <- methylation_all_chr[,colnames(methylation_all_chr) %in% metaData$lid_pid]
methylation_all_chr <- methylation_all_chr[,order(colnames(methylation_all_chr))]
metaData_in_relatednessMat <- metaData[metaData$lid_pid %in% colnames(coverage_all_chr),]
metaData_in_relatednessMat <- metaData_in_relatednessMat[order(rownames(metaData_in_relatednessMat)),]

#make sure nothing is messed up
if (! (all(colnames(coverage_all_chr) == colnames(methylation_all_chr)) & all(colnames(methylation_all_chr) == colnames(RelatednessMatrix_)) & all(colnames(methylation_all_chr) == rownames(metaData_in_relatednessMat)))
) { stop("troubleshoot your columns and lid_pids")
}

#make your design mat
# ##https://www.xzlab.org/software/pqlseq/PQLseqManual.pdf
design <- model.matrix(
  ~Age_at_sample + Sex_bool + predicted_height + Breed_Status_bool, data = metaData_in_relatednessMat
)

#loop through your variates
for (col_ in 2:length(colnames(design))){
  
  fit_colname <- colnames(design)[col_]
  new_colname <- paste0(fit_colname)
  
  pheno <- (design[,paste(fit_colname)])
  covariates <- as.matrix(design[,colnames(design)[colnames(design) != fit_colname & colnames(design) != "(Intercept)"]])
  
  #have to do this in batches bcs of how many samples we have
  i = 1
  while (i + 2000 < nrow(coverage_all_chr)) {
    print(i)
    j = i + 2000
    print(Sys.time())
    time_start <- Sys.time()
    fit__ = PQLseq2::pqlseq2(Y=methylation_all_chr[i:j,],
                             x=pheno,
                             K=RelatednessMatrix_,
                             lib_size=coverage_all_chr[i:j,],
                             model="BMM",
                             W = covariates,
                             ncores = 18)
    print(Sys.time() - time_start)
    if (i == 1) {fit_ <- fit__} else {fit_ <- rbind(fit_, fit__)}
    i = j + 1
  }
  j = nrow(methylation_all_chr)
  print(Sys.time())
  time_start <- Sys.time()
  fit__ = PQLseq2::pqlseq2(Y=methylation_all_chr[i:j,],
                           x=pheno,
                           K=RelatednessMatrix_,
                           lib_size=coverage_all_chr[i:j,],
                           model="BMM",
                           W = covariates,
                           ncores = 18)
  print(Sys.time() - time_start)
  fit_ <- rbind(fit_, fit__)
  fit_$padj <- stats::p.adjust(fit_$pvalue, method = "BH", n = length(fit_$pvalue))
  fit_$fdr <- stats::p.adjust(fit_$pvalue, method = "fdr", n = length(fit_$pvalue))
  
  colnames(fit_) <- paste0(colnames(fit_), "_", new_colname)
  if (col_ == 2){fit <- fit_} else {
    
    if (all(rownames(fit) != rownames(fit_))){
      warning("colnames do not match")
      break
    }
    
    fit <- cbind(fit, fit_)
    }
  
  # save it out.
  write_rds(paste0("pqlseq_res.rds"), x = fit)
}
