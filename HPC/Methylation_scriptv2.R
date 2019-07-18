setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)

# Load data
covars <- readRDS('Data/Covariates.rds')
methyl <- readRDS('Data/Methylation.rds')     # Very large

# # ------------------
# # SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
# s = sample(ncol(methyl), size = 10)
# methyl <- methyl[,s]
# # ------------------


###### Data preparation ####################
## Prepare data for mixed linear models

# Match OMICs data with covariates
index_methyl <- match(rownames(methyl), covars$subjectidp)
covars_met <- covars[index_methyl,]
all(rownames(covars_met) == rownames(methyl))   # Quick check that they match

# Log-transform data and view distribution
par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
set.seed(1)
for (index in 1:9) {
  plot(density(methyl[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
}

methyl <- log(methyl)     # Methylation data overwritten with log-transformed data

par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
set.seed(1)
for (index in 1:9) {
  plot(density(methyl[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
}


###### Obtain vector of ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points).
## Technical confounding and fixed effects need to be accounted for

tol = 10^(-4)    # set a tolerance parameter
icc_methyl <- rep(0, ncol(methyl))     # Zero vector to contain ICCs
names(icc_methyl) <- colnames(methyl)

# Run linear mixed models 
for (column in 1:ncol(methyl)) {
  model1 = lmer(methyl[, column] ~  (1 | id) + (1 | chip) + (1 | position) + age + gender + bmi, data = covars_met)
  vcov = as.data.frame(VarCorr(model1))$vcov
  names(vcov) = as.data.frame(VarCorr(model1))[, 1]
  ## For chip > tol, position > tol
  
  # Re-run for singular fits (locating and removing source of singular fit)
  if (vcov[2] < tol) {
    
    if(vcov[3] < tol) {
      model1 = lmer(methyl[, column] ~  (1 | id) + age + gender + bmi, data = covars_met)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip < tol, position < tol
    }
    
    else {
      model1 = lmer(methyl[, column] ~  (1 | id) + (1 | position) + age + gender + bmi, data = covars_met)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip < tol, position > tol
    }
  }
  
  else {
    if (vcov[3] < tol) {
      model1 = lmer(methyl[, column] ~  (1 | id) + (1 | chip) + age + gender + bmi, data = covars_met)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For chip > tol, position < tol
    }
  }
    
  icc_methyl[column] = vcov[1]/sum(vcov, na.rm = TRUE)
}


###### Create null distribution of ICCs for plates to compare to obtained ICCs ######
## Permute IDs to run LMM on broken structure (to mimic null distribution)

temp <- rep(0,ncol(methyl))     # temp vector to store 'null' ICCs after permutation
null_icc_matrix <- NULL     # Initiate matrix to store 'null' ICCs over all iterations
iterations = 100

# Run linear mixed models on permuted data
for (i in 1:iterations) {
  permuted <- sample(covars_met$id)     # Permute IDs 
  covars_met$permuted <- permuted
  
  for (column in 1:ncol(methyl)) {
    model1 = lmer(methyl[, column] ~  (1 | permuted) + (1 | chip) + (1 | position) + age + gender + bmi, data = covars_met)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    ## For chip > tol, position > tol
    
    # Re-run for singular fits (locating and removing source of singular fit)
    if (vcov[2] < tol) {
      
      if(vcov[3] < tol) {
        model1 = lmer(methyl[, column] ~  (1 | permuted) + age + gender + bmi, data = covars_met)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip < tol, position < tol
      }
      
      else {
        model1 = lmer(methyl[, column] ~  (1 | permuted) + (1 | position) + age + gender + bmi, data = covars_met)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip < tol, position > tol
      }
    }
    
    else {
      if (vcov[3] < tol) {
        model1 = lmer(methyl[, column] ~  (1 | permuted) + (1 | chip) + age + gender + bmi, data = covars_met)
        vcov = as.data.frame(VarCorr(model1))$vcov
        names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        ## For chip > tol, position < tol
      }
    }
    
    
    # Create vector of ICCs for each biomarker
    temp[column] = vcov[1]/sum(vcov, na.rm = TRUE)
  }
  
  # Collect matrix of ICCs by biomarker for each iteration
  null_icc_matrix <- cbind(null_icc_matrix, temp)
}

# Assign row/column names to ICC matrix
rownames(null_icc_matrix) <- colnames(methyl)
colnames(null_icc_matrix) <- paste0('iter ', seq(ncol(null_icc_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that OMICs profiles are STABLE. Test hypothesis by determining how probable
## it is under a suitable threshold (i.e. multiple correction) that observed ICCs were coincidental
## P-value is obtained from ranking null ICCs and determining where the observed ICC is among the
## null distribution

# Rank null ICC and obtain p-values for observed ICCs
p_val_methyl <- rep(0, length(icc_methyl))

for (column in 1:length(icc_methyl)) {
  p_val_methyl[column] = sum(null_icc_matrix[column,] >= icc_methyl[column])/length(null_icc_matrix[column,])

}

###### Save p-values ######
saveRDS(p_val_methyl, 'pvalues_methylation.rds')

# # Set threshold for p-values
# threshold <- 0.05
# stable_ids <- NULL