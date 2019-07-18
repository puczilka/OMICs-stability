setwd('C:/Users/marta/Desktop/IC/term2/Comp Epi/OMICS project')

# Packages
library(dplyr)
library(lme4)

# Load data
covars <- readRDS('Data/Covariates.rds')
proteins <- readRDS('Data/Proteins.rds')


###### Data preparation ####################
## Prepare data for mixed linear models

# Match OMICs data with covariates
index_proteins <- match(rownames(proteins), covars$subjectidp)
covars_pro <- covars[index_proteins,]
all(rownames(covars_pro) == rownames(proteins))   # Quick check that they match

# # Log-transform data and plot distribution
 par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
 set.seed(1)
 for (index in 1:9) {
   plot(density(proteins[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
 }

proteins <- log(proteins)     # Proteins data overwritten with log-transformed data

par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
set.seed(1)
for (index in 1:9) {
  plot(density(proteins[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
 }


###### Obtain vector of ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points).
## Technical confounding from plates and fixed effects need to be accounted for

tol = 10^(-4)    # set a tolerance parameter
icc_proteins <- rep(0, ncol(proteins))     # Zero vector to contain ICCs
names(icc_proteins) <- colnames(proteins)

# Run linear mixed models 
for (column in 1:ncol(proteins)) {
  model1 = lmer(proteins[, column] ~  (1 | id) + (1 | plate) + age + gender + bmi, data = covars_pro)
  vcov = as.data.frame(VarCorr(model1))$vcov
  names(vcov) = as.data.frame(VarCorr(model1))[, 1]
  
  # Check if intra-correlation due to plate is 0
  if (vcov[2] < tol) {
    model1 = lmer(proteins[, column] ~  (1 | id) + age + gender + bmi, data = covars_pro)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    }
  
  icc_proteins[column] = vcov[1]/sum(vcov)
}


###### Create null distribution of ICCs for plates to compare to obtained ICCs ######
## Permute IDs to run LMM on broken structure (to mimic null distribution)

temp <- rep(0,ncol(proteins))     # temp vector to store 'null' ICCs after permutation
null_icc_matrix <- NULL     # Initiate matrix to store 'null' ICCs over all iterations
iterations = 1000


# Run linear mixed models on permuted data
for (i in 1:iterations) {
  permuted <- sample(covars_pro$id)     # Permute IDs 
  covars_pro$permuted <- permuted
  
  for (column in 1:ncol(proteins)) {
    model1 = lmer(proteins[, column] ~ (1 | permuted) + (1 | plate) + age + gender + bmi, data = covars_pro)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    
    # Check if intra-correlation due to plate is 0  
    if ((vcov[2] < tol)) {
      model1 = lmer(proteins[, column] ~ (1 | permuted) + age + gender + bmi, data = covars_pro)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    }
    
    # Create vector of ICCs for each protein
    temp[column] = vcov[1]/sum(vcov)
  }
  
  # Collect matrix of ICCs by protein for each iteration
  null_icc_matrix <- cbind(null_icc_matrix, temp)
}

# Assign row/column names to ICC matrix
rownames(null_icc_matrix) <- colnames(proteins)
colnames(null_icc_matrix) <- paste0('iter ', seq(ncol(null_icc_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that OMICs profiles are UNSTABLE. Test hypothesis by determining how probable
## it is under a suitable threshold (i.e. multiple correction) that observed ICCs belong in the null distribution.
## P-value is obtained from ranking null ICCs and determining where the observed ICC is among the
## null distribution

# Rank null ICC and obtain p-values for observed ICCs
stab_pro <- rep(0, length(icc_proteins))
names(stab_pro) <- names(icc_proteins)

for (column in 1:length(icc_proteins)) {
  stab_pro[column] = sum(null_icc_matrix[column,] >= icc_proteins[column])/length(null_icc_matrix[column,])
  
}


###### Save p-values ######
ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(stab_pro, file = 'Results/stab_pro.rds')





