setwd('C:/Users/marta/Desktop/IC/term2/Comp Epi/OMICS project')

# Packages
library(dplyr)
library(lme4)

# Load data
covars <- readRDS('Data/Covariates.rds')
metabolites <- readRDS('Data/Metabolites.rds')


###### Data preparation ####################
## High amounts of missing data prevent ICCs from being calculated for particular biomarkers.
## Variables with too much missing data will show zero variance automatically and need to be removed

# Remove variables with too much missing data
missing_list <- NULL
for (i in 1:ncol(metabolites)) {
  if ( (sum(is.na(metabolites[,i]))/length(metabolites[,i])) > 0.95) {
    missing_list <- c(missing_list, i)
  }
}

metabolites <- metabolites[, -missing_list]


# # ------------------
# # SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
s = sample(ncol(metabolites), size = 10)

metabolites <- metabolites[,s]
# # ------------------


# Match OMICs data with covariates
index_met <- match(rownames(metabolites), covars$subjectidp)
covars_met <- covars[index_met,]
all(rownames(covars_met) == rownames(metabolites))   # Quick check that they match

# # Log-transform data and plot distribution
# par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
# set.seed(1)
# for (index in 1:9) {
#   plot(density(metabolites[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
# }

metabolites <- log(metabolites)     # OMIC data overwritten with log-transformed data

# par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
# set.seed(1)
# for (index in 1:9) {
#   plot(density(metabolites[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
# }


###### Obtain vector of ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points).
## Technical confounding from plates and fixed effects need to be accounted for

tol = 10^(-4)    # set a tolerance parameter
icc_met <- rep(0, ncol(metabolites))     # Zero vector to contain ICCs
names(icc_met) <- colnames(metabolites)


# Run linear mixed models 
for (column in 1:ncol(metabolites)) {
  model1 = lmer(metabolites[, column] ~  (1 | id) + age + gender + bmi, data = covars_met)
  vcov = as.data.frame(VarCorr(model1))$vcov
  names(vcov) = as.data.frame(VarCorr(model1))[, 1]
  
  icc_met[column] = vcov[1]/sum(vcov, na.rm = TRUE)
}


###### Create null distribution of ICCs for plates to compare to obtained ICCs ######
## Permute IDs to run LMM on broken structure (to mimic null distribution)

temp <- rep(0,ncol(metabolites))     # temp vector to store 'null' ICCs after permutation
null_icc_matrix <- NULL     # Initiate matrix to store 'null' ICCs over all iterations
iterations = 10

# Run linear mixed models on permuted data
for (i in 1:iterations) {
  permuted <- sample(covars_met$id)     # Permute IDs 
  covars_met$permuted <- permuted
  
  for (column in 1:ncol(metabolites)) {
    model1 = lmer(metabolites[, column] ~ (1 | permuted) + age + gender + bmi, data = covars_met)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    
    # Create vector of ICCs for each protein
    temp[column] = vcov[1]/sum(vcov, na.rm = TRUE)
  }
  
  # Collect matrix of ICCs by protein for each iteration
  null_icc_matrix <- cbind(null_icc_matrix, temp)
}

# Assign row/column names to ICC matrix
rownames(null_icc_matrix) <- colnames(metabolites)
colnames(null_icc_matrix) <- paste0('iter ', seq(ncol(null_icc_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that OMICs profiles are UNSTABLE. Test hypothesis by determining how probable
## it is under a suitable threshold (i.e. multiple correction) that observed ICCs belong in the null distribution.
## P-value is obtained from ranking null ICCs and determining where the observed ICC is among the
## null distribution

# Rank null ICC and obtain p-values for observed ICCs
stab_meta <- rep(0, length(icc_met))
names(stab_meta) <- names(icc_met)

for (column in 1:length(icc_met)) {
  stab_meta[column] = sum(null_icc_matrix[column,] >= icc_met[column])/length(null_icc_matrix[column,])

}


###### Save p-values ######
ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(stab_meta, file = 'Results/stab_meta.rds')
