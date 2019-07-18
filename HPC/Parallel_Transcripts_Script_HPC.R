## Parameters

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
m=as.numeric(args[2])

## Loading packages 

library(lme4)
library(parallel)
library(dplyr)

# Load data
covars <- readRDS('Data/Covariates.rds')
transcripts <- readRDS('Data/Transcripts.rds')     # Very large


###### Data preparation ####################
## High amounts of missing data prevent ICCs from being calculated for particular biomarkers.
## Variables with too much missing data will show zero variance automatically and need to be removed

# Remove variables with too much missing data
missing_list <- NULL
for (i in 1:ncol(transcripts)) {
  if ( (sum(is.na(transcripts[,i]))/length(transcripts[,i])) > 0.95) {
    missing_list <- c(missing_list, i)
  }
}

transcripts <- transcripts[transcripts!=missing_list]
# Match OMICs data with covariates
index_trans <- match(rownames(transcripts), covars$subjectidp)
covars_trans <- covars[index_trans,]

tol = 10^(-4)    # set a tolerance parameter
icc_trans <- rep(0, ncol(transcripts))     # Zero vector to contain ICCs
names(icc_trans) <- colnames(transcripts)

# Run linear mixed models 
foo1<- function(X) {
  for (column in 1:ncol(transcripts)) {
    model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | isolation) + (1| labeling) + (1| hybridization) + age + gender + bmi, data = covars_trans)
    vcov = as.data.frame(VarCorr(model1))$vcov
    names(vcov) = as.data.frame(VarCorr(model1))[, 1]
    ## For isolation, labelling, hybridization > tol
    
    # Re-run for singular fits (locating and removing source of singular fit)
    
    if (vcov[2] < tol) {
      
      if(vcov[3] < tol) {
        
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | id) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        ## For isolation, labelling, hybridization < tol
        
        else{
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        ## For isolation, labelling < tol | hybridization > tol
        
      }
      
      else {
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1| labeling) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        
        else {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      }
    }
    
    
    else {
      if(vcov[3] < tol) {
        
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | isolation) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        
        else {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | isolation) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      }
      
      else {
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | isolation) + (1| labeling) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
        
        else {
          model1 = lmer(transcripts[, column] ~  (1 | id) + (1 | isolation) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      }
    }
    
    
    icc_trans[column] = vcov[1]/sum(vcov, na.rm = TRUE)
  }
}

###### Create null distribution of ICCs for plates to compare to obtained ICCs ######
## Permute IDs to run LMM on broken structure (to mimic null distribution)

temp <- rep(0,ncol(transcripts))     # temp vector to store 'null' ICCs after permutation
null_icc_matrix <- NULL     # Initiate matrix to store 'null' ICCs over all iterations
iterations = 250

# Run linear mixed models on permuted data
foo2<- function(X) {
  for (i in 1:iterations) {
    permuted <- sample(covars_trans$id)     # Permute IDs 
    covars_trans$permuted <- permuted
    
    for (column in 1:ncol(transcripts)) {
      
      model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1| labeling) + (1| hybridization) + age + gender + bmi, data = covars_trans)
      vcov = as.data.frame(VarCorr(model1))$vcov
      names(vcov) = as.data.frame(VarCorr(model1))[, 1]
      ## For isolation, labelling, hybridization > tol
      
      # Re-run for singular fits (locating and removing source of singular fit)
      
      if (vcov[2] < tol) {
        
        if(vcov[3] < tol) {
          
          if(vcov[4] < tol) {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
          ## For isolation, labelling, hybridization < tol
          
          else{
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
          ## For isolation, labelling < tol | hybridization > tol
          
        }
        
        else {
          if(vcov[4] < tol) {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1| labeling) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
          
          else {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
        }
      }
      
      
      else {
        if(vcov[3] < tol) {
          
          if(vcov[4] < tol) {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
          
          else {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
        }
        
        else {
          if(vcov[4] < tol) {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1| labeling) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
          
          else {
            model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
            vcov = as.data.frame(VarCorr(model1))$vcov
            names(vcov) = as.data.frame(VarCorr(model1))[, 1]
          }
        }
      }
      
      
      # Create vector of ICCs for each biomarker
      temp[column] = vcov[1]/sum(vcov, na.rm = TRUE)
    }
    
    # Collect matrix of ICCs by biomarker for each iteration
    null_icc_matrix <- cbind(null_icc_matrix, temp)
  }
}

# Assign row/column names to ICC matrix
rownames(null_icc_matrix) <- colnames(transcripts)
colnames(null_icc_matrix) <- paste0('iter ', seq(ncol(null_icc_matrix)))


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that OMICs profiles are UNSTABLE. Test hypothesis by determining how probable
## it is under a suitable threshold (i.e. multiple correction) that observed ICCs belong in the null distribution.
## P-value is obtained from ranking null ICCs and determining where the observed ICC is among the
## null distribution

# Rank null ICC and obtain p-values for observed ICCs
stab_trans <- rep(0, length(icc_trans))
names(stab_trans) <- names(icc_trans)

foo3<- function(X) {
  for (column in 1:length(icc_trans)) {
    stab_trans[column] = sum(null_icc_matrix[column,] >= icc_trans[column])/length(null_icc_matrix[column,])
    
  }
}

ids=as.character(cut(1:ncol(transcripts), breaks = nchunks, labels = 1:nchunks))

t0=Sys.time()
no_cores=detectCores()
cl <- makeCluster(no_cores) 
clusterExport(cl, c("transcripts", "covars", "foo1","foo2","foo3"))
clusterEvalQ(cl, library(lme4))

ICC_matrix=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  X_chunk=transcripts[,ids==k]
  return(apply(X_chunk, 2, FUN = foo1))
})

Null_ICC_matrix=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  X_chunk=transcripts[,ids==k]
  return(apply(X_chunk, 2, FUN = foo2))
})

stab_trans=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  X_chunk=transcripts[,ids==k]
  return(apply(X_chunk, 2, FUN = foo3))
})

stopCluster(cl)
t1=Sys.time()
print(t1-t0)

###### Save p-values ######
pvalues=unlist(pvalues)
ifelse(dir.exists("Results"),"",dir.create("Results"))
saveRDS(pvalues, paste0("Results/univ_pvalues_one_node_several_cores_m", m, ".rds"))
