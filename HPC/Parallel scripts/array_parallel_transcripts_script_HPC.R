args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])

### Your working directory on the server
work_dir<-"/rds/general/user/mmk4218/home/OMICSproject"
setwd(work_dir)

# Packages
library(dplyr)
library(lme4)
library(parallel)


# Load data
covars <- readRDS('Data/Covariates.rds')
transcripts <- readRDS('Data/Transcripts.rds')     # Very large


# # ------------------
# # SUBSETTING DATA. FILE TOO LARGE. COMMENT OUT FOR HPC
s = sample(ncol(transcripts), size = 10)
transcripts <- transcripts[,s]
# # ------------------


# Match OMICs data with covariates
index_trans <- match(rownames(transcripts), covars$subjectidp)
covars_trans <- covars[index_trans,]
#all(rownames(covars_trans) == rownames(transcripts))   # Quick check that they match

# # View distribution
# par(mfrow = c(3,3), mar = c(4, 4, 1, 1))
# set.seed(1)
# for (index in 1:9) {
#   plot(density(transcripts[, index], na.rm = TRUE), lwd = 2, col = "navy", main = "", las = 1, xlab = "", ylab = "")
# }

###### Obtain vector of ICCs for ID ####################
## Use linear mixed models to isolate the variation within IDs (variability within time points).
## Technical confounding and fixed effects need to be accounted for

tol = 10^(-4)    # set a tolerance parameter
icc_trans <- rep(0, ncol(transcripts))     # Zero vector to contain ICCs
#names(icc_trans) <- colnames(transcripts)


foo1=function(column) {
  
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

return(vcov[1]/sum(vcov, na.rm = TRUE))
}


#Quick test to see if the function works : 
# foo1(1)
# test_foo1=sapply(1:ncol(transcripts), FUN = foo1)


foo2=function(column,iterations) {
  set.seed(1)
for (i in 1:iterations) {
  permuted <- sample(covars_trans$id)     # Permute IDs 
  covars_trans$permuted <- permuted
  
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
        
      } else {
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
    } else {
      if(vcov[3] < tol) {
        
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }  else {
          model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      } else {
        if(vcov[4] < tol) {
          model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1| labeling) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        } else {
          model1 = lmer(transcripts[, column] ~  (1 | permuted) + (1 | isolation) + (1| labeling) + (1 | hybridization) + age + gender + bmi, data = covars_trans)
          vcov = as.data.frame(VarCorr(model1))$vcov
          names(vcov) = as.data.frame(VarCorr(model1))[, 1]
        }
      }
    }
    
    
  # Create vector of ICCs for each biomarker
  temp= vcov[1]/sum(vcov, na.rm = TRUE)
  # Collect matrix of ICCs by biomarker for each iteration
  null_icc_matrix <- c(null_icc_matrix, temp)
}
  return(null_icc_matrix)
}

###### Create null distribution of ICCs for plates to compare to obtained ICCs ######
## Permute IDs to run LMM on broken structure (to mimic null distribution)

temp <- rep(0,ncol(transcripts))     # temp vector to store 'null' ICCs after permutation
null_icc_matrix <- NULL     # Initiate matrix to store 'null' ICCs over all iterations
iterations = 10

#QUICK TEST
# foo2(1,iterations)
# test_foo2=sapply(1:ncol(transcripts), FUN = function(x) {foo2(x,iterations)})
# test_foo22=lapply(1:ncol(transcripts), FUN = function(x) {foo2(x,iterations)})


### Parallelise loops
ids=as.character(cut(1:ncol(transcripts), breaks = nchunks, labels = 1:nchunks))

t0=Sys.time()
no_cores=detectCores()
cl <- makeCluster(no_cores) 
clusterExport(cl, c("transcripts", "covars","covars_trans", "foo1","foo2", "ids", "tol", "icc_trans","temp","iterations","null_icc_matrix"))
clusterEvalQ(cl, library(lme4, dplyr))

icc_trans=parSapply(cl=cl, 1:nchunks, FUN=function(k) {
    X_chunk=as.character(transcripts[,ids==k])
    return(apply(X_chunk, 2, FUN= foo1))
})

null_icc_matrix=t(parSapply(cl=cl, 1:nchunks, FUN=function(k) {
    X_chunk=as.character(transcripts[,ids==k])
    return(apply(X_chunk, 2, FUN =foo2))
}))

stopCluster(cl)
t1=Sys.time()
print(t1-t0)


# Assign row/column names to ICC matrix
#rownames(null_icc_matrix) <- colnames(transcripts)
#colnames(null_icc_matrix) <- paste0('iter ', seq(ncol(null_icc_matrix)))
null_icc_matrix=data.matrix(null_icc_matrix)


###### Calculate p-value for observed ICCs ######
## Null hypothesis is that OMICs profiles are UNSTABLE. Test hypothesis by determining how probable
## it is under a suitable threshold (i.e. multiple correction) that observed ICCs belong in the null distribution.
## P-value is obtained from ranking null ICCs and determining where the observed ICC is among the
## null distribution

# Rank null ICC and obtain p-values for observed ICCs
stab_trans <- rep(0, length(icc_trans))
#names(stab_trans) <- names(icc_trans)

for (column in 1:length(icc_trans)) {
  stab_trans[column] = sum(null_icc_matrix[column,] >= icc_trans[column])/length(null_icc_matrix[column,])
  
}


###### Save p-values ######
ifelse(dir.exists("Results"),"",dir.create("Results", showWarnings = FALSE))
saveRDS(stab_trans, file="Results/stab_trans.rds")
