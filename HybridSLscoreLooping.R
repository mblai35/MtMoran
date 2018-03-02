# HybridSLscoreLooping.R
# R version 3.4.3 (2017-11-30)
# March 1, 2018. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# The script is meant to loop through hybrid structure learning 
# algorithms with varying thresholds and scores to look for the most
# accurate structure learning method. Since the DREAM4 data is time-
# series data, we're also looking to see whether a static Bayesian
# network can recover the true network. 

#-----------------------------------------------------------------------
library(DREAM4)
library(SummarizedExperiment)
library(bnlearn)
library(parallel)
#-----------------------------------------------------------------------

# Load DREAM4 data and extract gold standard matrix. 
data(dream4_100_01)
names(assays(dream4_100_01))
expressionData <- assays(dream4_100_01)$simulated
names(metadata(dream4_100_01))
goldStandardMatrix <- metadata(dream4_100_01)$goldStandardAdjacencyMatrix

# Convert goldStandardMatrix into bn object for comparison. 
gs <- empty.graph(colnames(goldStandardMatrix))
amat(gs) <- goldStandardMatrix
plot(gs)

# Convert expression data to a data frame for structure learning. 
bnData <- t(expressionData)
bnData <- as.data.frame(bnData)

# Create vectors to loop through. 
methodVect <- c('quantile', 'interval')
breaksVect <- c(2, 3)
bootReps <- 250
thresh <- c(0.4, 0.5, 0.6)
maximize <- c("hc", "tabu")
restrict <- c("aracne", "si.hiton.pc")
score <- c("loglik", "aic", "bic", "bde", "bds", "bdj", "k2", "mbde", "bdla")

# Create dataframe for all combinations of the desired arguments within structure learning. 
strLearn <- expand.grid(method = methodVect, 
                        breaks = breaksVect,  
                        bootstrap = bootReps,
                        restrict = restrict,
                        max = maximize,
                        score = score)

# Add columns to store structure comparison information. 
threshNameVector <- expand.grid(c("tp", "fp", "fn"), thresh)
threshNameVector <- paste(threshNameVector$Var1, threshNameVector$Var2, sep = '')
threshCols <- matrix(NA, nrow = dim(strLearn)[1], ncol = length(threshNameVector))
colnames(threshCols) <- threshNameVector

# Combine with dataframe to store threshold information.
strLearn <- cbind(strLearn, threshCols)

# Create cluster. 
cl = makeCluster(detectCores(), "SOCK")

for (i in 1:(dim(strLearn)[1]))  {
  
  # Discretize data using current method and breaks. 
  bnDisc <- discretize(bnData[-1, ], 
                       method = as.character(strLearn[i, 'method']), 
                       breaks = strLearn[i, 'breaks'])
  
  # Run structure learning.
  arcs <- boot.strength(bnDisc, cluster = cl, 
                        R = strLearn[i, 'bootstrap'], m = nrow(data),
                        algorithm = "rsmax2", 
                        algorithm.args = list(restrict = as.character(strLearn[i, 'restrict']), maximize = as.character(strLearn[i, 'max']), 
                                              maximize.args = list(score = as.character(strLearn[i, "score"]))),
                        cpdag = TRUE, debug = FALSE)
  
  # Extract averages for varying thresholds. 
  
  for (j in 1:length(thresh)){
    
    # Find averaged network for each threshold.
    avg <- averaged.network(arcs, threshold = thresh[j])
    # Compare learned network to true network. 
    evaluation <- compare(skeleton(gs), skeleton(avg))
    
    # Extract true positives, false positives, and false negatives. 
    strLearn[i, paste("tp", thresh[j], sep = '')] <- evaluation$tp
    strLearn[i, paste("fp", thresh[j], sep = '')] <- evaluation$fp
    strLearn[i, paste("fn", thresh[j], sep = '')] <- evaluation$fn
    
  }
  
  print(i)
  
  if(i %% 10 == 0){write.csv(strLearn, "intermediate.csv")}
  
}

# Stop cluster
stopCluster(cl)

# Write csv for each method. 
write.csv(strLearn, "hybridScore.csv")


