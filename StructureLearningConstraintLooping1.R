# StructureLearningConstraintLooping1.R
# R version 3.4.3 (2017-11-30)
# January 26, 2018. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# The script is meant to loop through various structure learning 
# algorithms with varying thresholds and cutoffs to look for the most
# accurate structure learning method. Since the DREAM4 data is time-
# series data, we're also looking to see whether a static Bayesian
# network can recover the true network. 

#-----------------------------------------------------------------------
library(DREAM4)
library(SummarizedExperiment)
library(bnlearn)
#-----------------------------------------------------------------------

# Load DREAM4 data and extract gold standard matrix. 
data(dream4_010_01)
names(assays(dream4_010_01))
expressionData <- assays(dream4_010_01)$simulated
names(metadata(dream4_010_01))
goldStandardMatrix <- metadata(dream4_010_01)$goldStandardAdjacencyMatrix

# Convert goldStandardMatrix into bn object for comparison. 
gs <- empty.graph(colnames(goldStandardMatrix))
amat(gs) <- goldStandardMatrix
plot(gs)

# Convert expression data to a data frame for structure learning. 
bnData <- t(expressionData)
bnData <- as.data.frame(bnData)

# Create subsets of bnData to loop through. 
# steadyState <- bnData[1, ]
perturbation1 <- bnData[2:22, ]
perturbation2 <- bnData[23:43, ]
perturbation3 <- bnData[44:64, ]
perturbation4 <- bnData[65:85, ]
perturbation5 <- bnData[86:106, ]
knockout <- bnData[107: 126, ]
multifactorial <- bnData[127:136, ]

# Create vectors to loop through. 
bnDataSub <- list(bnData = bnData,  
                  perturbation1 = perturbation1, 
                  perturbation2 = perturbation2, 
                  perturbation3 = perturbation3, 
                  perturbation4 = perturbation4,
                  perturbation5 = perturbation5,
                  knockout = knockout, 
                  multifactorial = multifactorial)
methodVect <- c('quantile', 'interval')
breaksVect <- c(3, 4, 5)
algo <- c("gs", "iamb", "fast.iamb")
bootReps <- c(250, 500, 1000)
thresh <- c(0.4, 0.6, 0.7, 0.8)

# Create dataframe for all combinations of the desired arguments within structure learning. 
strLearn <- expand.grid(subs <- names(bnDataSub),
                        method = methodVect, 
                        breaks = breaksVect, 
                        algorithm = algo, 
                        bootstrap = bootReps, 
                        threshold = thresh)

# Add columns to store structure comparison information. 
strLearn$TruePos <- NA
strLearn$FalsePos <- NA
strLearn$FalseNeg <- NA

  
  for (i in 1:(dim(strLearn)[1]))  {
   
    # Discretize data using current method and breaks. 
    bnDisc <- discretize(bnDataSub[[strLearn[i, 1]]], 
                         method = as.character(strLearn[i, 'method']), 
                         breaks = strLearn[i, 'breaks'])
    
    # Run structure learning.
    arcs <- boot.strength(bnData, cluster = NULL, 
                          R = strLearn[i, 'bootstrap'], m = nrow(data),
                          algorithm = as.character(strLearn[i, 'algorithm']), 
                          cpdag = TRUE, debug = FALSE)
   
     # Extract averages. 
    avg <- averaged.network(arcs, threshold = strLearn[i, 'threshold'])
    # Compare learned network to true network. 
    evaluation <- compare(skeleton(gs), skeleton(avg))
    
    # Extract true positives, false positives, and false negatives. 
    strLearn$TruePos[i] <- evaluation$tp
    strLearn$FalsePos[i] <- evaluation$fp
    strLearn$FalseNeg[i] <- evaluation$fn
    
  }
  
  # Write csv for each method. 
  write.csv(strLearn, "strLearnCompConstraint1.csv")
  #write.csv(strLearn, paste(names(bnDataSub)[s], ".csv", sep = ''))

