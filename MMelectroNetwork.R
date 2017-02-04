# MMelectroNetwork.R
# R version 3.3.1 (2016-06-21)
# February 4, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating pheno network for Brassica data using 
# simone package. Data taken from Brassica control
# and droughted conditions. Focusing on restart for Mt. Moran. 

#-----------------------------------------------------------------------
library(bnlearn)
#-----------------------------------------------------------------------

# Read in middayCO file. 
middayCO <- read.csv(file = "middayCO.csv", row.names = 1)

#### Network analysis with bnlearn. 

# Function for creating blacklist. 
BlackList <- function(x){
  NameList <- colnames(x)
  BlackList <- data.frame(From = NameList[-4], 
                          To = rep("tempNorm", 10))
  BlackList <- rbind(BlackList, cbind(From = NameList[-6], 
                                      To = rep("VpdNorm", 10)))
  BlackList <- rbind(BlackList, cbind(From = NameList[-7], 
                                      To = rep("lightNorm", 10)))
}

# Create blacklist for middayCO. 
bl <- BlackList(middayCO)

# Structure learning with hill-climbing. 
PhenoStr <- hc(middayCO, blacklist = bl, 
               restart = 10, perturb = 1)

# Plot Bayesian Network. 
pdf("middayCObn.pdf")
plot(PhenoStr)  
dev.off()

# Write arcs to a csv file. 
write.csv(PhenoStr$arcs, file = "middayCOarcs.csv")

# Write number of tests to a csv file. 
write.csv(PhenoStr$learning$ntests, file = "ntests.csv")


