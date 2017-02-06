# discRnaBN9.R
# R version 3.3.1 (2016-06-21)
# February 5, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating genetic network for Brassica data using 
# bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

#### Preprocessing: 

# Read in phenotype file. 
rna <- read.csv(file = "BrassicaDEgenes.csv", row.names = 1)

# Remove cluster column. 
rna <- rna[, -49]

# Transpose data. 
rna <- t(rna)

# Convert to data.frame. 
rna <- as.data.frame(rna)

# Discretize rna data. 

# Discretize the data. 
discRNA <- discretize(rna, method = "interval", breaks = 5)

# Split rownames into Treatment/Timepoint and Replicate. 
rnaNames <- str_split_fixed(rownames(rna), '_', 2)

# Create a Timepoint column. 
discRNA$Timepoint <- as.numeric(str_split_fixed(rnaNames[, 1], '', 2)[, 2])

# Create a treatment column named INT. 
discRNA$INT <- as.factor(str_split_fixed(rnaNames[, 1], '', 2)[, 1])

# Subset only Timepoint 9. 
TP9 <- discRNA[discRNA$Timepoint == 9, ]

# Remove Timepoint column from TP9. 
TP9 <- TP9[, -231]

# Create list for blacklist.
tiers <- list("INT", names(TP9)[1:230])

# Create blacklist.
bl <- tiers2blacklist(nodes = tiers)

# Search for network.
bn.t9 <- tabu(TP9, blacklist = bl, score = "bde",
              iss = 10, tabu = 50)

# Write csv of network arcs. 
write.csv(bn.t9$arcs, file = "TP9gene.csv")
