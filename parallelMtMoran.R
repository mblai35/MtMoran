# parallelMtMoran.R
# R version 3.2.2 (2015-08-14)
# April 27th, 2017. Mallory B. Lai
# Reviewed by: TODO()

# Sources: 
# https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf

#------------------------------------------------------------------------------

library(parallel)

# Make a cluster by detecting cores. 
cl <- makeCluster(detectCores())

# Create a vector of repeating Hello World's that is the same length as 
# the number of cores. 
x <- rep("Hello World!", detectCores())

# Cluster version of lapply; applies function print to vector x.
y <- data.frame(parLapply(cl, x, print))

# Create data.frame. 
z <- data.frame(Core = rep(1:detectCores()), 
                Message = t(y))

# Write a .csv file. 
write.csv(z, "parCheck1.csv")

# End cluster. 
stopCluster(cl)
