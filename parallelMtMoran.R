# parallelMtMoran.R
# R version 3.2.2 (2015-08-14)
# April 27th, 2017. Mallory B. Lai
# Reviewed by: TODO()

# Sources: 
# https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf

#------------------------------------------------------------------------------

# Make a cluster by detecting cores. 
cl <- makeCluster(detectCores())

# Create a vector of repeating Hello World's that is the same length as 
# the number of cores. 
x <- rep("Hello World!", detectCores())

# Cluster version of lapply; applies function print to vector x.
parLapply(cl, x, print)

# Create a vector x of values 1 through the number of cores available.
x <- c(1:detectCores())

# Print values of x from each core. 
parLapply(cl, x, print)

# End cluster. 
stopCluster(cl)
