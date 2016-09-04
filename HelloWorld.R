# HelloWorld.R
# R version 3.1.2 (2015-08-14)
# September 4, 2016. Mallory B. Lai

# Sources: WyoARCC bootcamps 

#------------------------------------------------------------------------------

# Create character string "Hello World!"
out <- "Hello World!"

# Create matrix with 25 rows of "Hello World!"
out_matrix <- matrix(out, nrow = 25)

# Write to csv.
write.csv(out_matrix, file = "HelloWorldOutput.csv")



