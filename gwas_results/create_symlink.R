# Loading the required libraries

# Getting the names of the files to symlink
file_from <- basename(commandArgs(trailingOnly = TRUE)[1]) # the path to the file to be created
file_to <- commandArgs(trailingOnly = TRUE)[2] # the file that the link should point to

if(file.exists(file_to)) unlink(file_to)
file.symlink(file_from, file_to)

