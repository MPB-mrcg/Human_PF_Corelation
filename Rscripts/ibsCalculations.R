#!/usr/bin/env Rscript

# Function to calculate Identity by State (IBS) matrix
calculateIBS = function(dirname, input)
{
  library(data.table)
  
  # Set working directory
  setwd(dirname)
  
  # Read input file into a data.table
  inputfile = fread(input)
  
  # Remove first 3 columns from the input data
  file = as.data.frame(t(inputfile[,-1:-4]))
  
  library(tictoc)
  tic()  # Start timing
  
  # Initialize an empty matrix for the IBS results
  ibsDf = data.frame(matrix(NA, nrow=nrow(file), ncol=nrow(file)), row.names = row.names(file))
  colnames(ibsDf) <- row.names(file)
  
  # Loop through each row of the input file
  for (i in 1:nrow(file))
  {
    message("************************************* Processing sample row [", i , "] *****************************************")
    tic()  # Start timing for each row
    
    # Set the diagonal element of IBS matrix to 1
    ibsDf[i,i] <- 1
    
    # Extract the genotype data for the current sample
    sample1 = as.matrix(file[i,])
    
    # Initialize the counter for the next sample
    k = i+1
    
    # Loop through each subsequent sample
    while((i<k) & (k <= nrow(file)))
    {
      missingD=0
      s1_vs_s2 = numeric(length(sample1))
      
      # Extract genotype data for the next sample
      sample2 = as.matrix(file[k,])
      
      # Loop through each locus
      for (j in 1:length(sample1)) 
      {
        # Compare genotypes at each locus
        scol1 = sample1[j] 
        scol2 = sample2[j] 
        
        # Check for missing data
        if (is.na(scol1) | is.na(scol2))
        {
          missingD = missingD+1
        }
        else if (scol1 != scol2)
          s1_vs_s2[j] = 0
        else 
          s1_vs_s2[j] = 1
      }
      
      # Calculate IBS and update the matrix
      ibs = sum(s1_vs_s2)/(ncol(file)-missingD)
      ibsDf[i,k] <- ibs
      ibsDf[k,i] <- ibs
      
      # Move to the next sample
      k = k+1
    }
    
    toc()  # Stop timing for the current row
  }
  
  # Write the IBS matrix to a new file
  fwrite(ibsDf, gsub(".txt", "IBSMatrix.txt", input), sep="\t")
}

# Call the calculateIBS function with specified directory and input file
calculateIBS("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora", "Clair3AllMerged_PF_body_rectifiedvariantPos_Anno.snp3.recodefinal_IBSInput.txt")
