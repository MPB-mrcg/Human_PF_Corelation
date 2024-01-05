# ---------- Phased the vcf

# Load necessary libraries
library(data.table)
library(stringr)
library(tibble)

#==============================================
# Function to Calculate minor allele frequency
#==============================================

# Define a function to calculate MAF for a given variant
calculate_maf <- function(variant) {
  # Exclude missing values "./." from the variant
  h = as.character(variant)[which(as.character(variant) != "./.")]
  
  # Split the alleles and create a matrix
  hh = t(do.call(data.frame, strsplit(as.character(h), "/")))
  
  # Count alleles
  allele_counts <- table(hh)
  
  # Get the less common allele
  minor_allele <- min(as.numeric(names(allele_counts)))
  
  # Calculate MAF
  maf <- allele_counts[as.character(minor_allele)] / sum(allele_counts)
  
  return(maf)
}

#==================================
# Function to filter data
#==================================

FilterLoci = function(Genotype, allelicDepth, samples){
  # Read Genotype, allelicDepth, and samples files using fread
  genotypeData1 = fread(Genotype, header = F) # Read Genotype file
  allelicDepthData1 = fread(allelicDepth, header = F) # Read AD file
  samples1 = fread(samples, header = F)
  
  # Calculate MAF for each variant (column)
  maf_values <- apply(genotypeData1[,-1:-4], 1, calculate_maf)
  genotypeData2 = genotypeData1[which(maf_values >= 0.05),]
  allelicDepthData2 = allelicDepthData1[which(maf_values >= 0.05),]
  
  # Calculate missingness for each variant
  missingness = apply(genotypeData2[,-1:-4], 1, function(x) {
    length(which(x == "./."))/length(x)
  })
  
  # Filter variants with missingness <= 0.3
  genotypeData3 = as.data.frame(genotypeData2[which(missingness <= 0.3),])
  allelicDepthData3 = as.data.frame(allelicDepthData2[which(missingness <= 0.3),])
  
  # Calculate missingness for each sample
  missingnessSamp = apply(genotypeData3[,-1:-4], 2, function(x) {
    length(which(x == "./."))/length(x)
  })
  
  # Filter samples with missingness <= 0.1
  genotypeData = genotypeData3[,c(1:4, which(missingnessSamp <= 0.1)+4)]
  allelicDepthData = allelicDepthData3[,c(1:4, which(missingnessSamp <= 0.1)+4)]
  samples = samples1[which(missingnessSamp <= 0.1),]
  
  # Return a list containing filtered data
  return(list(genotypeData, allelicDepthData, samples))
}

#===============================
#Function to phase data
#===============================

# Function to phase VCF data
phase_Vcf = function(dataList, output) {
  genotypeData = dataList[[1]]
  allelicDepthData = dataList[[2]]
  isolates = dataList[[3]]  # read sample file
  isolates = isolates$V1  # get list of samples

  # Save variant info in variable
  if (ncol(genotypeData) == (length(isolates) + 1)) {
    first4Columns = t(do.call(data.frame, (strsplit(genotypeData$V1, ":"))))
    rownames(first4Columns) = NULL
    genotypeData = as.data.frame(subset(genotypeData, select = -1))
    allelicDepthData = as.data.frame(subset(allelicDepthData, select = -1))
  } else if (ncol(genotypeData) == (length(isolates) + 4)) {
    first4Columns = genotypeData[, 1:4]
    genotypeData = as.data.frame(subset(genotypeData, select = -1:-4))
    allelicDepthData = as.data.frame(subset(allelicDepthData, select = -1:-4))
  } else {
    message("number of samples in genotype differs from the number of samples you provided")
    break
  }

  # Write variant Info
  write.table(as.data.frame(first4Columns), "first4Columns.txt", quote = FALSE, row.names = FALSE, col.names = F)

  # Phase the genotype data
  phasedData = matrix(9, nrow = dim(genotypeData)[1], ncol = dim(genotypeData)[2])
  het = "0/0"
  het = unique(unlist(apply(genotypeData, 1, function(x) {
    h = unique(c(het, as.character(x)))
  })))

  # Loop through genotype data and phase the data using allelic depth
  for (j in 1:nrow(genotypeData)) {
    for (k in 1:ncol(genotypeData)) {
      genoDat = unlist(strsplit(as.character(genotypeData[j, k]), "/"))
      if (genoDat[1] == genoDat[2])
        phasedData[j, k] = genoDat[1]
      else {
        target = as.integer(unlist(strsplit(allelicDepthData[j, k], ',')))
        phasedData[j, k] = genoDat[which.max(target)]
      }
    }
  }

  # Convert character to numeric
  phasedData[] <- sapply(phasedData, as.numeric)

  phasedData = as.data.frame(cbind(first4Columns, phasedData))  # bind variant info to phased data
  names(phasedData) = c("CHROM", "POS", "REF", "ALT", isolates)  # rename columns

  phasedData[phasedData == '.'] = NA  # convert "." to NA
  phasedData = as.tibble(phasedData)

  write.table(phasedData, paste0(output, "PhasedData.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(phasedData)
}

# Function to create phase-like data for human
PhaseLike = function(dataList, output) {
  genotypeData = dataList[[1]]
  allelicDepthData = dataList[[2]]
  isolates = dataList[[3]]  # read sample file
  isolates = isolates$V1
  snpInfo = paste0(genotypeData$V1, "_", genotypeData$V2)
  Fourcols = genotypeData[, 1:4]
  genotypeData = genotypeData[, -1:-4]
  colnames(genotypeData) = isolates
  rownames(genotypeData) = snpInfo
  Fourcols = rbind(Fourcols, Fourcols)
  colnames(Fourcols) = c("CHROM", "POS", "REF", "ALT")

  library(stringr)
  xx = apply(genotypeData, 2, function(x) { as.numeric(str_split_fixed(x, "/", 2)) })
  snpInfo = c(snpInfo, paste0(snpInfo, ":1"))
  rownames(xx) = snpInfo
  xx1 = cbind.data.frame(Fourcols, xx)
  xx1 = xx1[order(rownames(xx1)), ]
  fwrite(xx1, paste0(output, "phasedlike.txt", sep = "\t"))
  return(xx1)
}

# Function to convert genotype to haplotype
Genotype2Haplotype = function(phasedData, output) {
  library(dplyr)
  phasedData = as.tibble(phasedData)

  # Get reference and alternative variables
  ref <- phasedData %>% select(REF) %>% pull()
  alt <- phasedData %>% select(ALT) %>% pull()

  phasedData = phasedData %>% select(-c(1:4))  # remove first 4 columns containing variant info
  phasedData = as.data.frame(phasedData)

  haplotype = matrix(NA, nrow = dim(phasedData)[1], ncol = dim(phasedData)[2])

  for (i in 1:nrow(phasedData)) {
    cat("Computing SNPs ", i, "\n")
    alts = unlist(strsplit(alt[i], ','))
    for (j in 1:ncol(phasedData)) {
      if (is.na(phasedData[i, j])) haplotype[i, j] = "?"
      else if (as.numeric(phasedData[i, j]) == 0) haplotype[i, j] = ref[i]
      else haplotype[i, j] = alts[as.numeric(phasedData[i, j])]
    }
  }

  haplotype = as.tibble(haplotype)
  names(haplotype) = names(phasedData)

  write.table(haplotype, paste0(output, ".haplotype.txt"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(haplotype)
}

#=============================
# CONVERT HAPLOTYPE TO FASTA
#============================
# Function to create FASTA file from haplotype data
fasta_from_haplotype <- function(haplotype, sample_names, outputFasta) {
  vec <- vector(mode = "character", length = length(sample_names))
  for (i in 1:length(sample_names)) {
    vec[i] = paste(">", sample_names[i])
  }

  transposer = t(haplotype)

  dat = file(outputFasta, open = "w")

  for (i in 1:nrow(transposer)) {
    cat(vec[i], file = dat, sep = '\n')
    sample_seqs <- as.character(transposer[i,])
    sample_seqs = paste0(sample_seqs, collapse = '')

    # Split long sequences into chunks of 500 characters
    if (nchar(sample_seqs) > 500) {
      cat(substring(sample_seqs, 1, 500), file = dat, sep = '\n')
      RestOfLine <- nchar(substring(sample_seqs, 501, nchar(sample_seqs)))

      if (RestOfLine > 500) {
        cat(substring(sample_seqs, 501, 1000), file = dat, sep = '\n')
        RestOfLine2 <- nchar(substring(sample_seqs, 1001, nchar(sample_seqs)))

        if (RestOfLine2 > 500) {
          cat(substring(sample_seqs, 1001, 1500), file = dat, sep = '\n')
          RestOfLine3 <- nchar(substring(sample_seqs, 1501, nchar(sample_seqs)))

          if (RestOfLine3 > 500) {
            cat(substring(sample_seqs, 1501, 2000), file = dat, sep = '\n')
            cat(substring(sample_seqs, 2001, nchar(sample_seqs)), file = dat, sep = '\n')
          } else {
            cat(substring(sample_seqs, 1501, nchar(sample_seqs)), file = dat, sep = '\n')
          }
        } else {
          cat(substring(sample_seqs, 1001, nchar(sample_seqs)), file = dat, sep = '\n')
        }
      } else {
        cat(substring(sample_seqs, 501, nchar(sample_seqs)), file = dat, sep = '\n')
      }
    } else {
      cat(sample_seqs, file = dat, sep = '\n')
    }
  }

  close(dat)
}

# Function to create trait file
creatTrait = function(hap, name) {
  # load libraries
  library(stringr)

  # set directory and read input files
  metadata = fread("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/Nora_Amplicon_Seq_Samples_information_30_05_2023.txt", header = T)
  haplotype = hap

  sampleID = names(haplotype)   # Get list of samples
  common <- metadata[metadata$SampleID %in% sampleID]  # Common samples between Haplotype and metadata

  common = common[match(sampleID, common$SampleID), ]

  trait.labels <- common %>% select(`Site`) %>% unique %>% pull() %>% sort()

  # Create a matrix of sampleids in the row and sites in the column
  Matrix <- data.frame(
    INDV = common$SampleID,
    Basse = 0,
    Brikama = 0,
    EFSTH = 0,
    Fajikunda = 0
  )

  for (t in 1:nrow(Matrix)) {
    cat("Computing SNPs ", t, "\n")
    for (c in 2:ncol(Matrix)) {
      if (names(Matrix)[c] == common$Site[t]) Matrix[t, c] <- 1
    }
  }

  trait.file = paste0(name, ".trait")
  cat("Begin Traits;", file = trait.file, sep = "\n")
  cat(paste0("Dimensions NTraits=", ncol(Matrix) - 1, ";"), file = trait.file, append = TRUE, sep = "\n")
  cat("Format labels=yes missing=? separator=Comma;", file = trait.file, append = TRUE, sep = "\n")
  cat("TraitLabels ", names(Matrix)[-1], ";", file = trait.file, append = TRUE)
  cat("\nMatrix ", file = trait.file, append = TRUE, sep = "\n")

  Matrix$New <- paste(Matrix[, 2], Matrix[, 3], Matrix[, 4], Matrix[, 5], sep = ',')

  write.table(Matrix[c(1, 6)], file = trait.file, append = TRUE, sep = "\t", col.names = F, quote = F, row.names = F)
  cat(";", file = trait.file, append = TRUE, sep = "\n")
  cat("END", file = trait.file, append = TRUE, sep = "\n")

  rownames(Matrix) = as.character(Matrix[, 1])
}

# Function to guide user to upload FASTA file in MABL database and set output format to Nexus
fasta2nexus = function() {
  message("Upload your fasta file in MABL database and set the output format to Nexus")
  # Define the URL you want to open
  url <- "http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi"
  # Open the URL in the default web browser
  browseURL(url)
}

#Pf Genes
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf/strictlyGenes2")
fileList=list.files(".", "GT.txt", full.names = T)
for (Genotypes in fileList) {
  h = fread(Genotypes, header = F)
  if(nrow(h) > 2){
    nam=gsub("_secondLociFilter_GT.txt", "", Genotypes)
    allelicDepth=paste0(nam, "_secondLociFilter_AD.txt")
    samplesFile=paste0(nam, "_samples.txt")
    outputFasta <- paste0(nam, "_Phased2.fasta")
    
    dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
    if(nrow(dataList[[1]]) > 2){
      phasedData=phase_Vcf(dataList, nam)
      haplotype=Genotype2Haplotype(phasedData, nam)
      fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
      creatTrait(haplotype, nam)
    }
    #fasta2nexus()
  }
}
  
###### Pf strictly Genes
# Set the working directory
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf/strictlyGenes2")

# Get a list of files with names ending in "GT.txt"
fileList = list.files(".", "GT.txt", full.names = TRUE)

# Iterate through each Genotypes file in the list
for (Genotypes in fileList) {
  # Read the Genotypes file using fread
  h = fread(Genotypes, header = FALSE)
  
  # Check if the number of rows is greater than 2
  if (nrow(h) > 2) {
    # Extract the base name from the file
    nam = gsub("_secondLociFilter_GT.txt", "", Genotypes)
    
    # Create file names for allelic depth and samples files
    allelicDepth = paste0(nam, "_secondLociFilter_AD.txt")
    samplesFile = paste0(nam, "_samples.txt")
    
    # Define the output FASTA file
    outputFasta <- paste0(nam, "_Phased2.fasta")
    
    # Filter loci using the FilterLoci function
    dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
    
    # Check if the number of rows in filtered data is greater than 2
    if (nrow(dataList[[1]]) > 2) {
      # Generate phase-like data using PhaseLike function
      phasedData = PhaseLike(dataList, nam)
      
      # Convert genotype data to haplotype using Genotype2Haplotype function
      haplotype = Genotype2Haplotype(phasedData, nam)
      
      # Create FASTA file from haplotype data
      fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
      
      # Create trait file from haplotype data
      creatTrait(haplotype, nam)
    }
  }
}

#HU Genes
# Set the working directory
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfHu/strictlyGenes2")

# Get a list of files with names ending in "GT.txt"
genoFiles = list.files(".", pattern = "GT.txt")

# Iterate through each Genotypes file in the list
for (Genotypes in genoFiles) {
    # Extract the base name from the file
    nam = gsub("_GT.txt", "", Genotypes)
    
    # Create file names for allelic depth and samples files
    allelicDepth = paste0(nam, "_AD.txt")
    samplesFile = "samples.txt" 
    
    # Define the output FASTA file
    outputFasta <- paste0(nam, "Phased.fasta")
    
    # Filter loci using the FilterLoci function
    dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
    
    # Generate phase-like data using PhaseLike function
    phasedData = PhaseLike(dataList, nam)
    
    # Convert genotype data to haplotype using Genotype2Haplotype function
    haplotype = Genotype2Haplotype(phasedData, nam)
    
    # Create FASTA file from haplotype data
    fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
    
    # Create trait file from haplotype data
    creatTrait(haplotype, nam)
}
