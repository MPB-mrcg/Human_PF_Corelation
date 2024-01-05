# ---------- Phased the vcf

library(data.table)
library(stringr)
library(tibble)

#==============================================
#Function to Calculate minor allele frequency
#==============================================

# Define a function to calculate MAF for a given variant
calculate_maf <- function(variant) {
  h = as.character(variant)[which(as.character(variant) != "./.")]
  hh=t(do.call(data.frame, strsplit(as.character(h), "/")))
  # Count alleles
  allele_counts <- table(hh)
  
  # Get the less common allele
  minor_allele <- min(as.numeric(names(allele_counts)))
  
  # Calculate MAF
  maf <- allele_counts[as.character(minor_allele)] / sum(allele_counts)
  
  return(maf)
}

#==================================
#Function to filter data
#==================================

FilterLoci = function(Genotype, allelicDepth, samples){
  genotypeData1 = fread(Genotype, header = F) #read Genotype file
  allelicDepthData1 = fread(allelicDepth, header = F) #read AD file
  samples1 = fread(samples, header = F)
  
  # Calculate MAF for each variant (column)
  maf_values <- apply(genotypeData1[,-1:-4], 1, calculate_maf)
  genotypeData2 = genotypeData1[which(maf_values >= 0.05),]
  allelicDepthData2 = allelicDepthData1[which(maf_values >= 0.05),]
  
  
  missingness = apply(genotypeData2[,-1:-4], 1, function(x) {
    length(which(x == "./."))/length(x)
  })
 
  genotypeData3 = as.data.frame(genotypeData2[which(missingness <= 0.3),])
  allelicDepthData3 = as.data.frame(allelicDepthData2[which(missingness <= 0.3),])
  
  missingnessSamp = apply(genotypeData3[,-1:-4], 2, function(x) {
    length(which(x == "./."))/length(x)
  })
  genotypeData = genotypeData3[,c(1:4, which(missingnessSamp <= 0.1)+4)]
  allelicDepthData = allelicDepthData3[,c(1:4, which(missingnessSamp <= 0.1)+4)]
  samples = samples1[which(missingnessSamp <= 0.1),]
  return(list(genotypeData, allelicDepthData, samples))
}


#===============================
#Function to phase data
#===============================

phase_Vcf=function(dataList, output){#(Genotype, allelicDepth, samples, output){
  genotypeData = dataList[[1]]
  allelicDepthData = dataList[[2]]
  isolates=dataList[[3]] #read sample file
  isolates = isolates$V1  #get list of samples
  
  #save variant info in variable
  if (ncol(genotypeData) == (length(isolates)+1)){
    first4Columns=t(do.call(data.frame, (strsplit(genotypeData$V1, ":")))) 
    rownames(first4Columns)=NULL
    #remove the column that contains the variant information from the datasets
    genotypeData = as.data.frame(subset(genotypeData, select=-1))
    allelicDepthData = as.data.frame(subset(allelicDepthData, select=-1))
  }
  else if(ncol(genotypeData) == (length(isolates)+4)){
    first4Columns=genotypeData[,1:4]
    #remove the column that contains the variant information from the datasets
    genotypeData = as.data.frame(subset(genotypeData, select=-1:-4))
    allelicDepthData = as.data.frame(subset(allelicDepthData, select=-1:-4))
  }
  else{
    message("number of samples in genotype differs from the number of samples you provided")
    break
  }
  
  #write variant Info
  write.table(as.data.frame(first4Columns), "first4Columns.txt", quote = FALSE, row.names = FALSE, col.names = F)  
  
  #FileName = gsub(".txt", "PhasedData.txt", Genotypes) #output file name
  
  #Phase the genotype data
  phasedData = matrix(9, nrow=dim(genotypeData)[1], ncol=dim(genotypeData)[2])
  het <- "0/0"
  het=unique(unlist(apply(genotypeData, 1, function(x){
    h = unique(c(het, as.character(x)))
  })))
  
  #loop througth genotype data and phase the data using allelic depth
  for(j in 1:nrow(genotypeData))
  {
    for(k in 1:ncol(genotypeData))
    {
      genoDat = unlist(strsplit(as.character(genotypeData[j,k]), "/"))
      if(genoDat[1] == genoDat[2])
        phasedData[j,k]= genoDat[1]
      else{
        target = as.integer(unlist(strsplit(allelicDepthData[j,k], ',')))
        phasedData[j,k]= genoDat[which.max(target)]
      }
    }
  }
  
  #remove all rows that only call variants at the refrence 
  #first4Columns1 = first4Columns[rowSums(is.na(phasedData)) != 0,]
  
  #phasedData1 = phasedData[rowSums(is.na(phasedData)) != 0,] 
  
  #convert character to numeric
  phasedData[] <- sapply(phasedData, as.numeric)
  
  phasedData = as.data.frame(cbind(first4Columns, phasedData)) #bind variant info to phased data
  names(phasedData) = c("CHROM", "POS", "REF", "ALT", isolates) #rename columns
  
  phasedData[phasedData =='.']=NA #convert "." to NA
 # phasedData <- sapply( phasedData, as.numeric )
  
  
  
  phasedData= as.tibble(phasedData)
  
  write.table(phasedData, paste0(output, "PhasedData.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(phasedData)
}

#phasedData=phase_Vcf(Genotypes, allelicDepth, samples, FileName)

#=============================================
#Function to create phase like data for human
#=============================================
PhaseLike = function(dataList, output)
{
  genotypeData = dataList[[1]]
  allelicDepthData = dataList[[2]]
  isolates=dataList[[3]] #read sample file
  isolates = isolates$V1 
  snpInfo = paste0(genotypeData$V1, "_", genotypeData$V2)
  Fourcols = genotypeData[,1:4]
  genotypeData = genotypeData[,-1:-4]
  colnames(genotypeData) = isolates
  rownames(genotypeData) = snpInfo
  Fourcols = rbind(Fourcols, Fourcols)
  colnames(Fourcols) = c("CHROM", "POS", "REF", "ALT")
  
  library(stringr)
  xx=apply(genotypeData, 2, function(x){as.numeric(str_split_fixed(x, "/", 2))})
  snpInfo = c(snpInfo, paste0(snpInfo, ":1"))
  rownames(xx) = snpInfo
  xx1 = cbind.data.frame(Fourcols, xx)
  xx1 = xx1[order(rownames(xx1)),]
  fwrite(xx1, paste0(output, "phasedlike.txt", sep="\t"))
  return(xx1)
}

#===============================
# CONVERT GENOTYPE TO HAPLOTYPE
#===============================

Genotype2Haplotype = function(phasedData, output){
  library(dplyr)
  phasedData=as.tibble(phasedData)
  
  #get reference and alternative variables 
  ref <- phasedData %>% select(REF) %>% pull()
  alt <- phasedData %>% select(ALT) %>% pull()
  
  phasedData <- phasedData %>% select(-c(1:4)) #remove first 4 columns containing variant info
  phasedData = as.data.frame(phasedData) #change class to dataframe
  #phasedData <- type_convert(phasedData) #target NA
  
  #subset.data <- phasedData[, sample(ncol(phasedData), 150)]
  
  haplotype <- matrix(NA, nrow=dim(phasedData)[1], ncol=dim(phasedData)[2])
  
  
  for (i in 1:nrow(phasedData)) {
    cat("Computing SNPs ", i, "\n")
    alts <- unlist(strsplit(alt[i], ','))
    for (j in 1:ncol(phasedData)) {
      if(is.na(phasedData[i,j])) haplotype[i,j] <-"?"
      else if(as.numeric(phasedData[i,j]) == 0) haplotype[i,j] <- ref[i]
      else haplotype[i,j] = alts[as.numeric(phasedData[i,j])]
    }
  }
  #names(haplotype)=str_extract(names(haplotype), "[^_]+")
  haplotype <- as.tibble(haplotype); names(haplotype) <- names(phasedData)
  write.table(haplotype, paste0(output,".haplotype.txt"), sep = '\t',
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(haplotype)
}
#haplotype=Genotype2Haplotype(phasedData, FileName)
#=============================
# CONVERT HAPLOTYPE TO FASTA
#============================
fasta_from_haplotype <- function(haplotype, sample_names, outputFasta)
{
  vec<- vector(mode="character", length(sample_names))
  for (i in 1: length(sample_names))
  {
    vec[i] = paste(">",sample_names[i])
  }
  
  transposer = t(haplotype)
  
  dat = file(outputFasta,open="w")
  
  for (i in 1:nrow(transposer))
  {
    cat(vec[i],file=dat,sep = '\n')
    sample_seqs <- as.character(transposer[i,])
    sample_seqs = paste0(sample_seqs,collapse ='')
    if(nchar(sample_seqs)>500)
    {
      
      cat(substring(sample_seqs, 1, 500), file=dat, sep = '\n')
      RestOfLine <- nchar(substring(sample_seqs, 501, nchar(sample_seqs)))
      if(RestOfLine >500)
      {
        cat(substring(sample_seqs, 501, 1000), file=dat, sep = '\n')
        RestOfLine2 <- nchar(substring(sample_seqs, 1001, nchar(sample_seqs)))
        if(RestOfLine2 >500)
        {
          cat(substring(sample_seqs, 1001, 1500), file=dat, sep = '\n')
          RestOfLine3 <- nchar(substring(sample_seqs, 1501, nchar(sample_seqs)))
          if(RestOfLine3 >500)
          {
            cat(substring(sample_seqs, 1501, 2000), file=dat, sep = '\n')
            cat(substring(sample_seqs, 2001, nchar(sample_seqs)), file=dat, sep = '\n')
          }
          else
          {
            cat(substring(sample_seqs, 1501, nchar(sample_seqs)), file=dat, sep = '\n')
          }
        }
        else
        {
          cat(substring(sample_seqs, 1001, nchar(sample_seqs)), file=dat, sep = '\n')
        }
        
        
      }
      else
      {
        cat(substring(sample_seqs, 501, nchar(sample_seqs)), file=dat, sep = '\n')
      }
      
    }
    else
      cat(sample_seqs, file = dat, sep = '\n')
  }
  
  close(dat)
}

#-------------------------create trait------------------------------
# CREATE TRAIT FILE
#--------------------------------------------------------------------
creatTrait=function(hap, name){
 
  #load libraries
  library(stringr)
  
  #set directory and read input files
  #setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora")
  metadata=fread("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/Nora_Amplicon_Seq_Samples_information_30_05_2023.txt", header = T)
  haplotype=hap#"Clair3AllMerged_PF_body_rectifiedvariantPos_Anno.snp3.recodefinal_GT.haplotype.txt"
  #haplotype = fread(hap)
  
  sampleID=names(haplotype)   #Get list of samples
  #remove samples
  #removeSamples=fread("popArtRemove.txt", header = F)
  #sampleID=sampleID[-which(sampleID %in% removeSamples$V1)]
  common <- metadata[metadata$SampleID %in% sampleID] #Common samples between Haplotype and metadata
  #controls=cbind.data.frame(ShipmentCode=sampleID[!(sampleID %in% metadata$ShipmentCode)], `Study site`=6)
  #common=rbind(as.data.frame(common), controls)
  common=common[match(sampleID, common$SampleID),]
  trait.labels <- common %>% select(`Site`) %>% 
    unique %>% pull() %>% sort()
  # common$`Site`[common$`Site`==1]="SSZ"
  # common$`Site`[common$`Site`==2]="SGZ"
  # common$`Site`[common$`Site`==3]="ECF"
  # common$`Site`[common$`Site`==4]="HEF"
  # common$`Site`[common$`Site`==5]="CEF"
  # common$`Site`[common$`Site`==6]="NPP"
  trait.labels
  
  #common$`Study site` = LETTERS[as.numeric(common$`Study site`)]
  #create a matrix of sampleids in the row and sites in te column
  Matrix <- data.frame(INDV = common$SampleID,
                       Basse = 0,
                       Brikama = 0,
                       EFSTH = 0,
                       Fajikunda = 0
  )
  for (t in 1:nrow(Matrix)) {
    cat("Computing SNPs ", t, "\n")
    for (c in 2:ncol(Matrix)) {
      if(names(Matrix)[c] == common$Site[t]) Matrix[t, c] <- 1
    }
  }
  
  trait.file = paste0(name, ".trait")#gsub(".haplotype.txt", ".trait", hap)
  
  cat("Begin Traits;", file=trait.file, sep="\n")
  cat(paste0("Dimensions NTraits=",ncol(Matrix)-1,";"), file=trait.file, append=TRUE, sep="\n")
  cat("Format labels=yes missing=? separator=Comma;", file=trait.file, append=TRUE, sep="\n")
  cat("TraitLabels ", names(Matrix)[-1], ";", file=trait.file, append=TRUE)
  cat("\nMatrix ", file=trait.file, append=TRUE, sep="\n")
  
  Matrix$New <- paste(Matrix[,2], Matrix[,3], Matrix[,4], Matrix[,5], sep = ',')
  
  write.table(Matrix[c(1,6)], file=trait.file, append=TRUE, sep="\t",
              col.names = F, quote = F, row.names = F)
  cat(";", file=trait.file, append=TRUE, sep="\n")
  cat("END", file=trait.file, append=TRUE, sep="\n")
  #close(trait.file)
  rownames(Matrix) = as.character(Matrix[,1])
  #write.table(Matrix[,1:5], "trait.txt", row.names = F)
}
#fasta_from_haplotype(haplotype, names(haplotype), outputFasta)

#--------------------------convert to nexus --------------------------------
#Upload your fasta file in MABL database and set the output format to Nexus
#---------------------------------------------------------------------------
fasta2nexus = function(){
  message("Upload your fasta file in MABL database and set the output format to Nexus")
  # Define the URL you want to open
  url <- "http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi"  # Replace with your desired URL
  
  # Open the URL in the default web browser
  browseURL(url)
}

#copy the result and save it in a file with a .nex extention
#-------------------------------------------------------------------------
# fileList=list.files("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf", "GT.txt", full.names = T)
# for (variable in c(2, 3, 4, 7)) {
#   setwd(dirname(fileList[variable]))
#   Genotypes=basename(fileList[variable])#"PF_0003.recode_snp_GT_Filtered10_2.txt"
#   nam=gsub("_secondLociFilter_GT.txt", "", fileList[variable])
#   allelicDepth=paste0(nam, "_secondLociFilter_AD.txt")#"PF_0003.recode_snp_AD_Filtered10_2.txt"
#   samples=paste0(nam, "_samples.txt")#"filteredSample_pf10_2.txt"
#   outputFasta <- paste0(nam, "Phased.fasta")
#   
#   #FileName = gsub(".txt", "", Genotypes)
#   dataList = FilterLoci(Genotypes, allelicDepth)
#   phasedData=phase_Vcf(dataList, samples, nam)
#   haplotype=Genotype2Haplotype(phasedData, nam)
#   fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
# }

fileList=list.files("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf", "GT.txt", full.names = T)
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora")
Genotypes="humanSubset_Clair3AllMerged_HU_body_rectifiedvariantPosSnpSort_firstLociFilter20.recode_GT.txt"

nam=gsub("_firstLociFilter20.recode_GT.txt", "", Genotypes)
allelicDepth=paste0(nam, "_firstLociFilter20.recode_AD.txt") #"_secondLociFilter_AD.txt") #"PF_0003.recode_snp_AD_Filtered10_2.txt"
samples=paste0("sampleList.txt") #"filteredSample_pf10_2.txt"
#FileName = gsub(".txt", "", Genotypes)
outputFasta <- paste0(nam, "Phased.fasta")

phasedData=phase_Vcf(Genotypes, allelicDepth, samples, nam)
haplotype=Genotype2Haplotype(phasedData, nam)
fasta_from_haplotype(haplotype, names(haplotype), outputFasta)


#fileList=list.files("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf", "GT.txt", full.names = T)
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfHu")
Genotypes="humanSubset_Pf3D7_02_v3_samples_firstLociFilter_GT.txt"
genes = "../pf_genes.txt"

nam=gsub("_firstLociFilter_GT.txt", "", Genotypes)
allelicDepth=paste0(nam, "_firstLociFilter_GT.txt")#"_secondLociFilter_AD.txt")#"PF_0003.recode_snp_AD_Filtered10_2.txt"
samplesFile="Pf3D7_02_v3_samples.txt"#"filteredSample_pf10_2.txt"
#FileName = gsub(".txt", "", Genotypes)
outputFasta <- paste0(nam, "Phased.fasta")

dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
phasedData=phase_Vcf(dataList, nam)#(Genotypes, allelicDepth, samples, nam)
haplotype=Genotype2Haplotype(phasedData, nam)
fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
creatTrait(haplotype, nam)
fasta2nexus()


#Human
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfHu")
Genotypes="humanSubset_Pf3D7_02_v3_samples_firstLociFilter_GT.txt"

nam=gsub("_firstLociFilter_GT.txt", "", Genotypes)
allelicDepth=paste0(nam, "_firstLociFilter_GT.txt")#"_secondLociFilter_AD.txt")#"PF_0003.recode_snp_AD_Filtered10_2.txt"
samplesFile="humanSubset_Pf3D7_02_v3_samples.txt"#"filteredSample_pf10_2.txt"
#FileName = gsub(".txt", "", Genotypes)
outputFasta <- paste0(nam, "Phased.fasta")

dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
phasedData=phase_Vcf(dataList, nam)#(Genotypes, allelicDepth, samples, nam)
haplotype=Genotype2Haplotype(phasedData, nam)
fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
creatTrait(haplotype, nam)
fasta2nexus()

#Pf strictly Genes
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf/strictlyGenes2")
fileList=list.files(".", "GT.txt", full.names = T)
for (Genotypes in fileList) {
  h = fread(Genotypes, header = F)
  if(nrow(h) > 2){
    nam=gsub("_secondLociFilter_GT.txt", "", Genotypes)
    allelicDepth=paste0(nam, "_secondLociFilter_AD.txt")#"_secondLociFilter_AD.txt")#"PF_0003.recode_snp_AD_Filtered10_2.txt"
    samplesFile=paste0(nam, "_samples.txt")
    #FileName = gsub(".txt", "", Genotypes)
    outputFasta <- paste0(nam, "_Phased2.fasta")
    
    dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
    if(nrow(dataList[[1]]) > 2){
      phasedData=phase_Vcf(dataList, nam)#(Genotypes, allelicDepth, samples, nam)
      haplotype=Genotype2Haplotype(phasedData, nam)
      fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
      creatTrait(haplotype, nam)
    }
    #fasta2nexus()
  }
}

#HU strictly Genes
#setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfHu/haplotype/chr4Genes")
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfHu/strictlyGenes/Split/pfHuHap")
#pattern="PF3D7_0220800"#PF3D7_0424400, "PF3D7_1035500"#"PF3D7_0424400"#"PF_PF3D7_1035700", PF3D7_1035500
#fileList=list.files(".", pattern, full.names = T)
fileList=list.files(".", pattern = "GT.txt", full.names = T)

#fileList1 = fileList[grep("GT.txt", fileList)][-1]

for (Genotypes in fileList) {
  h = fread(Genotypes, header = T)
  #genes = read.table("../../hu_genes.txt", header = T)
  #chr = unique(h$V1)
  #for (k in 1:nrow(genes)) {
    #huGeneGeno = h[which((h$V1 == genes$Chromosome[k]) & (h$V2 >= genes$Start_Position[k]) & (h$V2 <= genes$End_Position[k])),]
    if(nrow(h) > 2){
      nam=gsub("_GT.txt","", Genotypes)
      
      AD = gsub("_GT.txt","_AD.txt", Genotypes)
     # huGeneAD = AD[which(AD$V1 %in% huGeneGeno$V1 & AD$V2 %in% huGeneGeno$V2),]
      #fwrite(huGeneAD, paste0(nam, "_AD.txt"))  
      #samplesFile=list.files("../../../../genesVcfPf/haplotype", paste0(substr(nam, 15, 25), "_samples.txt"), full.names = T)#gsub("_GT.txt", "_samples.txt", Genotypes)
      samplesFile = gsub("firstLociFilter", "samples.txt", nam)
      outputFasta <- paste0(nam, "_Phased.fasta")
      #fwrite(huGeneGeno, paste0(nam,"_GT.txt"))
      
      dataList = FilterLoci(Genotypes, AD, samplesFile)
      if(nrow(dataList[[1]]) > 2){
        phasedData=PhaseLike(dataList, nam)#(Genotypes, allelicDepth, samples, nam)
        haplotype=Genotype2Haplotype(phasedData, nam)
        fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
        creatTrait(haplotype, nam)
      }
      #fasta2nexus()
    }
  }
  
#Pf strictly Genes
setwd("/Users/marthaanitademba/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Nora/genesVcfPf/strictlyGenes2")
fileList=list.files(".", "GT.txt", full.names = T)
for (Genotypes in fileList) {
  h = fread(Genotypes, header = F)
  if(nrow(h) > 2){
    nam=gsub("_secondLociFilter_GT.txt", "", Genotypes)
    allelicDepth=paste0(nam, "_secondLociFilter_AD.txt")#"_secondLociFilter_AD.txt")#"PF_0003.recode_snp_AD_Filtered10_2.txt"
    samplesFile=paste0(nam, "_samples.txt")
    #FileName = gsub(".txt", "", Genotypes)
    outputFasta <- paste0(nam, "_Phased2.fasta")
    
    dataList = FilterLoci(Genotypes, allelicDepth, samplesFile)
    if(nrow(dataList[[1]]) > 2){
      phasedData=PhaseLike(dataList, nam)#(Genotypes, allelicDepth, samples, nam)
     # haplotype=Genotype2Haplotype(phasedData, nam)
      #fasta_from_haplotype(haplotype, names(haplotype), outputFasta)
      #creatTrait(haplotype, nam)
    }
    #fasta2
  }
}
