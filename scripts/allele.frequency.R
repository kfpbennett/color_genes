# =============================================================================*
# ---------------------------------- Set up -----------------------------------
# =============================================================================*

setwd('grad/phd/dissertation/data/pigment')
library(tidyverse)
library(data.table)

# =============================================================================*
# --------------------------------- Functions ---------------------------------
# =============================================================================*

# Function that takes in VCF (w/o header) and separates genotypes into columns
getGenotypes <- function(data, range) {
  df <- data
  for(i in range) {
    df <- df %>% 
      separate(i, sep = '[^([:alnum:]|.)]+', extra = 'drop', remove = FALSE,
               into = c(paste0(colnames(data)[i], '.1'),
                        paste0(colnames(data)[i], '.2')))
    df <- df %>%
      select(1:i,(i+3):ncol(df), i+1, i+2)
  }
  df <- df %>% select(-c(ID, QUAL, FILTER, INFO, FORMAT, range))
  return(df)
}

# Function to be used in addFrequencies, calculates allele frequency by
# counting alleles. This version omits uncalled sites from the denominator.
alleleFreqBySite <- function(data, allele) {
  u <- sum(data == '.')
  if(allele == '.') {freq <- u/length(data)}
  else{
    sum <- sum(data == allele)
    freq <- sum/(length(data) - u)}
  return(freq)
}


# Function that takes in processed VCF and adds allele frequencies of
# ref (1) and alt (0) by population. Also calculates uncalled.
addFrequencies <- 
  function(processed_VCF, remove=FALSE,
           range2=c(5:6,47:48), range3=c(7:20,37:42), 
           range4=c(21:36,43:46), range12=49:56) {
  df <- as.data.frame(processed_VCF)
  freqs <- data.frame(
    pwfreq0 = apply(df[,range2], 1, alleleFreqBySite, '0'),
    pwfreq1 = apply(df[,range2], 1, alleleFreqBySite, '1'),
    pwfrequ = apply(df[,range2], 1, alleleFreqBySite, '.'),
    
    wfreq0 = apply(df[,range3], 1, alleleFreqBySite, '0'),
    wfreq1 = apply(df[,range3], 1, alleleFreqBySite, '1'),
    wfrequ = apply(df[,range3], 1, alleleFreqBySite, '.'),
    
    yfreq0 = apply(df[,range4], 1, alleleFreqBySite, '0'),
    yfreq1 = apply(df[,range4], 1, alleleFreqBySite, '1'),
    yfrequ = apply(df[,range4], 1, alleleFreqBySite, '.'),
    
    pyfreq0 = apply(df[,range12], 1, alleleFreqBySite, '0'),
    pyfreq1 = apply(df[,range12], 1, alleleFreqBySite, '1'),
    pyfrequ = apply(df[,range12], 1, alleleFreqBySite, '.')
  )
  out <- cbind(processed_VCF, freqs)
  
  if(remove == TRUE){
    out <- out %>% select(-c(range2, range3, range4, range12))
    }
  return(out)
  }

genos <- genos %>% rename(CHROM = '#CHROM')

# =============================================================================*
# --------------------------------- Do stuff ----------------------------------
# =============================================================================*

# Read in VCF, get genotypes in columns, calculate allele frequencies
genos <- fread('vcf/manakin_fil1_min5_nohead.vcf') %>%
  getGenotypes(range = 10:35) %>%
  addFrequencies(remove = FALSE)

freqs <- fread('vcf/manakin_fil1_min5_nohead.vcf') %>%
  getGenotypes(range = 10:35) %>%
  addFrequencies(remove = TRUE)

# Save as a .csv file
write.csv(freqs, 'freqs_fil1_min5.csv', row.names = FALSE)

