# =============================================================================*
# ---------------------------------- Set up -----------------------------------
# =============================================================================*

# setwd('C:/Users/kbenn/Documents/grad/phd/dissertation/data/pigment/')
library(tidyverse)
library(data.table)

# Files for 10k windows
abbababaFile <- 'abbababa/P1pop3.10kb.pip.csv'
popgenFile <- 'popgen3_10k.csv'

# Files for 25k windows
abbababaFile <- 'abbababa/P1pop3.25kb.pip.csv'
popgenFile <- 'popgenWindows_min21_25k.csv'


# =============================================================================*
# --------------------------------- Functions ---------------------------------
# =============================================================================*

# Sequence: assign everything a chromosome first, then take the reversed ones
# and give them new positions ("mid"). Then use addPos to assign new positions
# for plotting with a break in between each chromosome.

# Takes in a dataframe of windowed analysis and returns the dataframe with 
# an extra column for chromosome
assign.chrom <- function(wins, chrom_list) {
  
  s_ord <- c()
  for(i in 1:length(unique(chrom_list$chrom))) {
    chr <- chrom_list[chrom_list$chrom == unique(chrom_list$chrom)[i],]
    n <- nrow(chr)
    s_ord <- c(s_ord, 1:n)
  }
  chrom_list <- cbind(chrom_list, s_ord)
  
  scafs <- unique(wins$scaffold)
  chroms <- c()
  orders <- c()
  orients <- c()
  c_pos <- c()
  c_or <- c()
  for(i in 1:length(scafs)) {
    chroms[i] <- chrom_list$chrom[which(chrom_list$scaf == scafs[i])]
    orders[i] <- chrom_list$s_ord[which(chrom_list$scaf == scafs[i])]
    orients[i] <- chrom_list$orient[which(chrom_list$scaf == scafs[i])]
    c_pos[i] <- chrom_list$conf_pos[which(chrom_list$scaf == scafs[i])]
    c_or[i] <- chrom_list$conf_orient[which(chrom_list$scaf == scafs[i])]
  }
  rows <- c()
  for(i in 1:length(scafs)) {
    rows[i] <- nrow(filter(wins, scaffold == scafs[i]))
  }
  chrom <- rep(chroms, times = rows)
  order <- rep(orders, times = rows)
  orient <- rep(orients, times = rows)
  conf_pos <- rep(c_pos, times = rows)
  conf_orient <- rep(c_or, times = rows)
  
  output <- cbind(wins, chrom, order, orient, conf_pos, conf_orient) %>%
    mutate(#chrom = as.character(chrom), 
      orient = as.character(orient),
      conf_pos = as.numeric(conf_pos), 
      conf_orient = as.numeric(conf_orient)) %>%
    arrange(chrom, order, start)
  
  return(output)
}

# Takes in popgen data, uses the orient column to flip around reverse-oriented
# scaffolds, and returns the fixed dataframe
fix.orientation <- function(data) {
  df <- data %>% arrange(chrom, order, start) %>%
    group_by(scaffold) %>%
    mutate(rev_start = sort(start, decreasing = TRUE),
           rev_end = sort(end, decreasing = TRUE)) %>%
    mutate(newstart = case_when(orient == '-' ~ rev_start, TRUE ~ start),
           newmid = case_when(orient == '-' ~ mid - start + newstart, TRUE ~ mid),
           newend = case_when(orient == '-' ~ rev_end, TRUE ~ end)) %>%
    ungroup() %>%
    select(scaffold:conf_orient, newstart:newend) %>%
    arrange(chrom, order, newstart) %>%
    as.data.frame()
  return(df)
}

# Takes in dataframe with allele frequencies of each site and filters to
# specifications. Default values are 0.75 like vitellinus for pop4, 0.75
# like candei for pop3, fixed different between pops 2 and 12, and no more
# than 0.25 uncalled for pops 3 and 4, no more than 0.5 for parentals
filter.freqs <- 
  function(data, u=c(0.5, 0.3, 0.3, 0.5), pw=1, w=0.75, y=0.75, py=1) {
    output <- data %>%
    filter(pwfrequ <= u[1], 
           wfrequ <= u[2], 
           yfrequ <= u[3], 
           pyfrequ <= u[4],
           pwfreq1 >= pw, 
           wfreq1 >= w, 
           yfreq0 >= y, 
           pyfreq0 >= py)
    return(output)
}

# Adds to popgen data a column with the proportion of sites passing an allele
# frequency filter out of all sites in each genomic window
add.SNPs <- function(pg_data, freq_data){
  dSNPs <- vector()
  for(i in 1:nrow(pg_data)){
    dSNPs[i] <- nrow(
      freq_data[freq_data$CHROM == pg_data$scaffold[i] & 
                  freq_data$POS >= pg_data[i, 'start'] & 
                  freq_data$POS <= pg_data[i, 'end'], ]
      )}
  return(cbind(pg_data, dSNPs))
}

# Adds a position for easy plotting, with space in between each chromosome
# Slightly crummy function, but gets the job done
add.pos <- function(data){
  data <- arrange(data, chrom, order, newstart)
  # Vector of contig names
  scafs <- unique(data$scaffold)
  
  # Vector of number of rows occupied by each scaffold and last bin of each
  rows <- vector(mode = 'numeric', length = length(scafs))
  maxs <- vector(mode = 'numeric', length = length(scafs))
  for(i in 1:length(rows)){
    dat <- filter(data, scaffold == scafs[i])
    rows[i] <- nrow(dat)
    maxs[i] <- max(dat$newmid)}
  
  # Add up last bins
  addmaxs  <- vector(mode = 'numeric', length = length(maxs))
  for(i in 1:length(maxs)){
    addmaxs[i] <- sum(maxs[i], addmaxs[i-1])}
  
  # Fix the lengths of the vectors for the next step
  addmaxs <- addmaxs[-length(addmaxs)]
  rows <- rows[-1]
  
  # Vector containing the last bin length of each scaffold, repeated the number
  # of the next scaffold's bins times
  master <- rep((addmaxs), times = rows)
  master <- 
    append(master, 
           rep(0, times = length(data[data$scaffold == scafs[1],1])), 
           after = 0)
  
  # Vector of the new positions
  plotpos <- vector(mode = 'numeric', length = nrow(data))
  for(i in 1:length(plotpos)){
    plotpos[i] <- sum(data$newmid[i], master[i])}
  
  # Space out chromosomes
  chrom.maxs <- vector()
  for(i in 2:nrow(data)) {
    if(data[i,]$chrom != data[i-1,]$chrom) {chrom.maxs <- append(chrom.maxs,i)}
  }
  
  for(i in 1:length(chrom.maxs)) {
    plotpos[chrom.maxs[i]:length(plotpos)] <- 
      plotpos[chrom.maxs[i]:length(plotpos)] + 15000000
  }
  # The below creates cpos, the position from the beginning of each chromosome
  lastpos <- vector()
  for(i in 1:length(chrom.maxs)){
    lastpos[i] <- plotpos[chrom.maxs[i] - 1]
  }
  
  crows <- vector()
  for(i in 1:length(unique(data$chrom))) {
    crows[i] <- nrow(data[data$chrom == unique(data$chrom)[i],])
  }
  
  lastpos <- append(lastpos, 0, after = 0)
  lastpos <- rep(lastpos, times = crows)
  
  cpos <- vector()
  for(i in 1:length(plotpos)){
    cpos[i] <- plotpos[i] - lastpos[i] - 15000000
  }
  cpos[1:crows[1]] <- cpos[1:crows[1]] + 15000000
  
  return(cbind(data, cpos, plotpos))
}


# =============================================================================*
# --------------------------------- Do stuff ----------------------------------
# =============================================================================*

# Read in the data ------------------------------------------------------------

# Get list of files from Ragoo output
files <- list.files(path = './ragoo/orderings/', pattern = '*_orderings')

# Read the Ragoo output into a list
orderList <- lapply(paste0('./ragoo/orderings/', files), 
                     read.table, 
                     stringsAsFactors = FALSE)
names(orderList) <- map(str_split(files, '_'), 1)

# Convert the list into a dataframe
chromOrders <- bind_rows(orderList, .id = "column_label") %>%
  rename(chrom = column_label, 
         scaf = V1, 
         orient = V2, 
         conf_pos = V3, 
         conf_orient = V4)

co1 <- chromOrders[as.numeric(as.factor(chromOrders$chrom)) %in% c(1:30, 35),]
co2 <- chromOrders[as.numeric(as.factor(chromOrders$chrom)) %in% c(31:34, 36:82),]
co2$chrom <- 'unplaced'
chromOrders <- rbind(co1, co2)
rm(co1, co2)

# Read in scaffold lenths
lengths <- read.table(
  'GCF_001715985.3_ASM171598v3_assembly_report.txt',
  sep = '\t', header = FALSE, stringsAsFactors = FALSE) %>%
  select(scaffold = V7, length = V9)

# Read in ABBA-BABA results
abbababa <- fread('abbababa/P1pop3.10kb.pip.csv') %>%
  rename(ABmid = mid, ABsites = sites, ABsitesUsed = sitesUsed) %>%
  # fd values below 0 or above 1 are meaningless
  mutate(fd = replace(fd, fd < 0 | fd > 1, 0))

# Read in allele frequency data
freqs <- fread('freqs_min21_filtered.csv') %>%
  rename(CHROM = '#CHROM') %>%
  filter(wfreq1 >= 0.75, yfreq0 >= 0.75)

pg412 <- fread('popgen_4-12_min21_10k.csv') %>%
  mutate(Fst_4_12 = replace(Fst_4_12, Fst_4_12 < 0, 0))

# Read in popgen results, add ABBA-BABA, keep only scaffolds from Ragoo
# results, then add lengths
pg <- fread('popgen3_10k.csv') %>%
  # Negative fst values should be coerced to 0
  mutate(Fst_3_4 = replace(Fst_3_4, Fst_3_4 < 0, 0)) %>%
  left_join(pg412, by = c('scaffold', 'start', 'end', 'mid', 'sites', 'pi_4')) %>%
  left_join(lengths, by = 'scaffold') %>%
  mutate(rnd = dxy_3_4/dxy_4_12, 
         nd = dxy_3_4 - (pi_3 + pi_4) / 2)


# Make the final popgen analysis dataframe ------------------------------------

pgc <- assign.chrom(pg, chromOrders) %>%
  fix.orientation() %>%
  left_join(abbababa, by = c('scaffold', 'start', 'end')) %>%
  add.SNPs(freqs) %>%
  add.pos()

# Output file for 10k windows
write.csv(pgc, 'popgen_min21_rnd_10k.csv', row.names = FALSE)

# Output file for 25k windows
write.csv(pgc, 'popgen_min21_25k.csv', row.names = FALSE)




