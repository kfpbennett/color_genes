# =============================================================================*
# ---------------------------------- Set up -----------------------------------
# =============================================================================*

# Change path if you are not Kevin
setwd('C:/Users/kbenn/Documents/grad/phd/dissertation/data/pigment/')
library(tidyverse)
library(data.table)


# =============================================================================*
# --------------------------------- Functions ---------------------------------
# =============================================================================*

# Sequence: assign everything a chromosome first, then take the reversed ones
# and give them new positions ("mid"). Then order them by chromosome following
# the addTruePos method. Then plot chromosome by chromosome.



# Takes in a dataframe of windowed analysis and returns the df with an extra
# column for chromosome
assignChrom <- function(wins, chrom.list) {
  
  s_ord <- c()
  for(i in 1:length(unique(chrom.list$chrom))) {
    chr <- chrom.list[chrom.list$chrom == unique(chrom.list$chrom)[i],]
    n <- nrow(chr)
    s_ord <- c(s_ord, 1:n)
  }
  chrom.list <- cbind(chrom.list, s_ord)
  
  scafs <- unique(wins$scaffold)
  chroms <- c()
  orders <- c()
  orients <- c()
  c_pos <- c()
  c_or <- c()
  for(i in 1:length(scafs)) {
    chroms[i] <- chrom.list$chrom[which(chrom.list$scaf == scafs[i])]
    orders[i] <- chrom.list$s_ord[which(chrom.list$scaf == scafs[i])]
    orients[i] <- chrom.list$orient[which(chrom.list$scaf == scafs[i])]
    c_pos[i] <- chrom.list$conf_pos[which(chrom.list$scaf == scafs[i])]
    c_or[i] <- chrom.list$conf_orient[which(chrom.list$scaf == scafs[i])]
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

# For use in fixOrientation, takes in popgen data filtered to a reverse-
# oriented scaffold, reverses the positions
reverse <- function(data) {
  df <- data
  df[ , c('start', 'end', 'mid')] <- df[ , c('start', 'end', 'mid')] %>%
    arrange(-start)
  return(df)
}

# Takes in popgen data, uses the orient column to flip around reverse-oriented
# scaffolds, and returns the fixed original dataframe
fixOrientation <- function(data) {
  data <- data %>% arrange(chrom, order, start)
  pos.scafs <- data %>% filter(orient == '+')
  rev.scafs <- data %>% filter(orient == '-')
  
  fix.scafs <- data.frame()
  for(i in unique(rev.scafs$scaffold)) {
    fix.scafs <- fix.scafs %>%
      rbind(rev.scafs[rev.scafs$scaffold == i,] %>% reverse())
  }
  
  output <- rbind(pos.scafs, fix.scafs) %>% arrange(chrom, order, start)
  return(output)
}

# Takes in dataframe with allele frequencies of each site and filters to
# specifications. Default values are 0.75 like vitellinus for pop4, 0.75
# like candei for pop3, fixed different between pops 2 and 12, and no more
# than 0.25 uncalled for pops 3 and 4, no more than 0.5 for parentals
filterFreqs <- 
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
addFreq <- function(pg.data, freq.data){
  sites <- vector()
  for(i in 1:nrow(pg.data)){
    loopsites <- nrow(freq.data[freq.data$CHROM == pg.data$scaffold[i] & 
                                 freq.data$POS >= pg.data[i, 'start'] & 
                                 freq.data$POS <= pg.data[i, 'end'], ])
    sites <- append(sites, loopsites, after = length(sites))
  }
  freq <- sites/pg.data$sites
  return(cbind(pg.data, freq))
}

# Adds a position for easy plotting, with space in between each scaf (change)
addPos <- function(data){
  data <- arrange(data, chrom, order, start)
  # Vector of contig names
  scafs <- unique(data$scaffold)
  
  # Vector of number of rows occupied by each scaffold
  rows <- vector(mode = 'numeric', length = length(scafs))
  for(i in 1:length(rows)){
    dat <- filter(data, scaffold == scafs[i])
    rows[i] <- nrow(dat)}
  
  # Vector of last bins of each scaffold
  maxs <- vector(mode = 'numeric', length = length(scafs))
  for(i in 1:length(maxs)){
    dat <- filter(data, scaffold == scafs[i])
    maxs[i] <- max(dat$mid)}
  
  # Add up last bins
  addmaxs  <- vector(mode = 'numeric', length = length(maxs))
  for(i in 1:length(maxs)){
    addmaxs[i] <- sum(maxs[i], addmaxs[i-1])}
  
  # Fix the lengths of the vectors for the next step
  addmaxs <- addmaxs[-length(addmaxs)]
  rows <- rows[-1]
  
  # Vector containing the last bin length of each scaffold, repeated the number
  # of the next contig's bins times
  master <- rep((addmaxs), times = rows)
  master <- 
    append(master, 
           rep(0, times = length(data[data$scaffold == scafs[1],1])), 
           after = 0)
  
  # Vector of the new positions
  pos <- vector(mode = 'numeric', length = nrow(data))
  for(i in 1:length(pos)){
    pos[i] <- sum(data$mid[i], master[i])}
  
  # Space out chromosomes
  chrom.maxs <- vector()
  for(i in 2:nrow(data)) {
    if(data[i,]$chrom != data[i-1,]$chrom) {chrom.maxs <- append(chrom.maxs,i)}
  }
  
  for(i in 1:length(chrom.maxs)) {
    pos[chrom.maxs[i]:length(pos)] <- 
      pos[chrom.maxs[i]:length(pos)] + 15000000
  }
  
  return(cbind(data, pos))
}


# =============================================================================*
# --------------------------------- Do stuff ----------------------------------
# =============================================================================*

# Read in the data ------------------------------------------------------------

# Get list of files from Ragoo output
files <- list.files(path = './ragoo/orderings/', pattern = 'Chr*')

# Read the Ragoo output into a list
order.list <- lapply(paste0('./ragoo/orderings/', files), 
                     read.table, 
                     stringsAsFactors = FALSE)
names(order.list) <- map(str_split(files, '_'), 1)

# Convert the list into a dataframe
chrom.orders <- bind_rows(order.list, .id = "column_label") %>%
  rename(chrom = column_label, 
         scaf = V1, 
         orient = V2, 
         conf_pos = V3, 
         conf_orient = V4)

# Read in scaffold lenths
lengths <- read.table(
  'GCF_001715985.3_ASM171598v3_assembly_report.txt',
  sep = '\t', header = FALSE, stringsAsFactors = FALSE) %>%
  select(scaffold = V7, length = V9)

# Read in ABBA-BABA results
abbababa <- read.csv('abbababa/P1pop3.10kb.csv', stringsAsFactors = FALSE) %>%
  select(scaffold, start, end, ABsites = sites, ABUsed = sitesUsed, fd) %>%
  # fd values below 0 or above 1 are meaningless
  mutate(fd = replace(fd, fd < 0 | fd > 1, 0))

# Read in allele frequency data
freqs <- fread('freqs_fil1_min5.csv') %>%
  filterFreqs()

# Read in popgen results, add ABBA-BABA, keep only scaffolds from Ragoo
# results, then add lengths
pg <- read.csv('popgen_full.csv', header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(Fst_3_4 = replace(Fst_3_4, Fst_3_4 < 0, 0)) %>%
  filter(scaffold %in% chrom.orders$scaf) %>%
  left_join(lengths, by = 'scaffold')


# Make the final popgen analysis dataframe ------------------------------------

pg.c <- assignChrom(pg, chrom.orders) %>%
  fixOrientation() %>%
  left_join(abbababa, by = c('scaffold', 'start', 'end')) %>%
  addFreq(freqs) %>%
  filter(sites > 5) %>%
  addPos()

# Plot ------------------------------------------------------------------------

fstcols <- c(rep(c('black', 'royalblue'), 17), 'black')
dxycols <- c(rep(c('black', 'firebrick1'), 17), 'black')
fdcols <- c(rep(c('black', 'darkorchid1'), 17), 'black')
freqcols <- c(rep(c('black', 'seagreen2'), 17), 'black')

plot(pg.c$pos, 
     pg.c$Fst_3_4, 
     col = fstcols[as.numeric(pg.c$chrom)], 
     pch = 20)

plot(pg.c$pos, 
     pg.c$dxy_3_4, 
     col = dxycols[as.numeric(pg.c$chrom)], 
     pch = 20)

plot(pg.c$pos, 
     pg.c$fd, 
     col = fdcols[as.numeric(pg.c$chrom)], 
     pch = 20)

plot(pg.c$pos, 
     pg.c$freq, 
     ylim = c(0, 0.4),
     col = freqcols[as.numeric(pg.c$chrom)], 
     pch = 20)



# 
# plot(pg.c$pos, 
#      pg.c$Fst_3_4, 
#      pch = 20,
#      col = alpha(plotcols[as.numeric(pg.c$chrom)], pg.c$conf_pos))

