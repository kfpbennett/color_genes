# =============================================================================*
# ---------------------------------- Set up -----------------------------------
# =============================================================================*

# Change path if you are not Kevin
setwd('C:/Users/kbenn/Documents/grad/phd/dissertation/data/pigment/')
library(tidyverse)
library(ggplot2)


# =============================================================================*
# --------------------------------- Do stuff ----------------------------------
# =============================================================================*

# Read in the data ------------------------------------------------------------

pg.c <- read.csv('popgen_chrom.csv', stringsAsFactors = FALSE)


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


