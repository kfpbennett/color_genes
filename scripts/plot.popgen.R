# =============================================================================*
# ---------------------------------- Set up -----------------------------------
# =============================================================================*

# setwd('C:/Users/kbenn/Documents/grad/phd/dissertation/data/pigment/')

# Set directory and load packages ------------------------------------------

library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(data.table)


# Functions ----------------------------------------------------------------

make.bottom <- function(plotPars){
  stats_data <- get(plotPars[[1]])
  scaf <- plotPars[[2]]
  x <- plotPars[[3]]
  y <- plotPars[[4]]
  xlims <- plotPars[[5]]
  xbreaks <- plotPars[[6]]
  xlabs <- plotPars[[7]]
  ylims <- plotPars[[8]]
  ybreaks <- plotPars[[9]]
  ylabs <- plotPars[[10]]
  
  gg.bottom <- 
    ggplot(data = stats_data[stats_data$scaffold == scaf,], 
           mapping = aes(x = get(x), y = get(y)))+
    geom_point()+
    scale_x_continuous(
      breaks = xbreaks, labels = xlabs, limits = xlims, expand = c(0, 0))+
    scale_y_continuous(
      breaks = ybreaks, labels = ylabs, limits = ylims, expand = c(0, 0))+
    theme_bw()+
    theme(legend.position = 'none', 
          panel.background = element_rect(fill = 'gray97'),
          panel.border = element_rect(color = 'gray97'),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.line.y = element_blank())
  return(gg.bottom)
}

prep.gene.data <- function(gplot, anno, ggene = NULL){
  scaf <- as.character(gplot$data$scaffold[1])
  loc <- gplot$scales$scales[[1]]$limits
  genes <- anno %>%
    filter(feature == 'gene', scaffold == scaf, start > loc[1],end < loc[2]) %>%
    pull(gene)
  
  exongff <- anno[anno$gene %in% genes & anno$feature == 'exon',] %>%
    mutate(isGene = gene == ggene)
  
  isFirstTscpt <- vector()
  for(i in 1:length(genes)){
    isFirstTscpt <-
      c(isFirstTscpt,
        exongff[exongff$gene == genes[i],]$transcript_id ==
          exongff[exongff$gene == genes[i],]$transcript_id[1]
      )}
  exongff <- cbind(exongff, isFirstTscpt) %>%
    filter(isFirstTscpt == TRUE)
  
  return(exongff)
}

make.top <- function(dat, plotPars){
  glabel <- plotPars[[11]]
  xlims <- plotPars[[5]]
  ggene <- dat[dat$isGene == TRUE,]$gene %>% unique
  
  glstart <- c()
  glend <- c()
  glgene <- c()
  for(i in 1:length(unique(dat$gene))){
    newdf <- dat[dat$gene == unique(dat$gene)[i],]
    glstart[i] <- min(newdf$start)
    glend[i] <- max(newdf$end)
    glgene[i] <- newdf$gene %>% unique
  }
  glines <- data.frame(start = glstart, end = glend, glgene = glgene) %>%
    mutate(isGene = glgene == ggene)
  
  gg.top <- ggplot(data = dat)+
    geom_segment(data = glines, 
                 aes(x = start, xend = end, y = -4, yend = -4,
                     color = ifelse(isGene == TRUE, 'red', 'gray50')))+
    geom_rect(aes(xmin = start, xmax = end, ymin = -15, ymax = 9,
                  fill = ifelse(isGene == TRUE, 'red', 'gray50')),
              color = ifelse(dat$isGene == TRUE, 'red', 'gray50'),
              size = 0.6)+
    scale_fill_manual(values = c('gray50', 'red'))+
    scale_color_manual(values = c('gray50', 'red'))+
    scale_y_continuous(limits = c(-15, 30))+
    scale_x_continuous(limits = xlims)+
    theme_bw()+
    theme(legend.position = 'none',
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank())
  if(glabel == TRUE){
    gg.top <- gg.top+
      geom_text(mapping = aes(x = dat[dat$isGene == TRUE,]$start[1], y = 22), 
                label = ggene, color = 'red', size = 3.5, nudge_x = 8000)}
  return(gg.top)
}

plot.gene.stats <- function(
  annotation, gene, 
  stats_data, scaf, x, y,
  xlims, xbreaks = NULL, xlabs = NULL, 
  ylims, ybreaks = NULL, ylabs = NULL, 
  glabel = FALSE
){
  plotPars <- list(stats_data, scaf, x, y, xlims, xbreaks, xlabs, 
                   ylims, ybreaks, ylabs, glabel)
  bottom <- make.bottom(plotPars)
  exongff <- prep.gene.data(bottom, annotation, gene)
  top <- make.top(exongff, plotPars)
  
  gtop <- ggplotGrob(top)
  gbot <- ggplotGrob(bottom)
  
  g <- rbind(gtop, gbot, size = 'last')
  panels <- g$layout$t[grepl('panel', g$layout$name)]
  g$heights[panels] <- unit(c(0.85,1.5), 'cm')
  
  grid.newpage()
  grid.draw(g)
}


# =============================================================================*
# --------------------------------- Do stuff ----------------------------------
# =============================================================================*

# Read in the data ------------------------------------------------------------

# Input for 10kb windows
pg.c <- read.csv('popgen_min21_10k.csv') %>%
  drop_na()

# # Input for 25kb windows
# pg.c <- read.csv('popgen_min21_25k.csv') %>%
#   drop_na()

# Read in GFF modified to be R-friendly, keep important info
# Users other than KB will have to change the path
gff <- fread('../genomes/Mvitellinus/Mvit_gff_simplified.txt', sep = '\t') %>%
  separate(col = V9, into = c('other','gene'), sep = 'gene=', remove = FALSE) %>%
  select(scaffold = V1, feature = V3, start = V4, end = V5, gene, V9) %>%
  separate(col = gene, into = c('gene', 'remove1'), sep = ';') %>%
  separate(col = V9,into = c('remove2','transcript'),sep='transcript_id=') %>%
  separate(col = transcript, into = c('transcript_id', 'remove3'), sep = ';') %>%
  select(-c(remove1, remove2, remove3))


# Whole-geome plots ----------------------------------------------------------

# Plotting using ggplot --------

# Positions and labels for chromosomes in the plots
chrom.positions <- c(
  62000000, 169000000, 300000000, 448000000, 557900000, 620000000, 682000000, 
  746000000, 849200000, 929800000, 1004000000, 1104000000, 1215000000, 
  1308500000, 1525000000)
chrom.labels <- c('1', '1A', '2', '3', '4', '4A', '5', '6', '8','10', '12', 
                  '15', '20', '24', 'Z')

# Color schemes
palette0 <- c('royalblue', 'firebrick1', 'darkorchid1', 'springgreen3')
palette1 <- c('#5CBCFF', '#A696FF', '#E064E6', '#FF009F')
palette2 <- c('#FF2D3B', '#FFAA00', '#8CE830', '#00FFBF')
palette3 <- c('#FF2D3B', '#DD9400', '#9AD349', '#00FFC8')
palette4 <- c('#CE88B3', '#B0AAEE', '#6FD2FF', '#59EFFD')
palette5 <- c('#A400B1', '#FF4476', '#FFAA60', '#FFD775')
palette6 <- c('#f25a5a', '#cc5af2', '#5aa6f2', '#5af280')
palette7 <- c('#FF0047', '#E668DE', '#54B6FF', '#00E0FF')
palette8 <- c('#FF0080', '#8000FF', '#0080FF', '#00FF80')
palette9 <- c('#EA4752', '#A505D6', '#647CF2', '#6AB79E')

# Set the color scheme
pg.palette <- palette0

fstcols <- c(rep(c('gray30', pg.palette[1]), 17), 'gray30')
dxycols <- c(rep(c('gray30', pg.palette[2]), 17), 'gray30')
fdcols <- c(rep(c('gray30', pg.palette[3]), 17), 'gray30')
freqcols <- c(rep(c('gray30', pg.palette[4]), 17), 'gray30')


# Make the figures, adjust Y axes as necessary for different window sizes

# Base plot to be added on to with each popgen stat
baseplot <- ggplot(data = pg.c)+
  scale_x_continuous(
    breaks = chrom.positions,
    labels = chrom.labels,
    limits = c(-30000000, 1590000000),
    expand = c(0, 0)
  )+
  theme_bw()+
  theme(legend.position = 'none', 
        panel.background = element_rect(fill = 'gray97'),
        panel.border = element_rect(color = 'gray97'),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),# margin = margin(r = 3)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.y = element_line(color = 'black'))


fstplot <- baseplot+
  scale_color_manual(values = fstcols)+
  geom_point(mapping = aes(x = plotpos, y = Fst_3_4, color = chrom))+
  geom_segment(x = -7000000, y = -0.007, xend = 1572000000, yend = -0.007)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.4, by = 0.2),
    labels = c('0', '0.2', '0.4'),
    limits = c(-0.01, 0.54),
    expand = c(0, 0)
  )

dxyplot <- baseplot+
  scale_color_manual(values = dxycols)+
  geom_point(mapping = aes(x = plotpos, y = dxy_3_4, color = chrom))+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.4, by = 0.1),
    labels = c('0', '0.1', '0.2', '0.3', '0.4'),
    limits = c(-0.01, 0.41),
    expand = c(0, 0)
  )

fdplot <- baseplot+
  scale_color_manual(values = fdcols)+
  geom_point(mapping = aes(x = plotpos, y = fd, color = chrom))+
  geom_segment(x = -7000000, y = -0.007, xend = 1572000000, yend = -0.007)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.6, by = 0.2),
    labels = c('0', '0.2', '0.4', '0.6'),
    limits = c(-0.01, 0.7),
    expand = c(0, 0)
  )

freqplot <- baseplot+
  scale_color_manual(values = freqcols)+
  geom_point(mapping = aes(x = plotpos, y = dSNPs, color = chrom))+
  geom_segment(x = -7000000, y = -0.6, xend = 1572000000, yend = -0.6)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 40, by = 10),
    labels = c('0', '10', '20', '30', '40'),
    limits = c(-1, 45),
    expand = c(0, 0)
  )


# Plotting using base R --------

fstcols <- c(rep(c('black', 'royalblue'), 17), 'black')
dxycols <- c(rep(c('black', 'firebrick1'), 17), 'black')
fdcols <- c(rep(c('black', 'darkorchid1'), 17), 'black')
freqcols <- c(rep(c('black', 'seagreen2'), 17), 'black')
picols <- c(rep(c('black', 'darkorange2'), 17), 'black')

plot(pg.c[pg.c$sites > 45,]$plotpos, 
     pg.c[pg.c$sites > 45,]$Fst_3_4, 
     col = fstcols[as.numeric(pg.c[pg.c$sites > 45,]$chrom)], 
     pch = 20)


# Zoomed plots showing genes above -------------------------------------------

# Fst plot with BCO2 highlighted and labeled
plot.gene.stats(
  annotation = gff,
  gene = 'BCO2',
  stats_data = 'pg.c',
  scaf = 'NW_021940545.1',
  x = 'mid',
  y = 'Fst_3_4',
  xlims = c(5300000, 5900000),
  xbreaks = seq(from = 5400000, to = 5800000, by = 200000),
  xlabs = c('5.4 Mb', '5.6 Mb', '5.8 Mb'),
  ylims = c(-0.01, 0.55),
  ybreaks = c(0, 0.5),
  ylabs = c('0', '0.5'),
  glabel = TRUE
)

# Fst plot with RSPO2 highlighted but not labeled, no ticks on X-axis
plot.gene.stats(
  annotation = gff,
  gene = 'RSPO2',
  stats_data = 'pg.c',
  scaf = 'NW_021939396.1',
  x = 'mid',
  y = 'fd',
  xlims = c(16850000, 17550000),
  ylims = c(-0.1, 1.05),
  ybreaks = c(0, 1),
  ylabs = c('0', '1')
)

plot.gene.stats(
  annotation = gff,
  gene = 'BCO2',
  stats_data = 'pg.c',
  scaf = 'NW_021940545.1',
  x = 'mid',
  y = 'fd',
  xlims = c(5300000, 5900000),
  xbreaks = seq(from = 5400000, to = 5800000, by = 200000),
  xlabs = c('5.4 Mb', '5.6 Mb', '5.8 Mb'),
  ylims = c(-0.1, 0.7),
  ybreaks = c(0, 0.5),
  ylabs = c('0', '0.5'),
  glabel = TRUE
)
