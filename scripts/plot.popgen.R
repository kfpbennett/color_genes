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

# Input for 10kb windows
pg.c <- read.csv('popgen_chrom.csv') %>%
  drop_na()

# Input for 25kb windows
pg.c <- read.csv('popgen_chrom_25k.csv') %>%
  drop_na()

# Input for 50kb windows
pg.c <- read.csv('popgen_chrom_50k.csv') %>%
  drop_na()


# Plot ------------------------------------------------------------------------

# Plotting using ggplot --------

# Positions and labels for chromosomes in the plots
chrom.positions <- c(
  62000000, 169000000, 300000000, 448000000, 557900000, 620000000, 682000000, 
  746000000, 849200000, 929800000, 1004000000, 1104000000, 1215000000, 
  1308500000, 1525000000)
chrom.labels <- c('1', '1A', '2', '3', '4', '4A', '5', '6', '8','10', '12', 
                  '15', '20', '24', 'Z')

# Base plot to be added on to with each popgen stat
baseplot <- ggplot(data = pg.c %>% filter(sites > 5))+
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


# Color schemes
palette0 <- c('royalblue', 'firebrick1', 'darkorchid1', 'springgreen3')
palette1 <- c('#5CBCFF', '#A696FF', '#E064E6', '#FF009F')
palette2 <- c('#FF2D3B', '#F46A09', '#DD9400', '#BEB600')
palette3 <- c('#FF2D3B', '#DD9400', '#9AD349', '#00FFC8')
palette4 <- c('#CE88B3', '#B0AAEE', '#6FD2FF', '#59EFFD')
palette5 <- c('#A400B1', '#FF4476', '#FFAA60', '#FFD775')
palette6 <- c('#f25a5a', '#cc5af2', '#5aa6f2', '#5af280')
palette7 <- c('#f25a5a', '#b319e6', '#2273c3', '#2e9e4c')
palette8 <- c('#FF0080', '#8000FF', '#0080FF', '#00FF80')
palette9 <- c('#EA4752', '#A505D6', '#647CF2', '#6AB79E')

# Set the color scheme
pg.palette <- palette0

# Make the figures
fstcols <- c(rep(c('gray30', pg.palette[1]), 17), 'gray30')
dxycols <- c(rep(c('gray30', pg.palette[2]), 17), 'gray30')
fdcols <- c(rep(c('gray30', pg.palette[3]), 17), 'gray30')
freqcols <- c(rep(c('gray30', pg.palette[4]), 17), 'gray30')

fstplot <- baseplot+
  scale_color_manual(values = fstcols)+
  geom_point(mapping = aes(x = pos, y = Fst_3_4, color = chrom))+
  geom_segment(x = -7000000, y = -0.007, xend = 1572000000, yend = -0.007)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.4, by = 0.2),
    labels = c('0', '0.2', '0.4'),
    limits = c(-0.01, 0.5),
    expand = c(0, 0)
  )

dxyplot <- baseplot+
  scale_color_manual(values = dxycols)+
  geom_point(mapping = aes(x = pos, y = dxy_3_4, color = chrom))+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.4, by = 0.2),
    labels = c('0', '0.2', '0.4'),
    limits = c(-0.01, 0.56),
    expand = c(0, 0)
  )

fdplot <- baseplot+
  scale_color_manual(values = fdcols)+
  geom_point(mapping = aes(x = pos, y = fd, color = chrom))+
  geom_segment(x = -7000000, y = -0.007, xend = 1572000000, yend = -0.007)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.6, by = 0.2),
    labels = c('0', '0.2', '0.4', '0.6'),
    limits = c(-0.01, 0.64),
    expand = c(0, 0)
  )

freqplot <- baseplot+
  scale_color_manual(values = freqcols)+
  geom_point(mapping = aes(x = pos, y = freq, color = chrom))+
  geom_segment(x = -7000000, y = -0.005, xend = 1572000000, yend = -0.005)+
  scale_y_continuous(
    breaks = seq(from = 0, to = 0.1, by = 0.1),
    labels = c('0', '0.1'),
    limits = c(-0.01, 0.12),
    expand = c(0, 0)
  )


# Plotting using base R --------

fstcols <- c(rep(c('black', 'royalblue'), 17), 'black')
dxycols <- c(rep(c('black', 'firebrick1'), 17), 'black')
fdcols <- c(rep(c('black', 'darkorchid1'), 17), 'black')
freqcols <- c(rep(c('black', 'seagreen2'), 17), 'black')
picols <- c(rep(c('black', 'darkorange2'), 17), 'black')

plot(pg.c[pg.c$sites > 5,]$pos, 
     pg.c[pg.c$sites > 5,]$Fst_3_4, 
     col = fstcols[as.numeric(pg.c[pg.c$sites > 5,]$chrom)], 
     pch = 20)

plot(pg.c[pg.c$sites > 5,]$pos, 
     pg.c[pg.c$sites > 5,]$pi_4, 
     col = picols[as.numeric(pg.c[pg.c$sites > 5,]$chrom)], 
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
