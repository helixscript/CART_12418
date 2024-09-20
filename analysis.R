library(tidyverse)
library(GenomicRanges)
source('lib.R')

# Download intSite data, parse data, cleanup files.
system('wget -q -O a.zip https://microb120.med.upenn.edu/data/export/projects/exchange/neverMind.zip')
system('unzip a.zip')
d <- read_tsv('intSites.tsv')
invisible(file.remove(c('a.zip', 'intSites.tsv', 'readMe.txt'))); unlink('reports', recursive = TRUE)

d$posid <- paste0(d$chromosome, d$strand, d$position)

# Correct earlier shunting of D0 samples to separate patient ids.
d[d$subject == 'p12418-13-CART22-65s.G1',]$cellType <- 'CART22'
d[d$subject == 'p12418-13-CART22-65s.G1',]$subject <- 'p12418-13'

d[d$subject == 'p12418-13-huCART-19.G1',]$cellType <- 'CART19'
d[d$subject == 'p12418-13-huCART-19.G1',]$subject <- 'p12418-13'

# Standardize naming of D0 samples according to vector.
d[grepl('CART19', d$cellType),]$cellType <- 'CART-19'
d[grepl('CART22', d$cellType),]$cellType <- 'CART-22'

# Limit data to Whole blood samples.
d$timePoint <- toupper(d$timePoint)
d$cellType <- toupper(d$cellType)
d0 <- d[d$timePoint == 'D0',]
dx <- d[grepl('WHOLE\\s*BLOOD', d$cellType),]

# Cleanup time points.
dx[grepl('D247', dx$timePoint),]$timePoint <- 'M8'
dx[grepl('D28',  dx$timePoint),]$timePoint <- 'M1'
dx[grepl('M12',  dx$timePoint),]$timePoint <- 'Y1'
dx[grepl('M18',  dx$timePoint),]$timePoint <- 'Y1.5'
dx[grepl('M24',  dx$timePoint),]$timePoint <- 'Y2'
dx[grepl('M36',  dx$timePoint),]$timePoint <- 'Y3'
dx[grepl('M42',  dx$timePoint),]$timePoint <- 'Y3.5'

# Limit samples to those with 50 or more inferred cells.
dx <- group_by(dx, internalSampleID) %>% 
      mutate(totalCells = sum(estAbund)) %>% 
      ungroup() %>%
      filter(totalCells >= 50) %>%
      mutate(timePoint = factor(timePoint, levels = c("M1", "M3", "M6", "M8", "Y1", "Y1.5", "Y2", "Y3.5")))

# Build relative abundance plot data.
o <- bind_rows(lapply(split(dx, dx$internalSampleID), function(x){
       x <- filter(x, abs(nearestFeatureDist) <= 50000) %>%
            filter(relAbund >= 1) %>%
            arrange(desc(relAbund)) %>% 
            dplyr::slice(1:10) %>% 
            select(subject, timePoint, totalCells, internalSampleID, relAbund, nearestFeature, posid)

       x <- bind_rows(x[1,], x)
       x[1,]$relAbund <- 100 - sum(x[2:nrow(x),]$relAbund)
       x[1,]$nearestFeature <- 'LowAbund' 
       x
     })) %>% filter(! is.na(subject))

# Create a table of genes with two or more proximal integrations.
# Clones with integrations near these genes will be colored in the relative abundance plot.
g <- table(o$nearestFeature)
g <- sort(g, decreasing = TRUE)
g <- g[g > 1]

# Create a new column for plotting based on its inclusion in the table above.
o <- bind_rows(lapply(split(o, o$internalSampleID), function(x){
       x$gene <- x$nearestFeature
       x$gene <- ifelse(! x$gene %in% names(g), 'Other', x$gene)
       x
     }))

# Convert nearest gene to vector and define levels so that they are plotted in the prefered order.
o$gene <- factor(o$gene, levels = c('LowAbund', 'Other', names(g)[! names(g) %in% c('LowAbund', 'Other')]))

olabels <- group_by(o, subject, timePoint) %>% summarise(totalCells = totalCells[1]) %>% ungroup()

# Create relative abundance plot.
relAbundPLot <- ggplot(o, aes(timePoint, relAbund, fill = gene)) + 
  theme_bw() +
  scale_fill_manual(values = c('gray90', 'gray60', colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(g)-1))) +
  geom_col(color = 'black', linewidth = 0.1) +
  facet_wrap(subject~., ncol = 2) +
  labs(x = 'Time point', y = ' Percent relative abundance') +
  scale_y_continuous(limits = c(0,120), breaks = seq(0, 100, by = 25)) +
  geom_text(data = olabels, aes(x = timePoint, y = 115, label = totalCells), inherit.aes = FALSE, size = 3) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


# Build overlap table.
options(scipen=999)
r <- bind_rows(lapply(unique(dx$subject), function(p){
       a1 <- subset(d0, subject == p & cellType == 'CART-19') 
       a2 <- subset(d0, subject == p & cellType == 'CART-22')
       b  <- subset(dx, subject == p)
  
       if(nrow(a1)  == 0 | nrow(a2) == 0 | nrow(b) == 0) return(tibble())
  
       # Convert sites into genomic ranges by adding +/- 5 NT to each position.
       a1.g <- makeGRangesFromDataFrame(a1, seqnames.field = 'chromosome', start.field = 'position', end.field = 'position', keep.extra.columns = TRUE) + 5
       a2.g <- makeGRangesFromDataFrame(a2, seqnames.field = 'chromosome', start.field = 'position', end.field = 'position', keep.extra.columns = TRUE) + 5
       b.g  <- makeGRangesFromDataFrame(b,  seqnames.field = 'chromosome', start.field = 'position', end.field = 'position', keep.extra.columns = TRUE) + 5
  
       b_a1 <- findOverlaps(b.g, a1.g)
       b_a2 <- findOverlaps(b.g, a2.g)

       pval <- fisher.test(matrix(c(n_distinct(a1$posid), n_distinct(a2$posid), 
                                  n_distinct(b.g[queryHits(b_a1)]$posid), n_distinct(b.g[queryHits(b_a2)]$posid)),
                                  ncol = 2, byrow = FALSE))$p.val
       
       a1_malt1 <- paste0(n_distinct(subset(a1, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid), ' (',
                          round(n_distinct(subset(a1, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid) / n_distinct(a1$posid), 5), ')')
       
       a2_malt1 <- paste0(n_distinct(subset(a2, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid), ' (',
                          round(n_distinct(subset(a2, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid) / n_distinct(a2$posid), 5), ')')
       
       dx_malt1 <- paste0(n_distinct(subset(b, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid), ' (',
                          round(n_distinct(subset(b, nearestFeature == 'MALT1' & nearestFeatureDist == 0)$posid) / n_distinct(b$posid), 5), ')')
       
       tibble(subject = p, 
              product_CART19_sites = n_distinct(a1$posid),
              product_CART19_MALT1sites = a1_malt1,
              product_CART22_sites = n_distinct(a2$posid),
              product_CART22_MALT1sites = a2_malt1,
              postTransfusionSites = n_distinct(b$posid),
              postTransfusion_MALT1sites = dx_malt1,
              postTransfusion_overlapped_CART19_sites = n_distinct(b.g[queryHits(b_a1)]$posid),
              postTransfusion_overlapped_CART22_sites = n_distinct(b.g[queryHits(b_a2)]$posid),
              pVal = pval)
     }))

openxlsx::write.xlsx(r, 'overlapTable.xlsx')


# Build UCSC track.
#----------------------------------------------------------------------------------------------------------------------------------------------------
# http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&hgt.customText=http://microb120.med.upenn.edu/data/export/UCSC/CART_12418.ucsc
title <- 'CART_12418'
visibility <- 1
posColors  <- c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB")
negColors  <- c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000")
position   <- 'chr18:58671465-58754477'
abundCuts  <- c(3, 10, 20)
trackFileName <- 'CART_12418.ucsc'

# Convert Hex color codes to RGB color codes. 
posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')

# Bin relative abundances for tick mark colors on UCSC browser.
d$cuts <- cut(d$relAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
d$trackName <- ifelse(d$timePoint == 'D0', ifelse(d$subject == 'p12418-13', 'D0_p13', 'D0_not_p13'), ifelse(d$subject == 'p12418-13', 'Dx_p13', 'Dx_not_p13'))
d$siteLabel <- gsub('\\s', '', paste0(d$subject, '_', d$cellType, '_', d$timePoint))

d$score <- 0
d$color <- ifelse(d$strand == '+', posColors[d$cuts], negColors[d$cuts])

invisible(lapply(split(d, d$trackName), function(x){
  trackFileName <- paste0(x$trackName[1], '.ucsc')
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s", x$trackName[1], x$trackName[1], visibility)
  write(trackHead, file = trackFileName, append = FALSE)
  write.table(x[, c('chromosome', 'position', 'position', 'siteLabel', 'score', 'strand', 'position', 'position', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = trackFileName, append = TRUE, quote = FALSE)
}))

system('cat D0_not_p13.ucsc Dx_not_p13.ucsc D0_p13.ucsc Dx_p13.ucsc > CART_12418.ucsc')
invisible(file.remove(list.files(pattern = 'p13')))
# scp CART_12418.ucsc microb120.med.upenn.edu:/media/lorax/data/export/UCSC/CART_12418.ucsc
system('gzip CART_12418.ucsc')


# ScanStatistics
#----------------------------------------------------------------------------------------------------------------------------------------------------

d0$start <- d0$position
d0$end <- d0$position

dx$start <- dx$position
dx$end <- dx$position

# Create GRanges for early and late sites.
ag <- makeGRangesFromDataFrame(d0, keep.extra.columns = TRUE)
bg <- makeGRangesFromDataFrame(dx, keep.extra.columns = TRUE)

# Remove duplicate sites from each group across all subjects.
ag <- ag[! duplicated(ag$posid)]
bg <- bg[! duplicated(bg$posid)]

# Interpret value1 and value2 from the clusterSource flag.
scanStatsDF <- data.frame(scanStats(ag, bg))
scanStatsDF <- bind_rows(lapply(split(scanStatsDF, 1:nrow(scanStatsDF)), function(x){
  if(x$clusterSource == 'A'){
    x$earlySites = max(c(x$value1, x$value2))
    x$lateSites = min(c(x$value1, x$value2))
  } else {
    x$earlySites = min(c(x$value1, x$value2))
    x$lateSites = max(c(x$value1, x$value2))
  }
  
  x
})) %>% select(seqnames, start, end, width, earlySites, lateSites, target.min)

# Clean up scanStats result data.frame.
names(scanStatsDF) <- c('chromosome', 'start', 'end', 'width', 'earlySites', 'lateSites', 'target.min')
scanStatsDF$target.min <- round(scanStatsDF$target.min, 3)

# Create a data.frame corresponding to the data used for the scanStats analysis to draw nearest genes from.
g <- data.frame(c(ag, bg))

# Look within each cluster and tally the number of nearest genes for final report.
scanStatsDF <- bind_rows(lapply(split(scanStatsDF, 1:nrow(scanStatsDF)), function(x){
  tab <- sort(table(subset(g, seqnames == x$chromosome & start >= x$start & start <= x$end)$nearestFeature), decreasing = TRUE)
  x$nearestGenes <- paste0(paste0(names(tab), ' (', tab, ')'), collapse = ', ')
  x <- mutate(x, cluster = paste0(x$chromosome, ':', x$start, '-', x$end), .before = 'chromosome')
  #x$url <- paste0(config$scanStatsBEDfileURL, "&position=", x$chromosome, "%3A", x$start-10, "-", x$end+10)
  select(x, -chromosome, start, end)
}))

openxlsx::write.xlsx(scanStatsDF, 'scanStats.xlsx')
