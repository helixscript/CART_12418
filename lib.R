scanStats <- function(gr1, gr2, gr1.label = 'A', gr2.label = 'B', kvals = '15L:30L', nperm = 100,
                      cutpt.tail.expr = 'critVal.target(k, n, target = 5, posdiff = x)',
                      cutpt.filter.expr = 'as.double(apply(x, 2, median, na.rm = TRUE))'
                      #cutpt.filter.expr = paste0('apply(x, 2, quantile, probs = ', p, ', na.rm = TRUE)')
){
  library(GenomicRanges)
  library(geneRxCluster)
  
  gr1$clusterSource <- TRUE
  gr2$clusterSource <- FALSE
  
  gr3 <- c(gr1, gr2)
  
  df <- GenomicRanges::as.data.frame(gr3)
  df <- df[order(df$seqnames, df$start),]
  df$seqnames <- droplevels(df$seqnames)
  row.names(df) <- NULL
  df <- df[, c('seqnames', 'start', 'clusterSource')]
  
  comm <- paste0('gRxCluster(df$seqnames, df$start, df$clusterSource, ',  kvals, ', nperm=', nperm, ', cutpt.tail.expr=', cutpt.tail.expr, ', cutpt.filter.expr=', cutpt.filter.expr, ')')
  scan <- eval(parse(text = comm))
  
  if(length(scan)==0) return(GRanges())
  
  gr1.overlap <- GenomicRanges::findOverlaps(gr1, scan, ignore.strand=TRUE)
  gr2.overlap <- GenomicRanges::findOverlaps(gr2, scan, ignore.strand=TRUE)
  
  if(length(gr1.overlap) > 0 & length(gr2.overlap) == 0){
    scan$clusterSource <- gr1.label
    return(scan)
  }
  
  if(length(gr1.overlap) == 0 & length(gr2.overlap) > 0){
    scan$clusterSource <- gr2.label
    return(scan)
  }
  
  gr1.overlap <- stats::aggregate(queryHits ~ subjectHits, data=as.data.frame(gr1.overlap), FUN=length)
  gr2.overlap <- stats::aggregate(queryHits ~ subjectHits, data=as.data.frame(gr2.overlap), FUN=length)
  
  gr1.list <- list()
  gr1.list[gr1.overlap[,1]] <- gr1.overlap[,2]
  
  gr2.list <- list()
  gr2.list[gr2.overlap[,1]] <- gr2.overlap[,2]
  
  # index error protection
  gr1.list[seq(length(gr1.list)+1, (length(gr1.list) + length(scan)-length(gr1.list) + 1))] <- 0
  gr2.list[seq(length(gr2.list)+1, (length(gr2.list) + length(scan)-length(gr2.list) + 1))] <- 0
  
  scan$clusterSource <- '?'
  
  for(i in 1:length(scan))
  {
    if(! is.null(gr1.list[[i]]) & is.null(gr2.list[[i]])){
      scan[i]$clusterSource <- gr1.label
    } else if (is.null(gr1.list[[i]]) & ! is.null(gr2.list[[i]])){
      scan[i]$clusterSource <- gr2.label
    } else if (gr1.list[[i]] > gr2.list[[i]]){
      scan[i]$clusterSource <- gr1.label
    } else {
      scan[i]$clusterSource <- gr2.label
    }
  }
  
  scan
}



tuneScanStatistics <- function(A, B, CPUs=40, 
                               minWindow=5, maxWindow=100, windowStep=5, 
                               minProb=0.75, maxProb=0.95, probStep=0.025, 
                               minTarget=5, maxTarget=100, targetStep=5){
  
  t <- do.call(rbind, lapply(minWindow:(maxWindow - minWindow), function(i){
    r <- data.frame()
    
    for(i2 in (minWindow+windowStep):maxWindow){
      if(i2-i>=windowStep){
        for(p in seq(minProb, maxProb, by=probStep)){
          for(t in seq(minTarget, maxTarget, by=targetStep)){
            kvals <- paste0(i, 'L:', i2, 'L')
            cutpt.filter.expr <- paste0('apply(x, 2, quantile, probs = ', p, ', na.rm = TRUE)')
            cutpt.tail.expr   <- paste0('critVal.target(k, n, target = ', t, ', posdiff = x)')
            
            r <- rbind(r, data.frame(kvals=kvals,
                                     cutpt.filter.expr=cutpt.filter.expr,
                                     cutpt.tail.expr=cutpt.tail.expr))
          }
        }
      }
    }
    
    r
  }))
  
  t$s <- ceiling(seq_along(1:nrow(t))/(nrow(t)/CPUs))
  save(list = ls(all.names = TRUE), file='tuneData.RData', envir = environment())
  
  cluster <- parallel::makeCluster(CPUs)
  t2 <- do.call(rbind, parallel::parLapply(cluster, split(t, t$s), function(x){
    #t2 <- do.call(rbind, lapply(split(t, t$s), function(x){  
    load('tuneData.RData')
    
    do.call(rbind, lapply(1:nrow(x), function(i){
      message('chunk ', x$s[1], ' parameter row ', i)
      
      tryCatch({
        scan <- gt23::scanStats(A, B,
                                kvals=x[i,]$kvals, 
                                nperm='1000',
                                cutpt.filter.expr=x[i,]$cutpt.filter.expr,
                                cutpt.tail.expr=x[i,]$cutpt.tail.expr)
        
        data.frame(kvals=x[i,]$kvals,
                   cutpt.filter.expr=x[i,]$cutpt.filter.expr,
                   cutpt.tail.expr=x[i,]$cutpt.tail.expr,
                   FDR=gRxSummary(scan)$FDR,
                   totalClusters=gRxSummary(scan)$Clusters_Discovered,
                   Aclusters=length(subset(scan, clusterSource=='A')),
                   Bclusters=length(subset(scan, clusterSource=='B')),
                   avgClusterWidth=(sum(width(scan)) / length(scan)))
      }, 
      error = function(cond) {
        return(data.frame())
      })
    }))
  }))
  stopCluster(cluster)
  t2
}
