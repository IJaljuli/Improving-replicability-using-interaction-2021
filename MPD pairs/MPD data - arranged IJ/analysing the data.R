library(dplyr)
source('/Users/Iman/Documents/A cross laboratory investigation of timing endophenotypes in mouse behavior - original/Mollys data/analysis of multiple sets togather/general_functions.R')
setwd('~/Library/Mobile Documents/com~apple~CloudDocs/A cross laboratory investigation of timing endophenotypes in mouse behavior/MPD pairs/full sheets June 19')
compute.df.sater <- function(s2int,df.int,s2pooled,n1,n2){
  inv.v <- (((1/n1+1/n2)*s2pooled)^2*1/(n1+n2-2)+(2*s2int)^2*1/df.int)/
    (((1/n1+1/n2)*s2pooled+2*s2int)^2)
  df <- 1/inv.v
  return(df)
}

file.names <- list.files()

title <- lapply( strsplit(file.names , split = 'Lab_') , FUN = function(x)x[2] )
title <- unlist(title)

file.names <- file.names[!is.na(title)]

title <- title [!is.na(title)]

gxl.info <- read.csv('gxl_estimates.csv',header = T , stringsAsFactors = F)
gxl.fac <- gxl.info$gamma.estimator
gxl.L <- gxl.info$labs.int
gxl.S <- gxl.info$strains.int
# grip.strength.trans <- function(x) x^(1/3)
# percent.center.time <- function(x) log((x+0.1)/(100.1-x))
names(gxl.fac) <- names(gxl.S) <- names(gxl.L) <- gxl.info$File.name
keep.strains <- c("BALB/cJ","BTBR T+ tf/J","C3H/HeJ","C57BL/6J","DBA/2J","SWR/J")

# file.name = file.names[length(file.names)]
for (file.name in file.names){
  dd <- read.csv(file = file.name , header = T , stringsAsFactors = F)
  if (colnames(dd)[1] =='X') dd <- dd[,-1]
  variable.name <- colnames(dd)[2]
  'value' ->  colnames(dd)[2]
  if (file.name %in% file.names[3:4] ){
    dd$value <- dd$value^(1/3)
  }
  if (file.name %in% file.names[c(7,9)] ){
    dd$value <- log((dd$value+0.1)/(100.1-dd$value))
  }

  if (file.name == file.names[10] ){
    dd$value <- dd$value
  }

  # strain.list <- data.frame(strain = unique(dd$strain) , include = rep('YES' , strains))
  # write.csv(x = strain.list , file = paste0('Strains_' , file.name))
  dd <- dd %>% filter(!is.na(value), strain %in% keep.strains )
  
  dd.means <- dd %>% 
    group_by(strain, group) %>% 
    summarise(mean = mean(value,na.rm = T) , sd = sd(value,na.rm = T) , n = n() )
  strains <- length(unique(dd$strain)) 
  
  pairs.strainsT <- t(combn(unique(dd$strain) , m = 2)) 
  pairs.treatsT <- matrix(unique(dd$group),ncol = 2)
  
  pairs.strains <- expand.grid( strain1 = NA, strain2 = NA,
                                pair.num = 1:nrow(pairs.strainsT), 
                                treatment1 = unique(dd.means$group),
                                stringsAsFactors = F) %>%
    mutate( treatment2 = treatment1 ) 
  pairs.strains[,'strain1'] <- pairs.strainsT[pairs.strains$pair.num,1]
  pairs.strains[,'strain2'] <- pairs.strainsT[pairs.strains$pair.num,2]
  
  pairs.treats <- expand.grid( treatment1 = NA, treatment2 = NA,
                               treat.num = 1:nrow(pairs.treatsT), 
                                strain1 = unique(dd.means$strain), 
                               stringsAsFactors = F) %>%
    mutate( strain2 = strain1 ) 
  pairs.treats[,'treatment1'] <- pairs.treatsT[pairs.treats$treat.num,1]
  pairs.treats[,'treatment2'] <- pairs.treatsT[pairs.treats$treat.num,2]
  
  
  pairs <- rbind(pairs.strains[,-3] , pairs.treats[colnames(pairs.strains)[-3]])
  pairs <- pairs %>% 
    mutate(same.strain = (strain1 == strain2 ),
           same.treatment = (treatment1 == treatment2)) %>%
    filter((same.treatment + same.strain) == 1)
  pairs <- droplevels(pairs)
  pairs <- pairs[,c("strain1", "strain2", "treatment1", "treatment2")]

  pairs <- merge(x = pairs , y = dd.means , by.x = c('strain1', 'treatment1') , by.y = c('strain', 'group'))
  colnames(pairs)[5:7] <- paste0('strain1.' , colnames(pairs)[5:7]) 
  pairs <- merge(x = pairs , y = dd.means , by.x = c('strain2', 'treatment2') , by.y = c('strain', 'group'))
  colnames(pairs)[8:10] <- paste0('strain2.' , colnames(pairs)[8:10]) 
  pairs <- pairs %>% 
    mutate( diff = abs(strain1.mean - strain2.mean )) %>%
    mutate(s2.pooled = ( (strain1.n-1)*strain1.sd^2 + (strain2.n-1)*strain2.sd^2) / (strain1.n + strain2.n - 2)) %>%
    mutate( SE = sqrt(  s2.pooled*(1/strain1.n + 1/strain2.n )  ),
            SE.adjusted = sqrt(  s2.pooled*(1/strain1.n + 1/strain2.n + 2*gxl.fac[file.name] )  )) %>%
    mutate(satter.df = compute.df.sater(s2int = 2*s2.pooled*gxl.fac[file.name]  ,
                                        df.int = (gxl.L[file.name]-1)*(gxl.S[file.name]-1),
                                        s2pooled = s2.pooled,  n1 = strain1.n ,  n2 = strain2.n) ) %>% 
    mutate(gxl.adjusted.pv = 2*pt( diff/SE.adjusted  , df = satter.df , lower.tail = F) ) %>%
    mutate(common.pv = 2*pt( diff/SE  , df = (strain1.n + strain2.n - 2) , lower.tail = F) )
  
  pairs$FDR.common.pv <- p.adjust( pairs$common.pv , method = 'BH' ) 
  pairs$FDR.gxl.adjusted.pv <- p.adjust( pairs$gxl.adjusted.pv , method = 'BH' ) 
  ## saving in a comfort form , stupidly
  pairs.save <- pairs[,c("strain1" , "strain2", "gxl.adjusted.pv")]
  pairs.temp <- pairs[,c("strain1" , "strain2", "FDR.gxl.adjusted.pv")]
  colnames(pairs.temp)[3] <- colnames(pairs.save)[3] <- 'pvalue'
  pairs.save$FDR = F ; pairs.save$gxl = T
  pairs.temp$FDR = T ; pairs.temp$gxl = T
  pairs.save <- rbind(pairs.save , pairs.temp)
  pairs.temp <- pairs[,c("strain1" , "strain2", "common.pv")]
  colnames(pairs.temp)[3] <- 'pvalue'
  pairs.temp$FDR = F ; pairs.temp$gxl = F
  pairs.save <- rbind(pairs.save , pairs.temp)
  pairs.temp <- pairs[,c("strain1" , "strain2", "FDR.common.pv")]
  colnames(pairs.temp)[3] <- 'pvalue'
  pairs.temp$FDR = T ; pairs.temp$gxl = F
  pairs.save <- rbind(pairs.save , pairs.temp)
  # done 
  
  
  write.csv(x = pairs[, names(pairs)[c(3,1,4,2,11,17,16)]] , file = paste0('Analysed_' , title[which(file.names==file.name)]) , row.names = F)
  
  # write.csv(x = pairs[, names(pairs)[c(2,1,9,14:17)]] , file = paste0('Analysed_' , title[which(file.names==file.name)]) , row.names = F)
  write.csv(x = pairs.save , file = paste0('Heatmap_' , title[which(file.names==file.name)]) , row.names = F)
}



