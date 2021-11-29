setwd('~/Dropbox/A cross laboratory investigation of timing endophenotypes in mouse behavior/MPD pairs/full sheets June 19')

file.names <- list.files()

title <- lapply( strsplit(file.names , split = 'Heatmap_') , FUN = function(x)x[2] )
title <- unlist(title)

file.names <- file.names[!is.na(title)]

title <- title [!is.na(title)]
tables_names <- title



# pdf(file =  paste0(pdf.path , length(remove.strain) ,'_strains.pdf' ) )
pdf(file =  'hetamap_Crabbe4_grip_strength_Males.pdf ')

for( i in 1:length(tables_names)){
  
  measure.results.table <- read.csv(file.names[i],header = T,stringsAsFactors = F)
  
  measure.results.table <- measure.results.table %>% mutate( significant = (pvalue < 0.05) )
  
  plotiing.data <- measure.results.table %>% group_by( strain1 , strain2 ) %>%
    summarise( with.gxl =  sum( significant*gxl ,na.rm = T), # coloring
               without.gxl = sum( significant*(!gxl),na.rm = T  ) ,
               pvalue.missing = F) %>% # shape
    mutate( with.gxl = replace( with.gxl , pvalue.missing , 3 ) ,
            without.gxl = replace( without.gxl , pvalue.missing , 3 )  )
  
  plotiing.data$with.gxl <- as.character(plotiing.data$with.gxl)
  plotiing.data$with.gxl[plotiing.data$with.gxl=='0'] = 'Not significant'
  plotiing.data$with.gxl[plotiing.data$with.gxl=='1'] = 'Significant without FDR only'
  plotiing.data$with.gxl[plotiing.data$with.gxl=='2'] = 'Significant with FDR'
  plotiing.data$with.gxl[plotiing.data$with.gxl=='3'] = 'Missing'
  plotiing.data$without.gxl <- as.character(plotiing.data$without.gxl) 
  plotiing.data$without.gxl[plotiing.data$without.gxl=='0'] = 'Not significant'
  plotiing.data$without.gxl[plotiing.data$without.gxl=='1'] = 'Significant without FDR only'
  plotiing.data$without.gxl[plotiing.data$without.gxl=='2'] = 'Significant with FDR'
  plotiing.data$without.gxl[plotiing.data$without.gxl=='3'] = 'Missing' 
  
  
  specail_gen <- plotiing.data %>%
    mutate( yellow_dark =  ((with.gxl == 'Significant with FDR' )&( without.gxl == 'Significant with FDR' )) , 
            blank_white = (( with.gxl == 'Not significant' ) & ( without.gxl == 'Not significant' )),  
            missing.comp = ((with.gxl == 'Missing' )&(without.gxl == 'Missing' )) ) %>%
    group_by( strain1 ) %>% summarise( special = sum( !(yellow_dark | blank_white | missing.comp) ))
  sellect.cells <-  as.character( specail_gen$strain1 )
  #specail_genotypes_count[ sellect.cells ] <- specail_genotypes_count[sellect.cells] + specail_gen$special   
  strain.levels <- unique(c(plotiing.data$strain1,plotiing.data$strain2))
  strain.levels.sorted <- sort(strain.levels) #c(strain.levels[c(4,6,7,11,12,18)] , strain.levels[-c(4,6,7,11,12,18)])
  strain.levels1 <- strain.levels.sorted[strain.levels.sorted%in%unique(plotiing.data$strain1)]
  strain.levels2 <- strain.levels.sorted[strain.levels.sorted%in%unique(plotiing.data$strain2)]
  plotiing.data$strain1 <- factor(plotiing.data$strain1 , levels = strain.levels1)
  plotiing.data$strain2 <- factor(plotiing.data$strain2 , levels = strain.levels2)
  
  plotiing.data <- as.data.frame(plotiing.data)
  plotiing.data$with.gxl <- factor(plotiing.data$with.gxl , levels = c('Not significant' , 'Significant without FDR only' , 'Significant with FDR' , 'Missing') )
  plotiing.data$without.gxl <- factor(plotiing.data$without.gxl , levels = c('Not significant' , 'Significant without FDR only' , 'Significant with FDR' , 'Missing') )
  
  
  gp <-  ggplot( plotiing.data , aes(y = strain2, x = strain1 )) +        ## global aes
    geom_tile(aes(fill =   with.gxl  )) + ## to get the rect filled
    scale_fill_manual(breaks= c('Not significant' , 'Significant without FDR only' , 'Significant with FDR' , 'Missing') ,
                      values= c('white' , 'red', 'yellow' , 'grey'), drop=FALSE)+
    geom_point( aes( shape =  without.gxl ))  +
    scale_shape_manual(breaks= c('Not significant' , 'Significant without FDR only' , 'Significant with FDR' , 'Missing') ,
                       values = c(32,2,17,4) , drop=FALSE) +
    ggtitle(title[i]) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(gp)
}
dev.off()



