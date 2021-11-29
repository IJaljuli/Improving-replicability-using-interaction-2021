library(dplyr)
library(ggplot2)
library(XLConnect)
library(reshape2)
setwd(dir = '~/Dropbox/A cross laboratory investigation of timing endophenotypes in mouse behavior/MPD pairs/GxL app output/Pairwise Comparisons')
files.list <- sort(list.files())[-7] #[c(4:8,14:18)] 
# csv.files.list <- files.list[3:4]
# files.list <- files.list[-c(3:4)]
tables_names <- strsplit(files.list , split = 'FDR ')
tables_names <- unlist(lapply( tables_names , function(x) {x[2]})) %>% 
  gsub(pattern = ' ' ,replacement = '_',x = . ) %>% 
  strsplit( . , split = '.xlsx') %>% unlist(.)
tables_names[is.na(tables_names)] <- files.list[is.na(tables_names)] %>% strsplit( . , split = '.xlsx') %>% unlist(.)
measures <- tables_names
remove.strain <- c()
pdf.path <- '~/Dropbox/A cross laboratory investigation of timing endophenotypes in mouse behavior/MPD pairs/GxL app output/removing_'
keep.strains <- c('C3H/HeJ' , 'C57BL/6J' , 'DBA/2J' , 'FVB/NJ' , 'SJL/J') #, 'A/J' , 'CBA/J')
keep.strains <- c( keep.strains ,  paste0( keep.strains ,'fluoxetine') , paste0( keep.strains ,'control') )
# tables_names <- paste( tables_names , 1:length(tables_names ) ,sep = '_')
# I will only take the sets without FDR and perform by myself. 
while( length(remove.strain) < 4 ){
  for ( i in 1:length(files.list)){
    # with.FDR <- (length(unlist(strsplit( files.list[i] ,split = ' No ')))==1)
    with.FDR <- F
    
    read.sheet <- 'readWorksheetFromFile(file = files.list[i] , sheet = 1,header = T)'
    eval(parse ( text = paste ( tables_names[i] ,'<-' , read.sheet )))
    
    col.naming <- paste('colnames(',  tables_names[i] ,') <- as.character(', tables_names[i],'[1,])' ) 
    eval(parse ( text = col.naming ))
    
    col.naming <- paste('colnames(',  tables_names[i] , ')[c(3,6)] <- c( "F" , "T")')
    eval(parse ( text = col.naming ))
    eval(parse ( text = paste ( tables_names[i] ,'<-' , tables_names[i] ,'[-1,c(1:3,6)]' )))
    
    fdr.ind.col <- paste0( tables_names[i] ,' <- ', tables_names[i],
                           ' %>% mutate( FDR = ',with.FDR,')')
    eval(parse ( text = fdr.ind.col ))
    
    melting = paste0 ( tables_names[i], '<- melt(data =  ',tables_names[i],
                       ',variable.name = "GxL" , id.vars = c("Comparison","Difference","FDR"),
                       measure.vars = c("T","F") , value.name = "pvalue"  ) %>% 
                       mutate( measure = measures[i] , GxL = as.logical(GxL) )' )

    eval(parse ( text = melting ))
    
    strain.names <- paste0( 'strains <- strsplit(', tables_names[i] ,'$Comparison , split = " - ")')
    eval(parse ( text = strain.names ))
    strains <- t(unname(as.data.frame.list(strains)))
    colnames(strains) <- c('strain1' , 'strain2')
    strain.names <- paste0( tables_names[i] , ' <- cbind(  strains ,', tables_names[i] ,' )')
    eval(parse ( text = strain.names ))
    eval(parse ( text = 
                   paste0( tables_names[i] , ' <-' , tables_names[i] , 
                           '[, c("measure" , "strain1","strain2","Comparison","FDR","GxL","Difference","pvalue")]')))
    
    pval.text <-  paste0( tables_names[i] , ' <-' , tables_names[i] , 
                          '%>% mutate(pvalue = as.numeric( pvalue))%>% 
                          mutate(pvalue = replace( pvalue , is.na(pvalue) , 0 )) %>%
                          filter( (strain1 %in% keep.strains)&(strain2 %in% keep.strains) ) ' )
    eval(parse ( text = pval.text ) ) 
    
    fdr.text <-  paste0( 'temp <-' , tables_names[i] )
    eval(parse ( text = fdr.text ) ) 
    
    fdr.text <- paste0( 'temp$pvalue[ temp$GxL ] <- p.adjust(p = temp$pvalue[ temp$GxL ] , method = "BH" )' ) 
    eval(parse ( text = fdr.text ) ) 
    
    fdr.text <- paste0( 'temp$pvalue[ !temp$GxL ] <- p.adjust(p = temp$pvalue[ !temp$GxL ] , method = "BH" )' ) 
    eval(parse ( text = fdr.text ) ) 
    
    eval(parse ( text = paste0( 'temp$FDR <- T ' ) ) ) 
    
    eval(parse ( text = paste ( tables_names[i] , '<- rbind(', tables_names[i] , ' , temp )' ) ))
    
  }
  binding.tables <- 
    paste0( 'all.results.tables <- rbind(' , paste(tables_names,collapse = ' , '),')' )
  eval(parse ( text =  binding.tables ))
  
  strain.list <- unique(c(as.character( all.results.tables$strain1 ) , 
                          as.character( all.results.tables$strain2)))
  specail_genotypes_count <- rep(0 , length(strain.list))
  names(specail_genotypes_count) = strain.list
  
  
  
  # all.results.tables$pvalue  <- as.numeric(all.results.tables$pvalue) 
  
  all.results.tables <- all.results.tables %>% 
    # mutate( pvalue = replace( pvalue , is.na(pvalue), 0  ) ) %>% 
    mutate( significant = (pvalue < 0.05) )
  
  
  temp <- all.results.tables ; temp[,2:3] <- temp[,3:2]
  colnames(temp)[2:3] <- colnames(all.results.tables)[2:3]
  results.tables.not.yellow <- rbind( all.results.tables , temp )
  results.tables.not.yellow <- results.tables.not.yellow %>% 
    mutate( yellow =  (GxL&FDR&(pvalue<0.05)) , white_blank = (GxL&(pvalue>=0.05)) ) %>%
    group_by( strain1 ,GxL , FDR ) %>%
    summarise(  not_yellow = sum( !yellow ) )
  results.tables.not.yellow <- results.tables.not.yellow[ order(results.tables.not.yellow$not_yellow) , ]
  
  
  pdf(file =  paste0(pdf.path , length(remove.strain) ,'_strains.pdf' ) )
  measures <- unique(measures)
  
  for( i in 1:length(measures)){
    measure.results.table <- all.results.tables %>% 
      filter( measure == measures[i])
    
    measure.results.table2 <- measure.results.table[,c(1,3,2,4:9)]
    colnames(measure.results.table2)[2:3] <-colnames(measure.results.table2)[3:2] 
    measure.results.table <- rbind(measure.results.table,measure.results.table2)
    
    temp <- measure.results.table[,c(1,2,2,4:9)]
    colnames(temp)[3]<- 'strain2'
    temp <- temp[!duplicated( temp$strain1) , ] 
    temp$significant <- F
    temp$pvalue <- NA
    measure.results.table <- rbind(measure.results.table , temp)
    
    
    plotiing.data <- measure.results.table %>% group_by( strain1 , strain2 ) %>%
      summarise( with.gxl =  sum( significant*GxL ,na.rm = T), # coloring
                 without.gxl = sum( significant*(!GxL),na.rm = T  ) ,
                 #               pvalue.missing = all(is.na(pvalue)&(strain1!=strain2)) ) %>% # shape
                 pvalue.missing = all(is.na(pvalue))) %>% # shape
      
      mutate( with.gxl = replace( with.gxl , pvalue.missing , 3 ) ,
              without.gxl = replace( without.gxl , pvalue.missing , 3 )  )
    # plotiing.data <- measure.results.table %>% group_by( strain1 , strain2 ,GxL ) %>% 
    #   summarise( by.gxl =  sum( significant), n = n(),
    #              pvalue.missing = all(is.na(pvalue)&(strain1!=strain2)) ) %>% # shape
    #   mutate( with.gxl = replace( with.gxl , pvalue.missing , 3 ) , 
    #           without.gxl = replace( without.gxl , pvalue.missing , 3 )  )
    
    
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
    specail_genotypes_count[ sellect.cells ] <- specail_genotypes_count[sellect.cells] + specail_gen$special   
    
    strain.levels <- sort(unique(plotiing.data$strain1))
    plotiing.data$strain1 <- factor(plotiing.data$strain1 , levels = strain.levels)
    plotiing.data$strain2 <- factor(plotiing.data$strain2 , levels = strain.levels)
    
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
      ggtitle(measures[i]) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gp)
  }
  dev.off()
  remove.strain <-  c( remove.strain , names(specail_genotypes_count)[which.min(specail_genotypes_count)] )
  print( list( names(specail_genotypes_count)[which.min(specail_genotypes_count)]  , minimal.count = min(specail_genotypes_count)))
}

