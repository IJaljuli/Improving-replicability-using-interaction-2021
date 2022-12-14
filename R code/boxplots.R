### First run the script "full analysis code combined.R" 

file.names <- c("Lab_Crabbe4_grip_strength_Females.csv",
                "Lab_Crabbe4_grip_strength_Males.csv", 
                "Lab_Crowley2_body_weight_Males.csv",
                "Lab_Tarantino2_DistanceTraveled_Females.csv",
                "Lab_Tarantino2_PercentCenter_Females.csv",
                "Lab_Tordoff3_body_weight_Females.csv",
                "Lab_Tordoff3_body_weight_Males.csv",
                "Lab_Wiltshire2_DistanceTraveled_Males.csv" ,
                "Lab_Wiltshire2_PercentCenter_Males.csv",
                "Lab_Wiltshire2_TST_Males.csv")

dd.all <-  NULL

for (file.name in file.names){
  dd <- read_csv(file = paste0(wd,file.name),skip_empty_rows = T)
  if( ! 'sex' %in% names(dd) ){
    if (str_count(file.name, pattern = 'Female')>0){
      dd <- dd %>% mutate(sex = 'Females')
    }else{
      if (str_count(file.name, pattern = 'Male')>0){
        dd <- dd %>% mutate(sex = 'Males')
      }else{
        dd <- dd %>% mutate(sex =  'unspecified')
      } 
    }
  }
  
  if( ! 'treatment' %in% names(dd) ){
    if( 'group'  %in% names(dd)){
      dd <- dd %>% mutate(treatment = group) %>%
        dplyr::select(-c(group))
    }else{
      dd <- dd %>% mutate(treatment = 'Control')
    }
  }
  
  if( ! 'lab' %in% names(dd) ){
    ll <- strsplit(file.name,split = '_') %>% unlist(.) %>% .[[2]]
    dd <- dd %>% mutate(lab=ll)
  }
  if(any( c('X1') %in% names(dd))) { dd <-  dd %>% dplyr::select(-c('X1')) }
  if(any( c( '...1') %in% names(dd))) { dd <-  dd %>% dplyr::select(-c('...1')) }
  dd <- dd %>% 
    pivot_longer(cols = -c('sex','strain','lab','treatment'), names_to = 'y.name', values_to = 'y' ) %>% 
    mutate(y.name = recode(  file.name, 
                             'Lab_Wiltshire2_PercentCenter_Males.csv' = 'OFTsmall_centertime_10m_percent', 
                             'Lab_Wiltshire2_DistanceTraveled_Males.csv' = 'OFTsmall_dist_10m_sec',
                             'Lab_Tordoff3_body_weight_Males.csv' = 'Body.Weight' ,
                             'Lab_Tordoff3_body_weight_Females.csv' = 'Body.Weight', 
                             'Lab_Tarantino2_PercentCenter_Females.csv' = 'OFTlarge_centertime_10m_percent', 
                             'Lab_Tarantino2_DistanceTraveled_Females.csv' = 'OFTlarge_dist_10m_sec' ,
                             'Lab_Crowley2_body_weight_Males.csv' =  'Body.Weight'  ,
                             'Lab_Crabbe4_grip_strength_Males.csv' = 'grip.avg' , 
                             'Lab_Crabbe4_grip_strength_Females.csv' ='grip.avg'  ,
                             "Lab_Wiltshire2_TST_Males.csv" = 	'tst_7min_percent') ) %>%
    unite("y.name", c(y.name, treatment), sep = ':', remove = T) %>%
    dplyr::select(y.name,sex,strain,lab,y) %>% 
    mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR T<+> Itpr3<tf>/J'='BTBR',
                             'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                             'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                             'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                             'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J')) # %>%
  # filter(!(strain %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') )) # "C3H/HeJ" 
  
  if( dd$lab[1]=="Lab_Crabbe4_grip_strength_Females.csv"){
    dd <- dd %>%  filter(!is.na(y)) %>%  filter(y>0)
  }
  if( file.name=="Lab_Wiltshire2_TST_Males.csv"){
    dd <- dd %>%  mutate(y = y/100)
  }
  
  tds <- dd %>% group_split(sex,y.name)
  tds.keep <- lapply(tds, FUN = function(x){ x %>% filter(!is.na(y)) } ) %>%
    lapply(., FUN = function(x){ n_distinct(x$strain)>1} ) %>%
    unlist(.) %>% which(.)
  dd <-   lapply( as.list(tds[tds.keep]),
                  FUN = function(x) clean_dat(x, desig.variable = 'y' )  ) 
  
  dd <-  do.call(rbind, dd)
  
  if ( gxl.transformations[file.name] == "x^(1/3)" ){
    dd$y <- dd$y^(1/3)
  }else{
    dd.y.name <-  strsplit(dd$y.name[1], split = ':') %>% unlist(.)
    dd.y.name <- dd.y.name[1]
    
    if( (gxl.transformations[file.name] ==  "log((x+0.1)/(100.1-x))")|
        ( dd.y.name %in% c('tst_6min_percent', 'tst_7min_percent') )){
      t <- dd$y[dd$y>0] ; t <- 0.5*min(t, na.rm = T)
      dd <- dd %>%  mutate(y = log((y+t)/(1+t-y)  ))
    }
    
  }
  dd <- dd %>% mutate(sex = recode(sex, 'Males'='Male', 'Females'='Female'))
  
  dd <- dd %>%
    separate(col = y.name, into = c('y.name', 'treatment'), sep = ':' ) %>%
    mutate(treatment = recode(treatment, 'fluoxetine'='Fluoxetine', 'control'='Control')) %>%
    mutate(lab = as.character(lab))
  
  dd.all <- rbind(dd, dd.all)
  
  # dd.means <- dd %>%
  #   group_by(strain, lab, sex, y.name, treatment) %>% 
  #   summarise(mean = mean(y,na.rm = T) , sd = sd(y,na.rm = T) , n = n() )
  # 
  # 
  # pairs.right <- expand_grid( strain1 = as.character(unique( dd$strain)), 
  #                             y.name = unique(dd$y.name), 
  #                             treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
  #   filter( treatment1 < treatment2 , !is.na(strain1) )
  # 
  # pairs.right <- dd.means %>%
  #   rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment1','y.name', 'mean1R', 'sd1R', 'n1R')) %>%
  #   right_join(pairs.right)
  # 
  # pairs.right <- dd.means %>%
  #   rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment2', 'y.name', 'mean2R', 'sd2R', 'n2R')) %>%
  #   right_join(pairs.right)
  # 
  # 
  # pairs.left <- expand_grid( strain2 = as.character(unique( dd$strain)), 
  #                            y.name = unique(dd$y.name), 
  #                            treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
  #   filter( treatment1 < treatment2 , !is.na(strain2) )
  # 
  # pairs.left <- dd.means %>%
  #   rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment1', 'y.name', 'mean1L', 'sd1L', 'n1L')) %>%
  #   right_join(pairs.left)
  # 
  # pairs.left <- dd.means %>%
  #   rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment2', 'y.name','mean2L', 'sd2L', 'n2L')) %>%
  #   right_join(pairs.left)
  # 
  # pairs.right <- pairs.right %>% 
  #   mutate(diff.R = ( mean1R - mean2R ) )
  # 
  # pairs.left <- pairs.left %>% 
  #   mutate(diff.L = ( mean1L - mean2L ) )
  # 
  # pairs.all.WILTSHIRE2[[file.name]] <- full_join(pairs.right, pairs.left) %>%
  #   filter(strain1 > strain2) %>% 
  #   rename('measure.name'='y.name' ) %>%
  #   filter(!((strain1 %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') )&(strain2 %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') )) )
  # 
}

dd.all <- dd.all %>% rename('measure.name' = 'y.name') %>% select(strain, lab, y, measure.name, treatment, sex)

dd.all <- rbind(dd.all, all_data_sets) 
dd.add1 <- dd.all %>% mutate(y = NA, strain = ' ') %>%
  distinct()
dd.add2 <- dd.all %>% mutate(y = NA, strain = '  ')%>%
  distinct()
dd.all <- rbind(dd.all, dd.add1) %>%
  rbind(., dd.add2)

dd.all <- dd.all %>% 
  filter( strain %in% c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J", ' ', '  ')) %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))
 

cols.all <- c("TAUM"='#F8766D', "TAUL"='#00BA38', "JAX"='#619CFF',
          "Wiltshire2"= 'deepskyblue',
          "Tordoff3"='maroon3',
          "Tarantino2"='deeppink',
          "Crowley2"='darkgoldenrod', 
          "Crabbe4"='seagreen' )

# Body weight boxplot (1)
cols <- c("TAUM"='#F8766D', "TAUL"='#00BA38', "JAX"='#619CFF',
              "Tordoff3"='maroon3',
              "Crowley2"='darkgoldenrod')
(Bodyweight_boxplot <- dd.all %>%
    filter(measure.name%in%c("Body.Weight")) %>%
  mutate(measure.name='Body Weight (gr)',
         lab = factor(lab, levels = c('TAUM', 'TAUL', 'JAX', 'Crowley2', 'Tordoff3'))) %>%
  ggplot() +
  geom_boxplot(aes(x=strain, y=y, fill=lab))+
  facet_grid(measure.name ~ treatment + sex)+
  theme(axis.text.x = element_text(angle=90), 
        legend.position = 'top' ) +
  scale_fill_manual('',values = cols)+
  xlab('Genotype') + ylab(element_blank())) #'Body Weight (grams)'))


# Grip Strength 7 min (S1)
cols <- cols.all[c(1:3,8)]

(grip_boxplot <- dd.all %>%
    filter(measure.name%in%c( "grip.avg" )) %>%
    filter(is.na(y) | (y>0)) %>%
    mutate(y=y^3) %>%
    mutate(measure.name = "Grip Strength (gr)",
           lab = factor(lab, levels = c('TAUM', 'TAUL', 'JAX', 'Crabbe4'))) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))



# Tail suspention, 7  (S4)
cols <- cols.all[c(1:4)]
(TST7_boxplot <- dd.all %>%
    filter(measure.name%in%c( "tst_7min_percent" )) %>%
    mutate(measure.name = 'Tail Suspension Test\n(% of 7 min)',
           lab = factor(lab, levels = c('TAUM','TAUL','JAX', 'Wiltshire2'))) %>%
    mutate(y =  (exp(y)*100.1 - 0.1)/(100*( 1 + exp(y) )) ) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()) )




# Open Field center time, small 10 (S5)
cols <- cols.all[c(1:4)]
(CTsmall10_boxplot <- dd.all %>%
    filter(measure.name%in%c( "OFTsmall_centertime_10m_percent" )) %>%
    mutate(measure.name = "Center Time (%)\n10 Minutes in Small Arena",
           lab = factor(lab, levels = c('TAUM','TAUL','JAX', 'Wiltshire2'))) %>%
    mutate(y =  (exp(y)*100.1 - 0.1)/(100*( 1 + exp(y) )) ) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))



# Open Field center time, large 10 (S3)
cols <- cols.all[c(1:3,6)]
(CTlarge10_boxplot <- dd.all %>%
    filter(measure.name %in% c("OFTlarge_centertime_10m_percent")) %>%
    mutate(measure.name = "Center Time (%)\n10 Minutes in Large Arena",
           lab = factor(lab, levels = c('TAUM','TAUL','JAX', 'Tarantino2'))) %>%
    mutate(y =  (exp(y)*100.1 - 0.1)/(100*( 1 + exp(y) )) ) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))



# Distance Travelled, large 10 (2)
cols <- cols.all[c(1:3,6)]
(DTlarge10_boxplot <- dd.all %>%
    filter(measure.name %in% c("OFTlarge_dist_10m_sec")) %>%
    mutate(measure.name = "Distance Travelled (cm)\n10 Minutes in Large Arena",
           lab = factor(lab, levels = c('TAUM','TAUL','JAX', 'Tarantino2'))) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))


# Distance Travelled, large 10 (2)
cols <- cols.all[c(1:3,6)]
(DTlarge10_boxplot <- dd.all %>%
    filter(measure.name %in% c("OFTlarge_dist_10m_sec")) %>%
    mutate(measure.name = "Distance Travelled (cm)\n10 Minutes in Large Arena",
           lab = factor(lab, levels = c('TAUM','TAUL','JAX', 'Tarantino2'))) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))


# Distance Travelled, small 10 (S2)
cols <- cols.all[c(1,2,4)]
(DTsmall10_boxplot <- dd.all %>%
    filter(measure.name %in% c("OFTsmall_dist_10m_sec")) %>%
    mutate(measure.name = "Distance Travelled (cm)\n10 Minutes in Small Arena",
           lab = factor(lab, levels = c('TAUM','TAUL','Wiltshire2'))) %>%
    ggplot() +
    geom_boxplot(aes(x=strain, y=y, fill=lab))+
    facet_grid(measure.name ~ treatment + sex)+
    theme(axis.text.x = element_text(angle=90), 
          legend.position = 'top' ) +
    scale_fill_manual('',values = cols)+
    xlab('Genotype') + ylab(element_blank()))







library(gridExtra)


png('paper_figures_pg1.png')
# Figure 1:
grid.arrange(Bodyweight_boxplot +
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(angle=90), 
                     legend.position = 'top', 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             bodyweight_intplot + scale_color_discrete('') + 
               ylab(element_blank()) + theme(legend.position = 'none', 
                                             strip.text.x = element_blank(),
                                             strip.background.x = element_blank(),
                                             axis.text.y = element_text(angle=90), 
                                             plot.margin = unit(c(0,0,0,.1),'cm')))
dev.off()
png('paper_figures_pg2.png', width = 11)
# Figure 2:
grid.arrange(DTlarge10_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             DT10_large_intplot + theme(legend.position = 'none', 
                                        strip.text.x = element_blank(),
                                        strip.background.x = element_blank(), 
                                        plot.margin = unit(c(0,0,0,.1),'cm'),
                                        axis.text.y = element_text(angle=90)))
dev.off()
png('paper_figures_pg3.png', width = 11)
# Figure S1:
grid.arrange( grip_boxplot +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.y = element_text(angle=90), 
                      legend.position = 'top', 
                      plot.margin = unit(c(0,0,0,.1),'cm') ) ,
              gripavg_intplot + scale_color_discrete('') + 
                theme(legend.position = 'none', 
                      strip.text.x = element_blank(),
                      axis.text.y = element_text(angle=90),
                      strip.background.x = element_blank(), 
                      plot.margin = unit(c(0,0,0,.1),'cm')))
dev.off()
png('paper_figures_pg4.png', width = 11)
# Figure S2:
grid.arrange( DTsmall10_boxplot + 
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(), 
                      legend.position = 'top',
                      axis.text.y = element_text(angle=90), 
                      plot.margin = unit(c(0,0,0,.1),'cm') ) ,
              DT10_small_intplot + scale_color_manual('',values = cols) +
                theme(legend.position = 'none', 
                      strip.text.x = element_blank(),
                      strip.background.x = element_blank(), 
                      plot.margin = unit(c(0,0,0,.1),'cm'),
                      axis.text.y = element_text(angle=90)))
dev.off()
png('paper_figures_pg5.png', width = 11)
# Figure S3:
grid.arrange(CTlarge10_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             CT10_large_intplot + theme(legend.position = 'none', 
                                        strip.text.x = element_blank(),
                                        strip.background.x = element_blank(), 
                                        plot.margin = unit(c(0,0,0,.1),'cm'),
                                        axis.text.y = element_text(angle=90)))

dev.off()
png('paper_figures_pg6.png', width = 11)
# Figure S4:
grid.arrange(TST7_boxplot+
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             tst7_intplot +  theme(legend.position = 'none', 
                                   strip.text.x = element_blank(),
                                   strip.background.x = element_blank(), 
                                   plot.margin = unit(c(0,0,0,.1),'cm'),
                                   axis.text.y = element_text(angle=90)))
dev.off()
png('paper_figures_pg7.png', width = 11)
# Figure S5:
grid.arrange(CTsmall10_boxplot +
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             CT10_small_intplot + scale_color_manual('',values = cols) + 
               theme(legend.position = 'none', 
                     strip.text.x = element_blank(),
                     strip.background.x = element_blank(), 
                     plot.margin = unit(c(0,0,0,.1),'cm'),
                     axis.text.y = element_text(angle=90)))

dev.off()
