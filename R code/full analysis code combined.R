my.path <- '/Users/Iman/Documents/GitHub/Improving-replicability-using-interaction-2021/'

source(paste0(my.path, 'R code/Analysis_Functions.R'))

compute.df.sater <- function(s2int,df.int,s2pooled,n1,n2){
  inv.v <- (((1/n1+1/n2)*s2pooled)^2*1/(n1+n2-2)+(2*s2int)^2*1/df.int)/
    (((1/n1+1/n2)*s2pooled+2*s2int)^2)
  df <- 1/inv.v
  return(df)
}

compute.df.sater.2pairs <- function(s2int, df.int, s2pooled, n1, n2, n3, n4){
  inv.v <- (((1/n1+1/n2+1/n3+1/n4)*s2pooled)^2*1/(n1+n2+n3+n4-4)+(4*s2int)^2*1/df.int)/
    (((1/n1+1/n2+1/n3+1/n4)*s2pooled+4*s2int)^2)
  df <- 1/inv.v
  return(df)
}


alpha.threshold = 0.05
alpha.threshold2 = 0.05


#################################################
################## Data Import ##################
#################################################


### import gxl estimates from IMPC
gxl.info <- read.csv(paste0(my.path, 'gxl_estimates_from_IMPC.csv'),header = T , stringsAsFactors = F)



gxl2.fac <- gxl.info$gamma2.estimator
gxl.L <- gxl.info$labs.int
gxl.S <- gxl.info$strains.int
gxl.transformations <- gxl.info$transformation 
names(gxl2.fac) <-  names(gxl.S) <- names(gxl.L) <-
  names(gxl.transformations) <-  gxl.info$File.name


##### TAU-JAX import after manipulation
##############

OFT_combined_data <- read_csv(file = paste0(my.path, 'OFT_combined_data.csv'))


TST_combined_data <- read_csv(file = paste0(my.path, 'TST_combined_data.csv')) %>% select(-c('...1'))
bw_combined_data <- read_csv(file = paste0(my.path, 'bw_combined_data.csv')) %>% 
  rename('Body.Weight'= 'Body Weight')
grip_combined_data <- read_csv(file = paste0(my.path, 'grip_combined_data.csv'))

t.add <-  OFT_combined_data %>% 
  filter(  name %in% c('OFTlarge_centertime_20m_percent', 'OFTlarge_centertime_10m_percent', 'OFTsmall_centertime_10m_percent', 'OFTsmall_centertime_20m_percent'),
           value > 0 )  %>%
  group_by(name) %>% 
  summarise(t = 0.5*min(value))

OFT_combined_data1  <- OFT_combined_data %>% 
  filter(  name %in% c('OFTlarge_centertime_20m_percent', 'OFTlarge_centertime_10m_percent', 'OFTsmall_centertime_10m_percent', 'OFTsmall_centertime_20m_percent') )  %>% 
  left_join(t.add)
OFT_combined_data1$value <- log( (OFT_combined_data1$value + OFT_combined_data1$t )/
                                   ( 1 + OFT_combined_data1$t - OFT_combined_data1$value )) 

OFT_combined_data <- rbind(OFT_combined_data1 %>% select(-c('t')),
                           OFT_combined_data %>% filter( ! name %in% c('OFTlarge_centertime_20m_percent', 'OFTlarge_centertime_10m_percent', 'OFTsmall_centertime_10m_percent', 'OFTsmall_centertime_20m_percent') ))
# OFT_combined_data[1447,'value'] <- 100

grip_combined_data <- grip_combined_data %>% 
  rename('value'='grip.avg') %>%
  mutate(name = 'grip.avg', value = value^(1/3))

t6 <- TST_combined_data$tst_6min_percent
t6 <- t6[t6>0] ; t6 <- 0.5*min(t6, na.rm = T)
t7 <- TST_combined_data$tst_7min_percent
t7 <- t7[t7>0] ; t7 <- 0.5*min(t7, na.rm = T)
TST_combined_data <- TST_combined_data %>% 
  mutate( tst_6min_percent = log((tst_6min_percent+t6)/(1+t6-tst_6min_percent)),
          tst_7min_percent = log((tst_7min_percent+t7)/(1+t7-tst_7min_percent))) %>%
  pivot_longer(cols = c('tst_6min_percent', 'tst_7min_percent'), names_to = 'name', values_to = 'value')




data_sets_list <- c('OFT_combined_data', 'TST_combined_data', 'bw_combined_data', 'grip_combined_data')
all_data_sets_list <- list()
for( ds in data_sets_list ){
  tds <- eval(parse(text = ds ))
  tds <- tds[, (!names(tds) %in% c('X','X1','...1'))] %>%
    mutate(strain = as.character(strain) , lab = as.character(lab), sex = as.character(sex))
  
  if(  all(c('value', 'name') %in% names(tds) )){
    tds <- tds %>% rename('y.name' = 'name', 'y'='value')
  }else{
    tds <- tds %>% pivot_longer(cols = -c('lab', 'treatment', 'sex', 'strain'), 
                                names_to = 'y.name', values_to = 'y' )  
  }
  tds <- tds %>%  unite("y.name", c(y.name, treatment), sep = ':', remove = T) %>% group_split(y.name,sex)
  tds.names <- lapply(tds, function(x) paste(x[1,'y.name'],x[1,'sex'], sep = ':') ) %>%
    unlist(.)
  tds <- lapply(tds, function(x) {x %>% mutate( y.name = paste(y.name,sex, sep = ':') )} )
  
  tds.keep <- lapply(tds, FUN = function(x){ 
    x <- x %>% filter( !is.na(y), !is.nan(y) ) %>%
      group_by(lab,strain) %>% summarise(n=n()) 
    (nrow(x)>1)&all(x$n>2)} )  %>%
    unlist(.) %>% which(.)
  tds <- lapply(tds[tds.keep], FUN = function(x) clean_dat(x, desig.variable = 'y') )
  
  
  all_data_sets_list <- c(tds, all_data_sets_list) 
}

all_data_sets <- do.call(what = rbind, all_data_sets_list) %>%
  separate(col = 'y.name', into = c('measure.name' , 'treatment', 'sex'), sep = ':') %>%
  filter( ! (strain %in% c('C3H/HeJ'))) %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))

# all_data_gxlEstimates <- read_csv(my.path, 'all_data_gxlEstimates.csv') %>%
#   pivot_longer(cols = c(error_sd_value,  int_sd_value), names_to = 'value_type',values_to = 'value') %>%
#   mutate(x = recode(value_type,
#                     'int_sd_value'= ' '     , 'error_sd_value'='  ')) %>%
#   mutate(value_type = recode(value_type,"int_sd_value"='Interaction SD',  "error_sd_value"='Error SD'))



##############################################
################## Analyses ##################
##############################################


##################
# Analysis of TAU-JAX multi lab data  (Two way anova models)
all_data_sets_csvs <- all_data_sets %>% mutate(y = round(y, digits = 3)) %>%
  split(., list(all_data_sets$measure.name,all_data_sets$sex, all_data_sets$treatment)) 
all_data_sets_csvs <- Filter(function(x) dim(x)[1] > 0, all_data_sets_csvs)

entry.place <- 0
all_data_analyses_list <- all_data_altimumtruth_list <- dfs_TAUJAX <-
  all_data_means <- all_data_gxlEstimates <- all_data_errorEstimates <- 
  all_data_intEstimates <- list()
for( i in 1:length(all_data_sets_list)){
  temp_dat <- all_data_sets_list[[i]] %>%
    filter(! strain %in% c('C3H/HeJ')) #, 'CBA/J'))
  temp_dat <- droplevels(temp_dat)
  
  temp <- multi_lab_analysis(tidy_data = temp_dat, report.full.summary = T )
  
  entry.place <- entry.place + 1  
  entry.name <- ((all_data_sets_list)[[i]])$y.name[1]
  print(entry.name)
  print(i)
  
  all_data_analyses_list[[  entry.place ]] <- temp[["Post-Hoc Results"]]
  # all_data_altimumtruth_list[[ entry.place ]] <- temp[["Altimum Truth.Unadjusted"]]
  all_data_means[[ entry.place ]] <- temp[['summary_table']]
  all_data_gxlEstimates[[ entry.place ]] <- temp$gxl_factor
  all_data_errorEstimates[[ entry.place ]] <- temp$error_sd
  all_data_intEstimates[[ entry.place ]] <- temp$int_sd
  
  names(all_data_analyses_list)[entry.place] <-
    # names(all_data_altimumtruth_list)[entry.place] <- 
    names(all_data_means)[entry.place] <- 
    names(all_data_gxlEstimates)[entry.place] <- 
    names(all_data_errorEstimates)[entry.place] <- 
    names(all_data_intEstimates)[entry.place] <- 
    entry.name
  
  strain_TAUJAX <- all_data_sets_list[[i]] %>% 
    group_by(y.name,strain,sex) %>% summarise( n = 1) %>% 
    group_by(y.name,sex) %>% summarise( nstrain_TAUJAX = sum(n)-1 )
  
  lab_TAUJAX <- all_data_sets_list[[i]] %>% 
    group_by(y.name,lab,sex) %>% summarise( n = 1) %>% 
    group_by(y.name,sex) %>% summarise( nlab_TAUJAX = sum(n)-1 )
  
  df2_TAUJAX <- all_data_sets_list[[i]] %>% 
    group_by(y.name,strain, lab,sex) %>% summarise( n=n()) %>% 
    group_by(y.name,sex) %>% summarise( df2_TAUJAX = sum(n-1) )
  
  df1_TAUJAX <- all_data_sets_list[[i]] %>% 
    dplyr::select( y.name,strain, lab ,sex  ) %>% distinct() %>%
    group_by(y.name,sex) %>% summarise( df1_TAUJAX = n(), nlab = n_distinct(lab), nstrain = n_distinct(strain) ) %>%
    mutate( df1_TAUJAX = df1_TAUJAX - nlab - nstrain + 1 )
  
  lab_TAUJAX <- lab_TAUJAX %>% mutate(nlab_TAUJAX = nlab_TAUJAX + 1)
  strain_TAUJAX <- strain_TAUJAX %>% mutate(nlab_TAUJAX = nstrain_TAUJAX + 1)
  
  dfs_TAUJAX[[i]] <- merge(df1_TAUJAX, df2_TAUJAX) %>%
    left_join(strain_TAUJAX) %>% left_join(lab_TAUJAX) %>%
    mutate(gxl2.fac = temp$gxl_factor^2 ) %>%
    separate(col = 'y.name', into = c('measure.name.MPDname', 'treatment.MPDname'),sep = ':')
}

dfs_TAUJAX <- do.call(rbind, dfs_TAUJAX)
write_csv(dfs_TAUJAX %>% mutate(gxl2.fac = round(gxl2.fac, digits = 3)), 'TAUJAX_degrees_of_freedom.csv')

all_data_means <-  do.call(rbind, all_data_means)
rownames(all_data_means) <- NULL
rownames(dfs_TAUJAX) <- NULL

for (i in 1:length(all_data_analyses_list)){
  all_data_analyses_list[[i]] <- all_data_analyses_list[[i]] %>% mutate(measure.name = names(all_data_analyses_list)[i]) 
}
all_data_analyses_list <- do.call(rbind, all_data_analyses_list) %>%
  separate(measure.name, into = c('measure.name', 'treatment', 'sex'),sep = ':')


all_data_gxlEstimates <- do.call(cbind,all_data_gxlEstimates)
all_data_errorEstimates <- do.call(cbind,all_data_errorEstimates)

all_data_gxlEstimates <- data.frame(gxl_fac_value = as.numeric(all_data_gxlEstimates), 
                                    error_sd_value = as.numeric(all_data_errorEstimates), 
                                    int_sd_value = as.numeric(all_data_intEstimates), 
                                    measure = colnames(all_data_gxlEstimates)) %>%
  mutate( `TAU-JAX gxl factor`=round(int_sd_value/error_sd_value, digits = 4) ) %>%
  separate(col = 'measure', into = c('measure', 'treatment', 'sex'), sep = ':') %>%
  pivot_longer(cols = c(error_sd_value,  int_sd_value), names_to = 'value_type',values_to = 'value') %>%
  mutate(x = recode(value_type,
                    'int_sd_value'= ' '     , 'error_sd_value'='  ')) %>%
  mutate(value_type = recode(value_type,"int_sd_value"='  Interaction SD',  "error_sd_value"='   Error SD'))

write_csv(all_data_gxlEstimates %>%
            dplyr::select(c("measure", "treatment", "sex",'TAU-JAX gxl factor', 'value_type', 'value' )) %>%
            mutate(value = round(value, digits = 4)) %>%
            pivot_wider(names_from = 'value_type', values_from = 'value'),
          path = 'all_data_gxlEstimates.csv')

write_csv( dfs_TAUJAX  %>%
             mutate( gxl2.fac = round(gxl2.fac, digits = 5) ),
           path = 'TAUJAX_DFs_gxl.csv')
write_csv(all_data_analyses_list    , path = 'TAUJAX_analyses_list.csv')
write_csv(all_data_means %>% arrange(lab, sex, measure.name, treatment ), 'TAUJAX_means.csv')


gxl_IMPC_all <- read_csv(paste0(my.path, 'IMPC_gxl_estimates.csv')) %>%
  dplyr::select(-c('...1')) %>% rename('IMPC_code' = 'X1')

##################
## Three way ANOVA for treatment X lab X strain interaction
## Perform post-hoc to compare the treatment-effect across strains, i.e., Strain1(Flox-cont) Vs. Strain2(Flox-cont)

all_data_sets_separated <- all_data_sets %>% 
  split(., list(all_data_sets$measure.name,all_data_sets$sex )) 
all_data_sets_separated <- Filter(function(x) dim(x)[1] > 0, all_data_sets_separated)

keep_measures <- lapply(X = all_data_sets_separated, FUN = function(x) n_distinct(x$treatment)>1) %>% unlist(.)
all_data_sets_separated <- all_data_sets_separated[ keep_measures  ]

ctrl <- lmeControl(opt='optim')
threewya.aov.resuls <- tibble(measure.name = names(all_data_sets_separated),
                              sex = NA, 
                              nlab = NA, nstrain = NA, 
                              trt.str.lab.g = NA, trt.str.lab.s2 = NA, str.lab.g = NA, 
                              trt.str.lab.pv = NA, trt.str.pv = NA )

for( t in seq_along(names(all_data_sets_separated))){
  td <- names(all_data_sets_separated)[t]
  tds <- all_data_sets_separated[[t]]
  fit.lme.treat.interaction <- lme(y ~ strain*treatment, random=~1|lab/strain/treatment,
                                   data=tds, control=ctrl,
                                   method="REML",na.action='na.omit')
  fit.lme.MINUS.treat.fixed.interaction <- lme(y ~ strain+treatment, random=~1|lab/strain/treatment,
                                               data=tds, control=ctrl,
                                               method="REML",na.action='na.omit')
  
  (fit.lme.MINUS.treat.interaction <- lme(y ~ strain*treatment, random=~1|lab/strain,
                                          data=tds, control=ctrl,
                                          method="REML",na.action='na.omit'))
  
  threewya.aov.resuls$trt.str.lab.pv[t] <- as.data.frame(anova(fit.lme.treat.interaction, fit.lme.MINUS.treat.interaction))[2,'p-value']
  threewya.aov.resuls$trt.str.pv[t] <- as.data.frame(anova(fit.lme.treat.interaction, fit.lme.MINUS.treat.fixed.interaction))[2,'p-value']
  threewya.aov.resuls$str.lab.pv[t] <- as.data.frame(anova(fit.lme.treat.interaction, fit.lme.MINUS.treat.interaction))[2,'p-value']
  tem <- VarCorr(fit.lme.treat.interaction)
  tem1 <- as.matrix(tem)[6,2]
  tem2 <- as.matrix(tem)[7,2]
  tem3 <- as.matrix(tem)[4,2]
  
  threewya.aov.resuls$trt.str.lab.g[t] <- as.numeric(tem1)/as.numeric(tem2)
  threewya.aov.resuls$trt.str.lab.s2[t] <- as.numeric( as.matrix(tem)[6,2] )
  threewya.aov.resuls$str.lab.g[t] <- as.numeric(tem3)/as.numeric(tem2)
  threewya.aov.resuls$sex[t] <- tds$sex[1]
  threewya.aov.resuls$nlab[t] <- n_distinct(tds$lab)
  threewya.aov.resuls$nstrain[t] <- n_distinct(tds$strain)
}

threewya.aov.resuls <- threewya.aov.resuls %>% mutate( `gamma of lab-strain` = round(unname(str.lab.g), digits = 3),
                                                       `gamma of lab-strain-treatment` = round(unname(trt.str.lab.g), digits = 3),
                                                       `LRT pvalue of lab-strain-treatment` = round(unname(trt.str.lab.pv), digits = 3),
                                                       `LRT pvalue of strain-treatment` = trt.str.pv ) %>%
  arrange(`LRT pvalue of lab-strain-treatment`) %>%
  separate(measure.name, into = c('measure.name', 'discard'), sep = '.Male' ) %>%
  dplyr::select(-c('discard'))

write_csv(threewya.aov.resuls, file = 'results_of_three_way_anova_at_TAUJAX.csv',col_names = T)


## Begin: defining contrasts: 
v1 <- rep(0, 12); names(v1) <- names((fit.lme.treat.interaction$coefficients)$fixed) 
names(v1)[1] <- 'strainBALB/cJ'
places <- c(7:12); places <- combn(places, m = 2)
names(v1)[c(7:12)]
v.all <- c()
for (i in 1:ncol(places) ) {
  v2 <- rep(0, 12)
  v2[places[,i]] <- c(1,-1)
  v3 <- paste(v2, collapse = ', ')
  v3 <- c( ' = c(', v3, ')' )
  v3 <- paste(v3, collapse = ' ')
  nm <- paste(names(v1)[places[,i]], collapse = ' - ')
  nm <- c('`', nm, '`') ; nm <- paste(nm, collapse = '')
  v3 <- paste(c(nm, v3), collapse = '')
  v.all[i] <- v3 
}
paste(v.all, collapse = ', ')

c1 <- list( # text for contrasts comparing effects between strains
  `treatmentFluoxetine - strainBTBR:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0 ),
  `treatmentFluoxetine - strainC57BL/6J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0 ),
  `treatmentFluoxetine - strainDBA/2J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0 ), 
  `treatmentFluoxetine - strainSWR/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0 ),
  `treatmentFluoxetine - strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1 ),
  `strainBTBR:treatmentFluoxetine - strainC57BL/6J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0 ),
  `strainBTBR:treatmentFluoxetine - strainDBA/2J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0 ), 
  `strainBTBR:treatmentFluoxetine - strainSWR/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0 ), 
  `strainBTBR:treatmentFluoxetine - strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1 ),
  `strainC57BL/6J:treatmentFluoxetine - strainDBA/2J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0 ),
  `strainC57BL/6J:treatmentFluoxetine - strainSWR/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0 ),
  `strainC57BL/6J:treatmentFluoxetine - strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 ), 
  `strainDBA/2J:treatmentFluoxetine - strainSWR/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0 ),
  `strainDBA/2J:treatmentFluoxetine - strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 ), 
  `strainSWR/J:treatmentFluoxetine - strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1 )
)
c2 <- list( # text for treatment effect per strain contrasts
  `treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ),
  `strainBTBR:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ),
  `strainC57BL/6J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ),
  `strainDBA/2J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ),
  `strainSWR/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ),
  `strainCBA/J:treatmentFluoxetine` = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ))


# Begin: post hoc of Floxitine effect between strains
all.contrasts <- all.contrasts2 <- list()
for( t in seq_along(names(all_data_sets_separated))){
  td <- names(all_data_sets_separated)[t]
  tds <- all_data_sets_separated[[t]]
  fit.lme.treat.interaction <- lme(y ~ strain*treatment, random=~1|lab/strain/treatment,
                                   data=tds, control=ctrl,
                                   method="REML",na.action='na.omit')
  t.contrasts <-  tibble( contrast.name = names(c1), 
                          measure.name = names(all_data_sets_separated)[t], 
                          `F-value` = NA, `REML post-hoc` = NA )
  for( ci in 1:length(c1)){ # treatment effect between strains
    contrast.i <- anova(fit.lme.treat.interaction, L = c1[[ci]] )
    t.contrasts[ci, 'REML post-hoc'] <- as.data.frame(contrast.i)[1,4]
    t.contrasts[ci, 'F-value'] <- as.data.frame(contrast.i)[1,3]
  }
  all.contrasts[[t]] <- t.contrasts
  
  t.contrasts2 <-  tibble( contrast.name = names(c2), 
                           measure.name = names(all_data_sets_separated)[t], 
                           `F-value` = NA, `REML post-hoc` = NA )
  for( ci in 1:length(c2)){ # treatment effect per strain
    contrast.i <- anova(fit.lme.treat.interaction, L = c2[[ci]] )
    t.contrasts2[ci, 'REML post-hoc'] <- as.data.frame(contrast.i)[1,4]
    t.contrasts2[ci, 'F-value'] <- as.data.frame(contrast.i)[1,3]
  }
  all.contrasts2[[t]] <- t.contrasts2
  
}
all.contrasts <- do.call( rbind, all.contrasts )
all.contrasts2 <- do.call( rbind, all.contrasts2 )

all.contrasts <- all.contrasts %>% 
  separate( contrast.name, into = c('strain1', 'strain2'), sep = ' - '  ) %>%
  separate( strain2, into = c('strain2', 'delete'), sep = ':treatmentFluoxetine'  ) %>%
  separate( strain2, into = c('delete', 'strain2'), sep = 'strain'  ) %>%
  separate( strain1, into = c('delete', 'strain1'), sep = 'strain'  ) %>%
  separate( strain1, into = c('strain1', 'delete'), sep = ':treatmentFluoxetine'  ) %>%
  dplyr::select(-c('delete')) %>%
  mutate( strain1 = replace(strain1, is.na(strain1), 'BALB/cJ' ),
          strain2 = replace(strain2, is.na(strain2), 'BALB/cJ' ))  %>%
  separate( measure.name, into = c('measure.name', 'sex'), sep = '.Male'  ) %>%
  mutate( sex = 'Male' ) %>%
  dplyr::select( c(measure.name, strain1, strain2, sex, `REML post-hoc` ) )

reord <-  which(all.contrasts$strain1 < all.contrasts$strain2 )
all.contrasts[reord, c('strain1', 'strain2')] <- all.contrasts[reord, c('strain2', 'strain1')]


all.contrasts2 <- all.contrasts2 %>%
  separate( contrast.name, into = c('strain', 'treatment2'), sep = ':treatment'  ) %>%
  mutate(treatment1 = 'Control')%>%
  mutate( treatment2 = replace(treatment2, is.na(treatment2), 'Fluoxetine' ))  %>%
  mutate( strain = recode(strain, 'treatmentFluoxetine'= 'BALB/cJ' ))  %>%
  separate( measure.name, into = c('measure.name', 'sex'), sep = '.Male'  ) %>%
  mutate( sex = 'Male' ) %>%
  dplyr::select( c(measure.name, strain, treatment1, treatment2, sex, `REML post-hoc` ) ) %>%
  separate( strain, into = c('delete','strain'), sep = 'strain'  ) %>%
  mutate(strain = recode(strain, .missing = ''))

all.contrasts2$strain <- paste0(all.contrasts2$strain, all.contrasts2$delete) 
all.contrasts2 <- all.contrasts2 %>% select(-c('delete'))


write_csv(all.contrasts, file = 'Floxitine effect: TAUJAX_theewayANOVA_posthoc.csv')
write_csv(all.contrasts2, file = 'Floxitine effect per strain: TAUJAX_theewayANOVA_posthoc.csv')
### End: post hoc of Floxitine effect between strains

### begin : post hoc of Floxitine effect with gxl


### End: post hoc of Floxitine effect with gxl

### Begin: Tests of three-way interaction effect
gxl_TAUJAX <- read_csv( paste0(my.path, 'TAUJAX_DFs_gxl.csv')) %>% 
  rename( 'measure.name' = 'measure.name.MPDname', 'treatment' = 'treatment.MPDname', "gxl2.TAUJAX"="gxl2.fac") 
gxl_TAUJAX <- gxl_TAUJAX %>% arrange(measure.name,treatment, sex ) %>%
  # dplyr::select(measure.name,treatment,sex,gxl2.TAUJAX) %>% 
  mutate(gxl.TAUJAX=sqrt(gxl2.TAUJAX))
gxl_TAUJAX <- gxl_TAUJAX %>%
  dplyr::select(measure.name,treatment,sex,gxl.TAUJAX)
gxl_IMPC <- read_csv( paste0(my.path, 'gxl_estimates_from_IMPC.csv')) %>%
  rename("gxl2.IMPC" = "gamma2 estimator", 'measure.name'='MPD_code') %>%  dplyr::select(`File name`, measure.name, sex, IMPC_code, gxl2.IMPC)
gxl_IMPC <- gxl_IMPC %>% left_join(gxl_IMPC_all) %>% mutate(Sd_interaction = sqrt(S2_interaction))
gxl_IMPC <- gxl_IMPC %>% 
  rename('File.name'='File name') %>%
  mutate( `gamma of lab-strain (by IMPC)` = round(sqrt(gxl2.IMPC), digits=3) ) %>%
  dplyr::select(measure.name, sex, `gamma of lab-strain (by IMPC)`) %>%
  unique(.)


temp <- gxl_TAUJAX %>%  
  full_join(threewya.aov.resuls) %>% 
  filter(!is.na(`LRT pvalue of lab-strain-treatment`)) %>%
  mutate(`gamma of lab-strain (control only)` = round(gxl.TAUJAX, digits = 3)) %>%
  left_join(gxl_IMPC) %>% 
  dplyr::select(-c('gxl.TAUJAX')) %>%
  dplyr::select( "measure.name", "treatment","sex", 'gamma of lab-strain (by IMPC)',
                 "gamma of lab-strain (control only)", "gamma of lab-strain",
                 "gamma of lab-strain-treatment", "LRT pvalue of lab-strain-treatment") %>%
  filter(treatment == 'Control') 

write_csv(temp, file = 'three_way_interaction_test.csv')
### End: Tests of three-way interaction effect


##################################################################################################################
## Begin: Gxl adjustment for pairwise comparisns of treatment effect ( Floxitine only, which is only in Wiltshire)
##################################################################################################################

threewya.aov.resuls.TAUJAX <- read_csv('results_of_three_way_anova_at_TAUJAX.csv') %>% 
  dplyr::select(measure.name, sex, nlab, nstrain, trt.str.lab.g, trt.str.lab.s2, `LRT pvalue of lab-strain-treatment`, `LRT pvalue of strain-treatment`)

file.names <- c("Lab_Wiltshire2_DistanceTraveled_Males.csv",
                "Lab_Wiltshire2_PercentCenter_Males.csv",
                "Lab_Wiltshire2_TST_Males.csv" )

pairs.all <-  list()
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
                             'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J')) %>%
    filter(strain %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') ) # "C3H/HeJ" 
  
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
    mutate(treatment = recode(treatment, 'fluoxetine'='Fluoxetine', 'control'='Control'))
  
  dd.means <- dd %>%
    group_by(strain, lab, sex, y.name, treatment) %>% 
    summarise(mean = mean(y,na.rm = T) , sd = sd(y,na.rm = T) , n = n() )
  
  
  pairs.right <- expand_grid( strain1 = as.character(unique( dd$strain)), 
                              y.name = unique(dd$y.name), 
                              treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
    filter( treatment1 < treatment2 , !is.na(strain1) )
  
  pairs.right <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment1','y.name', 'mean1R', 'sd1R', 'n1R')) %>%
    right_join(pairs.right)
  
  pairs.right <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment2', 'y.name', 'mean2R', 'sd2R', 'n2R')) %>%
    right_join(pairs.right)
  
  
  pairs.left <- expand_grid( strain2 = as.character(unique( dd$strain)), 
                             y.name = unique(dd$y.name), 
                             treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
    filter( treatment1 < treatment2 , !is.na(strain2) )
  
  pairs.left <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment1', 'y.name', 'mean1L', 'sd1L', 'n1L')) %>%
    right_join(pairs.left)
  
  pairs.left <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment2', 'y.name','mean2L', 'sd2L', 'n2L')) %>%
    right_join(pairs.left)
  
  pairs.right <- pairs.right %>% 
    mutate(diff.R = ( mean1R - mean2R ) )
  
  pairs.left <- pairs.left %>% 
    mutate(diff.L = ( mean1L - mean2L ) )
  
  pairs.all[[file.name]] <- full_join(pairs.right, pairs.left) %>%
    filter(strain1 > strain2) %>% 
    rename('measure.name'='y.name' )
}
pairs.all <- do.call(rbind, pairs.all)  %>%
  mutate(sex = recode(sex, 'Males' = 'Male', 'Females' = 'Female')) %>%
  left_join(threewya.aov.resuls.TAUJAX) %>%
  mutate(diff = (diff.L - diff.R),
         s2.pooled = ( (n1L-1)*sd1L^2 + (n1R-1)*sd1R^2 + (n2L-1)*sd2L^2 + (n2R-1)*sd2R^2 ) / ( (n1L-1) + (n1R-1) + (n2L-1) + (n2R-1) ),
         n.op = ( (n1L^-1) + (n1R^-1) + (n2L^-1) + (n2R^-1) ) ) %>%
  mutate( diff.se2 = s2.pooled * (n.op),
          `diff.se2.gxladj (by TAUJAX)` = s2.pooled * (n.op + 4*trt.str.lab.g^2 )) %>%
  mutate( `t-test` = diff / sqrt(diff.se2),
          `t-test gxl ajd (by TAUJAX)` = diff / sqrt(`diff.se2.gxladj (by TAUJAX)`),
          `df.satt (by TAUJAX)` = compute.df.sater.2pairs(s2int = trt.str.lab.s2, df.int = (nstrain-1)*(nlab-1)*(2-1),
                                            s2pooled = s2.pooled, n1 = n1L, n2 = n1R, n3 = n2L, n4 = n2R ) ) %>%
  mutate( `p-value.common` = pt( abs(`t-test`), df = n1L + n1R + n2L + n2R - 4, lower.tail = F  )*2,
          `p-value gxl ajd (by TAUJAX)` = pt( abs(`t-test gxl ajd (by TAUJAX)`), df = `df.satt (by TAUJAX)`, lower.tail = F  )*2)






###################### End: Gxl adjustment for pairwise comparisns of treatment effect



##################
# Analysis of MPD: MPD data import and analysis

dd.means.filenames <- pairs.treatment.all <- list()
all.pairs.comps <- merging.info <- NULL
wd <- paste0(my.path, 'MPD pairs/MPD data - arranged IJ/')
(file.names <- list.files(path = wd, pattern = 'Lab_*'))
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
                             'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J')) %>%
    filter(strain %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') ) # "C3H/HeJ" 
  
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
  
  # Begin: treatment effect for each sole strain
  dd.treatment.effect <- dd %>% 
    separate(col = y.name, into = c('measure.name', 'treatment'), sep = ':' ) %>%
    mutate(treatment = recode(treatment, 'fluoxetine'='Fluoxetine', 'control'='Control'))%>%
    group_by(strain, lab, sex, measure.name, treatment) %>%
    summarise(y.mean = mean(y, na.rm=T), sd = sd(y, na.rm = T), n=n()) 
  pairs.treatment <- dd.treatment.effect %>%
    rename_at(vars(c('strain', 'measure.name', 'treatment','y.mean', 'sd', 'n')), ~ c('strain', 'measure.name', 'treatment2','mean2', 'sd2', 'n2')) 
  pairs.treatment <- dd.treatment.effect %>%
    rename_at(vars(c('strain', 'measure.name', 'treatment','y.mean', 'sd', 'n')), ~ c('strain', 'measure.name', 'treatment1','mean1', 'sd1', 'n1')) %>%
    right_join(pairs.treatment) %>%
    filter(treatment1 < treatment2)
  

  pairs.treatment <-  pairs.treatment %>%
    left_join(threewya.aov.resuls) %>% # start fixing here
    mutate( diff = (mean1 - mean2 )) %>%
    mutate(s2.pooled = ( (n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2- 2)) %>%
    mutate( SE = sqrt(  s2.pooled*(1/n1 + 1/n2 )  ),
            SE.TAUJAX = sqrt(  s2.pooled*(1/n1 + 1/n2 ) + 2*trt.str.lab.s2 )) %>%
    mutate(satter.df.TAUJAX = compute.df.sater(s2int = trt.str.lab.s2,
                                               df.int = (nlab-1)*(nstrain-1),
                                               s2pooled = s2.pooled,  n1 = n1 ,  n2 = n2) ) %>%
    mutate(gxl.adjusted.pv.TAUJAX = 2*pt( abs(diff)/SE.TAUJAX  , df = satter.df.TAUJAX, lower.tail = F) ,
           common.pv       = 2*pt( abs(diff)/SE ,          df = (n1 + n2 - 2) , lower.tail = F) ,
           FDR.common.pv = NA, FDR.gxl.adjusted.pv = NA)
  
  pairs.treatment.all[[file.name]] <- pairs.treatment
  # End: treatment effect for each sole strain
  
  dd.means.filenames[[file.name]] <-  dd %>% group_by(lab, sex, y.name) %>% 
    summarise(mean = mean(y,na.rm = T) , sd = sd(y,na.rm = T) , n = n() )
  dd.means <- dd %>%
    group_by(strain, lab, sex, y.name) %>% 
    summarise(mean = mean(y,na.rm = T) , sd = sd(y,na.rm = T) , n = n() )
  n.strains <- dd %>% summarise(n_distinct(strain) ) %>% unname(.)  
  
  
  pairs <- expand_grid( strain1 = as.character(unique( dd$strain)), strain2 = as.character(unique( dd$strain)),
                        y.name1 = unique(dd$y.name), y.name2 = unique(dd$y.name) ) %>%
    filter( (strain1 > strain2)&(y.name1==y.name2),
            !is.na(strain1) , !is.na(strain2) )
  pairs <- dd.means %>%
    rename_at(vars(c('strain', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'y.name1','mean1', 'sd1', 'n1')) %>%
    right_join(pairs)
  pairs <- dd.means %>%
    rename_at(vars(c('strain', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'y.name2','mean2', 'sd2', 'n2')) %>%
    right_join(pairs)
  pairs <- pairs %>%
    dplyr::select(strain1, strain2, sex, lab, y.name1 , y.name2, mean1, mean2, sd1, sd2, n1, n2)
  pairs <- pairs %>%
    mutate( diff = (mean1 - mean2 )) %>%
    mutate(s2.pooled = ( (n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2- 2),
           gxl2.IMPC = as.numeric(gxl2.fac[file.name]),
           gxl.L.IMPC = gxl.L[file.name],
           gxl.S.IMPC = gxl.S[file.name]) %>%
    mutate( SE = sqrt(  s2.pooled*(1/n1 + 1/n2 )  ),
            SE.IMPC = sqrt(  s2.pooled*(1/n1 + 1/n2 + 2*gxl2.IMPC )  ),
            SE.TAUJAX = NA ) %>%
    mutate(satter.df.IMPC = compute.df.sater(s2int = s2.pooled*gxl2.IMPC,
                                             df.int = (gxl.L.IMPC-1)*(gxl.S.IMPC-1),
                                             s2pooled = s2.pooled,  n1 = n1 ,  n2 = n2),
           satter.df.TAUJAX = NA ) %>%
    mutate(gxl.adjusted.pv.IMPC = 2*pt( abs(diff)/SE.IMPC  , df = satter.df.IMPC ,     lower.tail = F) ,
           gxl.adjusted.pv.TAUJAX = NA, 
           common.pv       = 2*pt( abs(diff)/SE  ,          df = (n1 + n2 - 2) , lower.tail = F) ,
           FDR.common.pv = NA, FDR.gxl.adjusted.pv = NA)
  
  
  pairs$FDR.common.pv <- p.adjust( pairs$common.pv , method = 'BH' )
  pairs$FDR.gxl.adjusted.pv.IMPC <- pairs$FDR.gxl.adjusted.pv.TAUJAX <- NA # p.adjust( pairs$gxl.adjusted.pv , method = 'BH' )
  
  temp.pairs <-  pairs %>% mutate( File.name = file.name) %>%
    dplyr::select(File.name, strain1, strain2, sex,  lab,  y.name1, y.name2, mean1, mean2, sd1, sd2, n1, n2, diff, s2.pooled,
                  SE, SE.IMPC, satter.df.IMPC, satter.df.TAUJAX, 
                  gxl.adjusted.pv.IMPC, gxl.adjusted.pv.TAUJAX,
                  common.pv, FDR.common.pv, FDR.gxl.adjusted.pv ) %>%
    separate(col = 'y.name1', into = c('measure.name', 'treatment1'),sep = ':') %>%
    separate(col = 'y.name2', into = c('measure.name0', 'treatment2'),sep = ':')
  
  all.pairs.comps <- rbind(temp.pairs, all.pairs.comps)
  temp <-  pairs %>% mutate(File.name = file.name)
  temp <- temp[,c('File.name', 'sex', 'y.name2')] %>%
    separate(col = 'y.name2', into = c('measure.name.MPDname', 'treatment.MPDname'),sep = ':') %>%
    distinct()
  merging.info <- rbind(temp, merging.info)
  
}

pairs.treatment.all <- do.call(rbind, pairs.treatment.all) # within lab t-test of MPD 


dd.means.filenames <- do.call(what = rbind, dd.means.filenames) %>%
  separate(col = 'y.name', into = c('measure.name.MPD', 'treatment'), sep = ':') 
dd.means.filenames <-  dd.means.filenames %>%
  mutate( mean = round(mean, digits = 5),  sd = round(sd, digits = 5))
write_csv(dd.means.filenames)


all.pairs.comps <- all.pairs.comps %>%
  filter(treatment1==treatment2) %>%
  rename('treatment' = 'treatment1', 'lab.MPD'='lab' ) %>%
  dplyr::select(-c('treatment2', 'measure.name0'))
names(all.pairs.comps)[-c(1:7)] <- paste('MPD:', names(all.pairs.comps)[-c(1:7)] ) 
write_csv(all.pairs.comps, 'MPD_all_pairs_analysis.csv')
merging.info <- merging.info %>% mutate(measure.name.MPDname = 
                                          recode(measure.name.MPDname, 
                                                 "Grip.AVG" ="grip.avg",      "OFT1 %C 10m"="OFTsmall_centertime_10m_percent",
                                                 "OFT1 Dist 10m"="OFTsmall_dist_10m_sec", "TST 6m"="tst_7min_percent" ),
                                        sex = recode(sex, 'Females'='Female', 'Males'='Male'),
                                        treatment.MPDname = recode(treatment.MPDname, "control"="Control",    "fluoxetine"="Fluoxetine"))



####################################################
################## Visualisations ##################
####################################################

###### Interaction plots of TAU-JAX data
################### Splitting by Sex and treatment: #############

# pdf('interactions_plot_TAUJAX.pdf', width = 14)

all_data_gxlEstimates <- all_data_gxlEstimates %>%
  mutate(value = round(value, digits = 1))

temp_dat <- subset(all_data_gxlEstimates, measure %in% c("Body.Weight") )

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("Body.Weight")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))

temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin)%>%
  mutate(measure.name='Body Weight')
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

bodyweight_intplot <- all_data_sets %>% filter(measure.name%in%c("Body.Weight")) %>%
  mutate(measure.name='Body Weight') %>%
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')







temp_dat <- subset(all_data_gxlEstimates, measure %in% c("grip.avg") )
set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("grip.avg")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))

temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%
  mutate(measure.name = "Grip Strength (cube root)")
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

gripavg_intplot <- all_data_sets %>% 
  filter(measure.name%in%c( "grip.avg" )) %>%
  mutate(measure.name = "Grip Strength (cube root)") %>%
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')


# dev.off()
# pdf('interactions_plot_TAUJAX_page2.pdf', width = 14)

temp_dat <- subset(all_data_gxlEstimates, measure %in% c("tst_7min_percent") )

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("tst_7min_percent")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))
temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%   
  mutate(measure.name = 'Tail Suspension (logit), 7 min')
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

tst7_intplot <- all_data_sets %>%
  filter(measure.name %in% c("tst_7min_percent")) %>%
  mutate(measure.name = 'Tail Suspension (logit), 7 min') %>%
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free') +
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')

# dev.off()
# pdf('interactions_plot_TAUJAX_page3.pdf', width = 14)

temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTsmall_centertime_10m_percent") ) %>%
  rename('measure.name' = 'measure')
set_bar_begin <- all_data_sets %>% filter(measure.name%in%c( "OFTsmall_centertime_10m_percent")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))
temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin)%>%
  mutate( measure.name = recode(measure.name ,  "OFTsmall_centertime_10m_percent"= 'CT% (logit), small arena, 10 min' )  )
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

CT10_small_intplot <- all_data_sets %>%
  filter(measure.name %in% c( "OFTsmall_centertime_10m_percent")) %>%
  mutate( measure.name = recode(measure.name , "OFTsmall_centertime_10m_percent"= 'CT% (logit), small arena, 10 min' )  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')

temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTlarge_centertime_10m_percent") ) %>%
  rename('measure.name' = 'measure')
set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("OFTlarge_centertime_10m_percent")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))
temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin)%>%
  mutate( measure.name = recode(measure.name , "OFTlarge_centertime_10m_percent" = 'CT% (logit), large arena, 10 min')  )
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

CT10_large_intplot <- all_data_sets %>%
  filter(measure.name %in% c("OFTlarge_centertime_10m_percent")) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_centertime_10m_percent" = 'CT% (logit), large arena, 10 min')  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')

temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTsmall_dist_10m_sec") )%>%
  rename('measure.name' = 'measure')

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c( "OFTsmall_dist_10m_sec")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))

temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%
  mutate( measure.name = recode(measure.name , "OFTsmall_dist_10m_sec" = 'DT, small arena, 10 min' )  )

temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

DT10_small_intplot <- all_data_sets %>%
  filter(measure.name %in% c("OFTsmall_dist_10m_sec")) %>%
  mutate( measure.name = recode(measure.name ,"OFTsmall_dist_10m_sec" = 'DT, small arena, 10 min' )  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')



temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTlarge_dist_10m_sec") )%>%
  rename('measure.name' = 'measure')

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("OFTlarge_dist_10m_sec")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))

temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_dist_10m_sec" = 'DT (Cube root), large arena, 10 min' )  )

temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

DT10_large_intplot <- all_data_sets %>%
  filter(measure.name %in% c("OFTlarge_dist_10m_sec")) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_dist_10m_sec" = 'DT (Cube root), large arena, 10 min' )  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')


temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTsmall_centertime_20m_percent") ) %>%
  rename('measure.name' = 'measure')
set_bar_begin <- all_data_sets %>% filter(measure.name%in%c( "OFTsmall_centertime_20m_percent")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))
temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin)%>%
  mutate( measure.name = recode(measure.name ,  "OFTsmall_centertime_20m_percent"= 'CT% (logit), small arena, 20 min' )  )
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

CT20_small_intplot <- all_data_sets %>%
  filter(measure.name %in% c( "OFTsmall_centertime_20m_percent")) %>%
  mutate( measure.name = recode(measure.name , "OFTsmall_centertime_20m_percent"= 'CT% (logit), small arena, 20 min' )  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')


temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTlarge_centertime_20m_percent") ) %>%
  rename('measure.name' = 'measure')
set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("OFTlarge_centertime_20m_percent")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))
temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin)%>%
  mutate( measure.name = recode(measure.name , "OFTlarge_centertime_20m_percent" = 'CT% (logit), large arena, 20 min')  )
temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

CT20_large_intplot <- all_data_sets %>%
  filter(measure.name %in% c("OFTlarge_centertime_20m_percent")) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_centertime_20m_percent" = 'CT% (logit), large arena, 20 min')  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')


temp_dat <- subset(all_data_gxlEstimates, measure %in% c("OFTlarge_dist_20m_sec") )%>%
  rename('measure.name' = 'measure')

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("OFTlarge_dist_20m_sec")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))


temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_dist_20m_sec" = 'DT, large arena, 20 min')  )

temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

DT20_large_intplot <- all_data_sets %>%
  filter(measure.name %in% c("OFTlarge_dist_20m_sec")) %>%
  mutate( measure.name = recode(measure.name , "OFTlarge_dist_20m_sec" = 'DT, large arena, 20 min')  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')


temp_dat <- subset(all_data_gxlEstimates, measure %in% c( "OFTsmall_dist_20m_sec") )%>%
  rename('measure.name' = 'measure')

set_bar_begin <- all_data_sets %>% filter(measure.name%in%c("OFTsmall_dist_20m_sec")) %>%
  group_by(measure.name, treatment, sex, strain) %>%
  summarise(bar_begin = mean(y)) %>%
  group_by(measure.name, treatment, sex) %>%
  summarise(bar_begin = quantile(bar_begin,p=0.1))

temp_dat <- temp_dat %>% left_join(set_bar_begin) %>%
  mutate(bar_end = value+bar_begin) %>%
  mutate( measure.name = recode(measure.name , "OFTsmall_dist_20m_sec" = 'DT, small arena, 20 min' )  )

temp_dat2 <- temp_dat %>%
  pivot_longer(cols = c(bar_begin, bar_end), names_to = 'bars', values_to = 'bars_lims')

DT20_small_intplot <- all_data_sets %>%
  filter(measure.name %in% c( "OFTsmall_dist_20m_sec")) %>%
  mutate( measure.name = recode(measure.name ,  "OFTsmall_dist_20m_sec" = 'DT, small arena, 20 min' )  ) %>% 
  ggplot(data = .) +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="point") +
  stat_summary(aes(x = strain, y = y, colour = lab, group=lab),fun=mean, geom="line")+
  # geom_col(data = temp_dat, aes(x,value), width = 0.4, alpha=0.5)+
  geom_point(data = temp_dat, aes(x=x, y=bar_begin), shape = 25, fill='black',alpha=1)+
  geom_point(data = temp_dat, aes(x=x, y=bar_end), shape = 24, fill='black',alpha=1)+
  geom_line(data = temp_dat2, aes(x=x, y=bars_lims, group=x), size=1, alpha=1)+
  geom_text(data = temp_dat,aes(x,bar_end*0.98,label = value_type),angle=90, size = 2, hjust = -0.1)+
  geom_text(data = temp_dat,aes(x,bar_begin+value/2,label = value),angle=90, size = 2, vjust= -1)+
  facet_grid(measure.name~treatment+sex, scales = 'free')+
  scale_x_discrete(limits=c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  '))+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Mean measure value')









library(gridExtra)
cols <- c("TAUM"='#F8766D', "TAUL"='#00BA38')



pdf('int_and_boxplots_pg1.pdf', width = 10)
grid.arrange(Bodyweight_boxplot +
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(angle=90), 
                     legend.position = 'top', 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             bodyweight_intplot + scale_color_discrete('') + theme(legend.position = 'none', 
                                                                   strip.text.x = element_blank(),
                                                                   strip.background.x = element_blank(),
                                                                   axis.text.y = element_text(angle=90), 
                                                                   plot.margin = unit(c(0,0,0,.1),'cm')))

dev.off()

pdf('int_and_boxplots_pg2.pdf', width = 10)
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

pdf('int_and_boxplots_pg3.pdf', width = 10)
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

pdf('int_and_boxplots_pg4.pdf', width = 10)
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

pdf('int_and_boxplots_pg5.pdf', width = 10)

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

pdf('int_and_boxplots_pg6.pdf', width = 10)
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

pdf('int_and_boxplots_pg7.pdf', width = 10)
grid.arrange(CTlarge20_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             CT20_large_intplot + theme(legend.position = 'none', 
                                        strip.text.x = element_blank(),
                                        strip.background.x = element_blank(), 
                                        plot.margin = unit(c(0,0,0,.1),'cm'),
                                        axis.text.y = element_text(angle=90)))

dev.off()

pdf('int_and_boxplots_pg8.pdf', width = 10)
grid.arrange(CTsmall20_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             CT20_small_intplot + scale_color_manual('',values = cols) +
               theme(legend.position = 'none', 
                     strip.text.x = element_blank(),
                     strip.background.x = element_blank(), 
                     plot.margin = unit(c(0,0,0,.1),'cm'),
                     axis.text.y = element_text(angle=90)))

dev.off()

pdf('int_and_boxplots_pg9.pdf', width = 10)
cols <- c("TAUM"='#F8766D', "TAUL"='#00BA38', "JAX"='#619CFF')

grid.arrange(DTlarge20_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             DT20_large_intplot + scale_color_manual('',values = cols) +
               theme(legend.position = 'none', 
                     strip.text.x = element_blank(),
                     strip.background.x = element_blank(), 
                     plot.margin = unit(c(0,0,0,.1),'cm'),
                     axis.text.y = element_text(angle=90)))



dev.off()

pdf('int_and_boxplots_pg10.pdf', width = 10)
grid.arrange(DTsmall20_boxplot + 
               theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(), 
                     legend.position = 'top',
                     axis.text.y = element_text(angle=90), 
                     plot.margin = unit(c(0,0,0,.1),'cm') ) ,
             DT20_small_intplot + scale_color_manual('',values = cols) +
               theme(legend.position = 'none', 
                     strip.text.x = element_blank(),
                     strip.background.x = element_blank(), 
                     plot.margin = unit(c(0,0,0,.1),'cm'),
                     axis.text.y = element_text(angle=90)))



dev.off()

pdf('int_and_boxplots_pg11.pdf', width = 10)
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



##########################################################
################## Combine the results ###################
##########################################################



MPD_all_pairs <- all.pairs.comps %>% # read_csv( 'MPD_all_pairs_analysis.csv') %>% 
  rename('strain.1'='strain1', 'strain.2'='strain2')
# gxl_IMPC_all <- read_csv(my.path, 'IMPC_gxl_estimates.csv') %>%
#   dplyr::select(-c('...1')) %>% rename('IMPC_code' = 'X1')

gxl_IMPC <- read_csv(paste0(my.path, 'gxl_estimates_from_IMPC.csv')) %>%
  rename("gxl2.IMPC" = "gamma2 estimator", 'measure.name'='MPD_code') %>%
  dplyr::select(`File name`, measure.name, sex, IMPC_code, gxl2.IMPC)
gxl_IMPC <- gxl_IMPC %>% left_join(gxl_IMPC_all) %>% mutate(gxl2.IMPC = S2_interaction/S2_error )
gxl_IMPC <- gxl_IMPC %>% rename("S2_interaction IMPC"="S2_interaction",     "S2_error IMPC"="S2_error", 'IMPC_lab'='lab') %>%
  # rename('File.name'='File name')%>%
  dplyr::select(# File.name, 
    measure.name, sex, IMPC_code, gxl2.IMPC, `S2_interaction IMPC`, `S2_error IMPC`, df1_IMPC, df2_IMPC)%>%
  distinct() %>%
  mutate(treatment='Control') # %>% filter(!is.na(IMPC_code))

TAUJAX_all_data_analyses <-  all_data_analyses_list %>% # read_csv(my.path, 'TAUJAX_analyses_list.csv') %>%
  # select(-c('strain.1...2','strain.2...3'))%>% rename('strain.1'='strain.1...8', 'strain.2'='strain.2...9') %>%
  dplyr::select(measure.name,treatment,sex,strain.1,strain.2, `Diff Multi`,`t-value REML`,`pval REML`)%>%
  distinct() %>%
  rename('t-value REML TAUJAX'='t-value REML', 'pval REML TAUJAX'='pval REML', 'Diff Multi TAUJAX' = 'Diff Multi')

gxl_TAUJAX <- dfs_TAUJAX %>% # read_csv(my.path, 'TAUJAX_DFs_gxl.csv') %>% 
  rename( 'measure.name' = 'measure.name.MPDname', 'treatment' = 'treatment.MPDname', "gxl2.TAUJAX"="gxl2.fac") %>%
  dplyr::select(measure.name, treatment, sex, df1_TAUJAX, df2_TAUJAX, gxl2.TAUJAX)



gxl_TAUJAX2 <- all_data_gxlEstimates %>% # read_csv(my.path, 'all_data_gxlEstimates.csv') %>% 
  rename('measure.name' = 'measure') %>%
  distinct()

# all_data_gxlEstimates <- all_data_gxlEstimates %>% # read_csv(my.path, 'all_data_gxlEstimates.csv') %>% 
#   rename('measure.name' = 'measure') %>%
#   select(-c("value_type" ,"value",'x')) %>%
#   distinct()

gxl_TAUJAX <- gxl_TAUJAX %>% left_join(gxl_TAUJAX2) # %>% mutate( gxl2.TAUJAX = (int_sd_value/error_sd_value)^2 )
gxl_TAUJAX <- gxl_TAUJAX %>% pivot_wider(id_cols = c("measure.name", "treatment", "sex" ,"df1_TAUJAX", "df2_TAUJAX","gxl2.TAUJAX", "gxl_fac_value","TAU-JAX gxl factor"),
                                         names_from = 'value_type', values_from = 'value') %>%
  mutate( `S2_interaction TAUJAX`= `  Interaction SD`^2, `S2_error TAUJAX`=`   Error SD`^2) %>%
  dplyr::select(measure.name, treatment, sex, df1_TAUJAX, gxl2.TAUJAX,`S2_interaction TAUJAX` ,`S2_error TAUJAX` )

IMPC_all_data_analyses <- rbind(read_csv(my.path, 'IMPC_analyses_list - Males.csv'),
                                read_csv(my.path, 'IMPC_analyses_list - Females.csv')) %>%
  rename(IMPC_code = measure.name)%>%
  dplyr::select(strain.1, strain.2,`Diff Multi`, `t-value REML`, `pval REML` ) %>%
  rename('t-value REML IMPC'='t-value REML', 'pval REML IMPC'='pval REML', 'Diff Multi IMPC' = 'Diff Multi') %>%
  distinct()


IMPC_significance_of_strains <- rbind(read_csv(my.path, 'IMPC_strain_significance_Male.csv') ,
                                      read_csv(my.path, 'IMPC_strain_significance_Female.csv') ) %>%
  rename('signif strain by IMPC'='signif strain')%>% mutate(treatment = 'Control')

gxl_TAUJAX <- gxl_TAUJAX %>% # select(-c('value_type', 'value','x')) %>% 
  distinct()

all_comps_merged_FLOXandCONT <- MPD_all_pairs %>% mutate(sex = recode(sex,  'Males'="Male",   'Females'="Female"),
                                                         treatment = recode(treatment,  "control"="Control","fluoxetine"="Fluoxetine"))%>%
  left_join(TAUJAX_all_data_analyses) %>% 
  left_join(gxl_TAUJAX) %>%
  left_join(gxl_IMPC) %>% 
  left_join(IMPC_significance_of_strains)

all_comps_merged <-  all_comps_merged_FLOXandCONT # %>% filter(  measure.name != "tst_6min_percent" )

all_comps_merged <- all_comps_merged %>% mutate( `MPD: SE.TAUJAX` = sqrt(`MPD: s2.pooled`*(1/`MPD: n1` + 1/`MPD: n2` + 2*gxl2.TAUJAX) ) ,
                                                 `MPD: satter.df.TAUJAX` =  compute.df.sater(s2int = `MPD: s2.pooled`*gxl2.TAUJAX, 
                                                                                             df.int = df1_TAUJAX, 
                                                                                             s2pooled = `MPD: s2.pooled`,n1=`MPD: n1`, n2=`MPD: n2` ),
                                                 `MPD: t-TAUJAX` = abs(`MPD: diff`)/`MPD: SE.TAUJAX` ) %>%
  mutate(`MPD: gxl.adjusted.pv.TAUJAX` = 2*pt(abs(`MPD: t-TAUJAX`) , df = `MPD: satter.df.TAUJAX`, lower.tail = F ) )

all_comps_merged <- all_comps_merged %>% dplyr::select(measure.name, strain.1, strain.2, sex, lab.MPD,  treatment, `MPD: diff`, `MPD: satter.df.IMPC`, `MPD: satter.df.TAUJAX`, `MPD: gxl.adjusted.pv.IMPC`, `MPD: gxl.adjusted.pv.TAUJAX`,
                                                       `MPD: common.pv`, `Diff Multi TAUJAX`, `MPD: t-TAUJAX` , `t-value REML TAUJAX`, `pval REML TAUJAX`, IMPC_code, `signif strain by IMPC`, File.name)
# write_csv(all_comps_merged,file = 'MPD_crossed_analysis.csv')

all_comps_merged$`MPD: gxl.adjusted.pv.IMPC`[is.na(all_comps_merged$`MPD: gxl.adjusted.pv.IMPC`)] <-
  all_comps_merged$`MPD: gxl.adjusted.pv.TAUJAX`[is.na(all_comps_merged$`MPD: gxl.adjusted.pv.IMPC`)]

all_comps_merged$`MPD: satter.df.IMPC`[is.na(all_comps_merged$`MPD: satter.df.IMPC`)] <-
  all_comps_merged$`MPD: satter.df.TAUJAX`[is.na(all_comps_merged$`MPD: satter.df.IMPC`)]

all_comps_merged2 <- all_comps_merged

all_comps_merged <- all_comps_merged %>% 
  mutate(`True Significance (by TAUJAX)` = recode( (`Diff Multi TAUJAX` >= 0)*1 , '1'='Positive difference', '0'='Negative Difference' ),
         `Significance in MPD` = recode( ( `MPD: diff`>=0 )*1, '1'='Positive difference',  '0'='Negative Difference' )) %>%
  mutate(`Significance in MPD with GXL (from IMPC)` = `Significance in MPD`,
         `Significance in MPD with GXL (from TAUJAX)` = `Significance in MPD`) %>%
  
  mutate(`True Significance (by TAUJAX)` = replace(`True Significance (by TAUJAX)`, (`pval REML TAUJAX`> alpha.threshold), 'No difference' ),
         `Significance in MPD` = replace(`Significance in MPD`, (`MPD: common.pv`> alpha.threshold), 'No difference' ),
         `Significance in MPD with GXL (from TAUJAX)` = replace(`Significance in MPD with GXL (from TAUJAX)` , 
                                                                (`MPD: gxl.adjusted.pv.TAUJAX`> alpha.threshold),
                                                                'No difference' ),
         `Significance in MPD with GXL (from IMPC)` = replace(`Significance in MPD with GXL (from IMPC)` , 
                                                              (`MPD: gxl.adjusted.pv.IMPC`> alpha.threshold),
                                                              'No difference' ) ) 

all_comps_merged <- all_comps_merged %>% 
  mutate(`Significance in MPD with GXL (from TAUJAX)` = recode(`Significance in MPD with GXL (from TAUJAX)` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ),
         `Significance in MPD with GXL (from IMPC)` = recode(`Significance in MPD with GXL (from IMPC)` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ),
         `True Significance (by TAUJAX)` = recode(`True Significance (by TAUJAX)` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ),
         `Significance in MPD` = recode(`Significance in MPD` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ))
all_comps_merged <- all_comps_merged %>% 
  mutate(`True Significance (by TAUJAX)` = factor(`True Significance (by TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from TAUJAX)` = factor(`Significance in MPD with GXL (from TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from IMPC)` = factor(`Significance in MPD with GXL (from IMPC)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD` = factor(`Significance in MPD`, levels = c("Sig difference",'No difference' )))

temp <- all_comps_merged %>% filter(measure.name!="tst_7min_percent", treatment == "Control") %>% 
  group_by(`True Significance (by TAUJAX)`,
           `Significance in MPD`,
           `Significance in MPD with GXL (from IMPC)`) %>%
  summarise(n=n())
temp.template <- temp[,1:3]

temp <- temp.template %>% left_join(temp) %>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write.csv(temp , 'rejections (gxl by IMPC, no TST, Control).csv')


temp <- all_comps_merged %>% filter(measure.name!="tst_7min_percent", treatment ==  "Fluoxetine") %>%
  group_by(`True Significance (by TAUJAX)`,
           `Significance in MPD`,
           `Significance in MPD with GXL (from IMPC)`) %>%
  summarise(n=n()) 
temp <- temp.template %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write.csv(temp, 'rejections (gxl by IMPC, no TST, Floxitine).csv')

temp <- all_comps_merged %>% filter(measure.name!="tst_7min_percent") %>%
  group_by(`True Significance (by TAUJAX)`,
           `Significance in MPD`,
           `Significance in MPD with GXL (from IMPC)`) %>%
  summarise(n=n()) 
temp <- temp.template %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)

write.csv(temp,'rejections (gxl by IMPC, no TST, Both treatments).csv')


temp <- all_comps_merged %>%
  group_by( `True Significance (by TAUJAX)`,
            `Significance in MPD`,
            `Significance in MPD with GXL (from IMPC)`) %>% 
  summarise( n = n() ) 

temp <- temp.template %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write.csv(temp, 'rejections (gxl by IMPC, with TST, Both treatments).csv')


temp <- all_comps_merged %>% filter(measure.name!="tst_7min_percent") %>%
  group_by(`True Significance (by TAUJAX)`,
           `Significance in MPD`,
           `Significance in MPD with GXL (from TAUJAX)`) %>%
  summarise(n=n())  
temp.template <- temp[,1:3]

temp <- temp.template.tj %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from TAUJAX)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)

write.csv(temp, 'rejections (gxl by TAUJAX, no TST).csv')


temp <- all_comps_merged %>% filter(measure.name!="tst_7min_percent", treatment ==  "Control") %>%
  group_by(`True Significance (by TAUJAX)`,
           `Significance in MPD`,
           `Significance in MPD with GXL (from TAUJAX)`) %>%
  summarise(n=n())  

temp <- temp.template.tj %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from TAUJAX)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)

write.csv(temp, 'rejections (gxl by TAUJAX, no TST, Control).csv')


##### replacing adjustment by lower alpha threshold

all_comps_merged2 <- all_comps_merged2 %>% 
  mutate(`True Significance (by TAUJAX)` = recode( (`Diff Multi TAUJAX` >= 0)*1 , '1'='Positive difference', '0'='Negative Difference' ),
         `Significance in MPD` = recode( ( `MPD: diff`>=0 )*1, '1'='Positive difference',  '0'='Negative Difference' )) %>%
  mutate(`Significance in MPD with rejection at 0.5%` = `Significance in MPD`) %>%
  
  mutate(`True Significance (by TAUJAX)` = replace(`True Significance (by TAUJAX)`, (`pval REML TAUJAX`> alpha.threshold), 'No difference' ),
         `Significance in MPD` = replace(`Significance in MPD`, (`MPD: common.pv`> alpha.threshold), 'No difference' ),
         `Significance in MPD with rejection at 0.5%` = replace(`Significance in MPD with rejection at 0.5%` , 
                                                                (`MPD: gxl.adjusted.pv.TAUJAX`> alpha.threshold2),
                                                                'No difference' )) 

all_comps_merged2 <- all_comps_merged2 %>% 
  mutate(`True Significance (by TAUJAX)` = recode(`True Significance (by TAUJAX)` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ),
         `Significance in MPD` = recode(`Significance in MPD` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ),
         `Significance in MPD with rejection at 0.5%` = recode(`Significance in MPD with rejection at 0.5%` , "Negative Difference"="Sig difference", "Positive difference"="Sig difference" ))

all_comps_merged2 <- all_comps_merged2 %>% 
  mutate(`True Significance (by TAUJAX)` = factor(`True Significance (by TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with rejection at 0.5%` = factor(`Significance in MPD with rejection at 0.5%`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD` = factor(`Significance in MPD`, levels = c("Sig difference",'No difference' )))


temp <- all_comps_merged2 %>%
  group_by( `True Significance (by TAUJAX)`,
            `Significance in MPD`,
            `Significance in MPD with rejection at 0.5%`) %>% 
  summarise( n = n() ) 

temp.template.tj <- temp.template <- temp[,1:3]
names(temp.template.tj)[3] <- 'Significance in MPD with rejection at 0.5%'
temp <- temp.template.tj %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with rejection at 0.5%` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_nominal_alpha = NA, typeIError_10nt_nominal_alpha = NA,
                        Power_nominal_alpha = NA, Power_10nt_nominal_alpha = NA,
                        FDP_nominal_alpha = NA, FDP_10nt_nominal_alpha = NA)%>%
  as.data.frame()
temp[1,c('typeIError_nominal_alpha','Power_nominal_alpha','FDP_nominal_alpha')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_10nt_nominal_alpha','Power_10nt_nominal_alpha','FDP_10nt_nominal_alpha')] <- c(tIe2,Pwr2,FDP2)
write.csv(temp, 'rejections (gxl by  TAUJAX, with TST: adjustment replaced by lower alpha).csv')

# sum(temp$n)

all_comps_merged %>%
  dplyr:: select(measure.name, strain.1, strain.2, sex, lab.MPD, treatment, `MPD: diff`, 
                 `MPD: gxl.adjusted.pv.IMPC`, `MPD: gxl.adjusted.pv.TAUJAX`, `MPD: common.pv`, `Diff Multi TAUJAX`, # `MPD: t-TAUJAX`,
                 `pval REML TAUJAX`, `signif strain by IMPC`, File.name, `True Significance (by TAUJAX)`,
                 `Significance in MPD`, `Significance in MPD with GXL (from IMPC)`, `Significance in MPD with GXL (from TAUJAX)`) %>%
  mutate(`MPD: diff` = round(`MPD: diff`, digits = 3),
         `MPD: common.pv` = round(`MPD: common.pv`, digits = 3),
         `MPD: gxl.adjusted.pv.IMPC`= round(`MPD: gxl.adjusted.pv.IMPC`, digits = 3),
         `MPD: gxl.adjusted.pv.TAUJAX`= round(`MPD: gxl.adjusted.pv.TAUJAX`, digits = 3),
         `pval REML TAUJAX`=round(`pval REML TAUJAX`, digits = 3) ) %>%
  rename( 'pvalue original' = 'MPD: common.pv', 
          'pvalue with gxl (from IMPC)'='MPD: gxl.adjusted.pv.IMPC', 
          'pvalue with gxl (from TAUJAX)'='MPD: gxl.adjusted.pv.TAUJAX', 
          'MPD lab name'='lab.MPD')


################## 
## Crossing treatment related contrasts
# treatment effect within strain: 
pairs.treatment.all <- # MPD pairwise comparisons of treatment effect per strain
  pairs.treatment.all %>% 
  left_join(all.contrasts2) %>% # post hoc comparisons of treatment effect per strain
  left_join(gxl_IMPC) %>%
  mutate(SE.IMPC = sqrt(  s2.pooled*(1/n1 + 1/n2 ) + 2*`S2_interaction IMPC` )) %>%
  mutate(satter.df.IMPC = compute.df.sater(s2int = `S2_interaction IMPC`,
                                           df.int = df2_IMPC,
                                           s2pooled = s2.pooled,  n1 = n1 ,  n2 = n2) ) %>%
  mutate(gxl.adjusted.pv.IMPC = 2*pt( abs(diff)/SE.IMPC  , df = satter.df.IMPC, lower.tail = F)) %>% 
  mutate(`True Significance (by TAUJAX)` = recode( (`REML post-hoc` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
         `Significance in MPD` = recode( (common.pv <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
         `Significance in MPD with GXL (from IMPC)` = recode( (gxl.adjusted.pv.IMPC <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
         `Significance in MPD with GXL (from TAUJAX)` = recode( (gxl.adjusted.pv.TAUJAX <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' )) %>%
  mutate(`True Significance (by TAUJAX)` = factor(`True Significance (by TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from TAUJAX)` = factor(`Significance in MPD with GXL (from TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from IMPC)` = factor(`Significance in MPD with GXL (from IMPC)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD` = factor(`Significance in MPD`, levels = c("Sig difference",'No difference' )))

temp <- pairs.treatment.all %>% group_by(`True Significance (by TAUJAX)`, `Significance in MPD`,
                                 `Significance in MPD with GXL (from IMPC)`) %>% summarise(n=n())
temp <- temp.template %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write_csv(temp, 'rejections_treatmen_per_strain(gxl by IMPC).csv')

temp <- pairs.treatment.all %>% group_by(`True Significance (by TAUJAX)`, `Significance in MPD`,
                                 `Significance in MPD with GXL (from TAUJAX)` ) %>% summarise(n=n()) 

temp <- temp.template.tj %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write_csv(temp, 'rejections_treatmen_per_strain(gxl by TAUJAX).csv')


# treatment effect between strains: 
pairs.all <- # MPD pairwise comparisons of treatment effect between pairs of strains
  pairs.all %>%
  left_join(all.contrasts) %>% # post hoc comparisons of treatment effect between strains
  left_join(gxl_IMPC) %>%
  mutate( `diff.se2.gxladj (by IMPC)` = s2.pooled * n.op + 4*`S2_interaction IMPC`^2 ) %>%
  mutate(`t-test gxl ajd (by IMPC)` = diff / sqrt(`diff.se2.gxladj (by IMPC)`),
          `df.satt (by IMPC)` = compute.df.sater.2pairs(s2int = `S2_interaction IMPC`, df.int = df1_IMPC,
                                            s2pooled = s2.pooled, n1 = n1L, n2 = n1R, n3 = n2L, n4 = n2R ) ) %>%
  mutate( `p-value gxl ajd (by IMPC)` = pt( abs(`t-test gxl ajd (by IMPC)`), df = `df.satt (by IMPC)`, lower.tail = F  )*2) %>%
  mutate(`True Significance (by TAUJAX)` = recode( (`REML post-hoc` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
       `Significance in MPD` = recode( (`p-value.common` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
       `Significance in MPD with GXL (from IMPC)` = recode( (`p-value gxl ajd (by IMPC)` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
       `Significance in MPD with GXL (from TAUJAX)` = recode( (`p-value gxl ajd (by TAUJAX)` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' )) %>%
  mutate(`True Significance (by TAUJAX)` = factor(`True Significance (by TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from TAUJAX)` = factor(`Significance in MPD with GXL (from TAUJAX)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD with GXL (from IMPC)` = factor(`Significance in MPD with GXL (from IMPC)`, levels = c("Sig difference",'No difference' )),
         `Significance in MPD` = factor(`Significance in MPD`, levels = c("Sig difference",'No difference' )))

temp <- pairs.all %>% group_by(`True Significance (by TAUJAX)`, `Significance in MPD`,
                                 `Significance in MPD with GXL (from IMPC)`) %>% summarise(n=n()) 
temp <- temp.template %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write_csv(temp,'rejections_treatment_between_strains(gxl by IMPC).csv')


temp <- pairs.all %>% group_by(`True Significance (by TAUJAX)`, `Significance in MPD`,
                                 `Significance in MPD with GXL (from TAUJAX)` ) %>% summarise(n=n())
temp <- temp.template.tj %>% left_join(temp)%>%
  mutate(n = replace(n, is.na(n) , 0 ))

x1 <- (temp$`True Significance (by TAUJAX)` == 'Sig difference')
x2 <- (temp$`Significance in MPD` == 'Sig difference')
x3 <- (temp$`Significance in MPD with GXL (from IMPC)` == 'Sig difference')

A <- temp$n[x1 & x2 & x3 ]
B <- temp$n[(x1 & x2 & (!x3))]
C <- temp$n[(x1 & (!x2) & (!x3))]
D <- temp$n[((!x1) & (x2) & (x3))]
E <- temp$n[(!x1) & (x2) & (!x3)]
F1 <- temp$n[(!x1) & (!x2) & (!x3)]

A <- ifelse( length(A) == 0 , 0, A )
B <- ifelse( length(B) == 0 , 0, B )
C <- ifelse( length(C) == 0 , 0, C )
D <- ifelse( length(D) == 0 , 0, D )
E <- ifelse( length(E) == 0 , 0, E )
F1 <- ifelse( length(F1) == 0 , 0, F1 )

tIe1 <- (D+E)/(D+E+F1) ; tIe2 <- (D)/(D+E+F1)
Pwr1 <- (A+B)/(A+B+C)  ; Pwr2 <- (A)/(A+B+C)
FDP1 <- (D+E)/(D+E+A+B); FDP2 <- (D)/(D+A)


temp <- temp %>% mutate(BLANC=' ', 
                        typeIError_no_GxL = NA, typeIError_with_GxL = NA,
                        Power_no_GxL = NA, Power_with_GxL = NA,
                        FDP_no_GxL = NA, FDP_with_GxL = NA)%>%
  as.data.frame()
temp[1,c('typeIError_no_GxL','Power_no_GxL','FDP_no_GxL')] <- c(tIe1,Pwr1,FDP1)
temp[1,c('typeIError_with_GxL','Power_with_GxL','FDP_with_GxL')] <- c(tIe2,Pwr2,FDP2)
write_csv(temp,'rejections_treatment_between_strains(gxl by TAUJAX).csv')











##################### Perform analysis of Wiltshire2 with 28 strains, without ground truth, but inter the gxl from tau-jax. 
# report the full table/endpoints separate by affected by taujax

file.names <- c("Lab_Wiltshire2_DistanceTraveled_Males.csv",
                "Lab_Wiltshire2_PercentCenter_Males.csv",
                "Lab_Wiltshire2_TST_Males.csv" )

pairs.all.WILTSHIRE2 <-  list()
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
    mutate(treatment = recode(treatment, 'fluoxetine'='Fluoxetine', 'control'='Control'))
  
  dd.means <- dd %>%
    group_by(strain, lab, sex, y.name, treatment) %>% 
    summarise(mean = mean(y,na.rm = T) , sd = sd(y,na.rm = T) , n = n() )
  
  
  pairs.right <- expand_grid( strain1 = as.character(unique( dd$strain)), 
                              y.name = unique(dd$y.name), 
                              treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
    filter( treatment1 < treatment2 , !is.na(strain1) )
  
  pairs.right <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment1','y.name', 'mean1R', 'sd1R', 'n1R')) %>%
    right_join(pairs.right)
  
  pairs.right <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain2', 'treatment2', 'y.name', 'mean2R', 'sd2R', 'n2R')) %>%
    right_join(pairs.right)
  
  
  pairs.left <- expand_grid( strain2 = as.character(unique( dd$strain)), 
                             y.name = unique(dd$y.name), 
                             treatment1 = unique(dd$treatment), treatment2 = unique(dd$treatment) ) %>%
    filter( treatment1 < treatment2 , !is.na(strain2) )
  
  pairs.left <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment1', 'y.name', 'mean1L', 'sd1L', 'n1L')) %>%
    right_join(pairs.left)
  
  pairs.left <- dd.means %>%
    rename_at(vars(c('strain','treatment', 'y.name','mean', 'sd', 'n')), ~ c('strain1', 'treatment2', 'y.name','mean2L', 'sd2L', 'n2L')) %>%
    right_join(pairs.left)
  
  pairs.right <- pairs.right %>% 
    mutate(diff.R = ( mean1R - mean2R ) )
  
  pairs.left <- pairs.left %>% 
    mutate(diff.L = ( mean1L - mean2L ) )
  
  pairs.all.WILTSHIRE2[[file.name]] <- full_join(pairs.right, pairs.left) %>%
    filter(strain1 > strain2) %>% 
    rename('measure.name'='y.name' ) %>%
    filter(!((strain1 %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') )&(strain2 %in% c('BTBR', 'BALB/cJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'SWR/J') )) )
  
}
pairs.all.WILTSHIRE2 <- do.call(rbind, pairs.all.WILTSHIRE2)  %>%
  mutate(sex = recode(sex, 'Males' = 'Male', 'Females' = 'Female')) %>%
  left_join(threewya.aov.resuls.TAUJAX) %>%
  mutate(diff = (diff.L - diff.R),
         s2.pooled = ( (n1L-1)*sd1L^2 + (n1R-1)*sd1R^2 + (n2L-1)*sd2L^2 + (n2R-1)*sd2R^2 ) / ( (n1L-1) + (n1R-1) + (n2L-1) + (n2R-1) ),
         n.op = ( (n1L^-1) + (n1R^-1) + (n2L^-1) + (n2R^-1) ) ) %>%
  mutate( diff.se2 = s2.pooled * (n.op),
          `diff.se2.gxladj (by TAUJAX)` = s2.pooled * (n.op + 4*trt.str.lab.s2^2 )) %>%
  mutate( `t-test` = diff / sqrt(diff.se2),
          `t-test gxl ajd (by TAUJAX)` = diff / sqrt(`diff.se2.gxladj (by TAUJAX)`),
          `df.satt (by TAUJAX)` = compute.df.sater.2pairs(s2int = trt.str.lab.s2, df.int = (nstrain-1)*(nlab-1)*(2-1),
                                                          s2pooled = s2.pooled, n1 = n1L, n2 = n1R, n3 = n2L, n4 = n2R ) ) %>%
  mutate( `p-value.common` = pt( abs(`t-test`), df = n1L + n1R + n2L + n2R - 4, lower.tail = F  )*2,
          `p-value gxl ajd (by TAUJAX)` = pt( abs(`t-test gxl ajd (by TAUJAX)`), df = `df.satt (by TAUJAX)`, lower.tail = F  )*2) %>%
  mutate(`Significance in MPD` = recode( (`p-value.common` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ),
         `Significance in MPD with GXL (from TAUJAX)` = recode( (`p-value gxl ajd (by TAUJAX)` <= alpha.threshold)*1 , '1'='Sig difference', '0'='No difference' ))

pairs.all.WILTSHIRE2 %>% group_by( `Significance in MPD`,
                                   `Significance in MPD with GXL (from TAUJAX)`) %>% summarise(n=n()) %>%
  write_csv('rejections_treatment_between_strains (30 strains at WILTSHIRE2).csv')



pairs.all.WILTSHIRE2 %>% group_by( `Significance in MPD`,
                                   `Significance in MPD with GXL (from TAUJAX)`) %>% summarise(n=n()) %>%
  write_csv('rejections_treatment_between_strains (30 strains at WILTSHIRE2).csv')

pairs.all.WILTSHIRE2 %>% group_by( measure.name, `Significance in MPD`,
                                   `Significance in MPD with GXL (from TAUJAX)`) %>% summarise(n=n()) %>%
  write_csv('rejectionsPERendpoint_treatment_between_strains (30 strains at WILTSHIRE2).csv')



pairs.all.WILTSHIRE2 %>% select(strain1, strain2, measure.name, `LRT pvalue of lab-strain-treatment`,
                                `LRT pvalue of strain-treatment`, diff, `p-value.common`,`p-value gxl ajd (by TAUJAX)`,
                                `Significance in MPD`,`Significance in MPD with GXL (from TAUJAX)`) %>%
  write_csv('analysis_treatment_between_strains (30 strains at WILTSHIRE2).csv')

  



# 4point2 plot
threewya.aov.resuls.TAUJAX <- read_csv('results_of_three_way_anova_at_TAUJAX.csv') %>% 
  mutate(treatment = 'Either') %>%
  full_join(gxl_IMPC) %>%
  full_join(gxl_TAUJAX) %>%
  select(measure.name, sex, `gamma of lab-strain`,  gxl2.IMPC, gxl2.TAUJAX, treatment) %>%
  mutate(`gamma of lab-strain (3-way analysis)` = `gamma of lab-strain`,
         `gamma of lab-strain (2-way analysis, TAU-JAX)` =sqrt(gxl2.TAUJAX),
         `gamma of lab-strain (2-way analysis, IMPC)` =sqrt(gxl2.IMPC) ) %>%
  pivot_longer(cols = c('gamma of lab-strain (3-way analysis)','gamma of lab-strain (2-way analysis, TAU-JAX)','gamma of lab-strain (2-way analysis, IMPC)'), 
               names_to = 'Source', values_to = 'Gamma estimate') %>%
  mutate( Source = recode(Source ,
                          'gamma of lab-strain (3-way analysis)' = '3-way analysis',
                          'gamma of lab-strain (2-way analysis, TAU-JAX)' = '2-way analysis, TAU-JAX (Control)',
                          'gamma of lab-strain (2-way analysis, IMPC)' = '2-way analysis, IMPC')) %>%
  mutate(measure.name = paste( measure.name , '('),
         sex = paste(sex, ')'),
         Source = replace(Source, treatment == "Fluoxetine", "2-way analysis, TAU-JAX (Fluoxetine)")) %>%
  mutate(measure.name = paste(measure.name, sex))%>%
    filter(!is.na(`Gamma estimate`))

pdf('4point2.pdf', width = 12)
ggplot(threewya.aov.resuls.TAUJAX, aes(x = measure.name , y = `Gamma estimate`, 
                                       color =  Source , shape = Source, size = Source ) ) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'bottom') +
  scale_shape_manual(values = c(0,1,2,5)) +
  scale_size_manual(values = 2*c(1.5,2,2.5,0.8))+
  scale_color_manual(values = c(6,3,2,4))+ 
  coord_flip()+
  xlab(NULL)
dev.off()
