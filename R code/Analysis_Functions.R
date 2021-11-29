#### Functions for 3-lab replicability trial ####
require('dplyr')
require('tidyverse')
require('tidyr')
require('lmerTest')
require('lme4')
require('readxl')
library('emmeans')
library('nlme')



## Function 1: Tidy
# Arguments:  data file, the demographic variables and the name of the analysed measure. 
# value: A list of tibbles split by test name and (mainly) treatment.  
clean_dat <- function( messy.data = NULL, strain = 'strain', lab = 'lab', sex = 'sex', desig.variable = NULL ){
  if( is.numeric(desig.variable) ){ 
    message('The designated variable name was supplied as column number. It was turned to the actual name.')
    desig.variable <- names(messy.data)[desig.variable]
  }
  
  if( 'gender' %in% names(messy.data) ){
    messy.data <- messy.data %>% mutate(sex = gender)
  }
  if( ! sex %in% names(messy.data) ){
    messy.data <- messy.data %>% mutate(sex = 'unspecified')
  }
  
  
  oldnames <- c(strain,   lab,    sex,   'y.name',   desig.variable)
  newnames <- c('strain', 'lab', 'sex', 'y.name', 'y')
  
  sub_data <- messy.data %>% rename_at(vars(oldnames), ~ newnames) %>%
    dplyr::select(strain, lab, sex, y.name, y ) 
  if( is.factor(sub_data$lab)){ 
    sub_data$lab <- droplevels(sub_data$lab)
  }else{
    sub_data$lab <- parse_factor(sub_data$lab)
  }
  
  sub_data$strain <- as.character(sub_data$strain)
  sub_data <- sub_data %>% filter(!is.na(strain))
  sub_data$strain <- as.factor(sub_data$strain)
  if( ! is.character(sub_data$sex)){ 
    sub_data$sex <- as.character(sub_data$sex, trim_ws = T )
  }
  
  
  if( ! is.numeric(sub_data$y)|is.double(sub_data$y) ){
    sub_data$y <- as.numeric(sub_data$y)
  }
  
  sub_data <- sub_data %>% 
    filter_all(.vars_predicate = all_vars(!is.na(.))) %>%
    filter( ! strain %in% c( '', ' '), ! lab%in% c( '', ' '), ! sex %in% c( '', ' '), ! y.name%in% c( '', ' ') ) 
  
  if( nrow(sub_data)==0 ){
    stop('Check the data structure. Are strain, lab, sex and the measure variable supplied  propperly?')
  }
  if( nrow(messy.data) < nrow(sub_data) ){
    message( paste0('Number of rows removed due to missing values: ', nrow(messy.data) - nrow(sub_data)  ) )
    message( paste0('Active number of samples: ', nrow(sub_data)  ) )
  }
  
  if(n_distinct(sub_data[,'sex']) > 2 ) {
    warning('Sex variable contains over 2 sex y.names.')
  }
  if(n_distinct(sub_data[,'y.name']) > 1 ) {
    message('There are at least two treatments')
    sub_data <- sub_data %>% group_split(y.name)
  }
    return(sub_data)
}


## Function 2: Perform Mixed model analysis
# Arguments:  
# value: 
multi_lab_analysis <- function( tidy_data = NULL,
                                desig.variable = 'y',
                                trans = function(y){y},
                                mult.comp.method = 'none', 
                                report.full.summary = FALSE) {
  # fit.lme <- lme(y ~ strain, random = ~ 1|lab/strain,
  #                data = tidy_data,
  #                method="REML",na.action="na.omit")
  if(is.list(tidy_data)&(!is.tbl(tidy_data)) ){
    stop('The data must be tibble, table or data.frame')
  }
  if( all (c(n_distinct( tidy_data$y.name ) > 1, n_distinct( tidy_data$strain) > 1) )){
    stop('Multiple treatments and multiple strains are reported. make sure that one of them must have a single level.')
  }
  
  tidy_data <- tidy_data %>% mutate(strain.orig = strain)
  switchstraintreat <- ( n_distinct( tidy_data$y.name ) > 1)
  
  means_table <- tidy_data %>% 
    separate(col = 'y.name', into = c('measure.name', 'treatment', 'sex'),sep = ':',   remove = T ) %>%
    group_by(treatment, lab, sex, strain, measure.name) %>%
    summarise(meann = mean(y, na.rm = T), sdd = sd(y, na.rm = T), n=n()) %>% 
    mutate(meann = round(meann, digits = 5), sdd = round(sdd, digits = 5)) %>% 
    dplyr::select(lab, measure.name,sex,strain, treatment, meann,   sdd, n)

  if( switchstraintreat ) {
    tidy_data <- tidy_data %>% 
      separate(col = 'y.name', into = c('measure.name', 'strain'),sep = ':',   remove = F )
  }
  
  
  fit.lmer <-  lmer(y ~ strain + (1|lab/strain),
                 data = tidy_data, REML = T)
  vars.mat <- as.data.frame(VarCorr(fit.lmer))
  sigma_e_hat <- vars.mat[3,5]
  gamma_hat  <- vars.mat[1,5]/ vars.mat[3,5]
  sigma_int_hat <- vars.mat[1,5]

  signif_pairs <- multi_lab_PostHoc(Multi.lab.fit = fit.lmer,  
                                    mult.comp.method = mult.comp.method,
                                    data.summary = means_table )
  

  mlf <- multi_lab_alti_truth(Multi.lab.fit.summary = signif_pairs)
  lsp <- signif_pairs
  
  if( switchstraintreat) {
    lsp <- signif_pairs %>% mutate( treatment1 = strain1, treatment2 = strain2 ) %>% 
      mutate( strain1 = tidy_data$strain.orig[1], strain2 = tidy_data$strain.orig[1])
    mlf <- mlf %>% mutate( treatment1 = strain1, treatment2 = strain2 ) %>% 
      mutate( strain1 = tidy_data$strain.orig[1], strain2 = tidy_data$strain.orig[1])
  }
  
  if(report.full.summary){
    return( c(`Altimum Truth` = mlf ,
                  `Post-Hoc Results` = list(lsp), 
              `summary_table` = list(means_table),
              gxl_factor = gamma_hat,
              error_sd = sigma_e_hat,
              int_sd = sigma_int_hat) )
  }
  return( list( `Altimum Truth` = mlf,
                `summary_table` = means_table,
                gxl_factor = gamma_hat,
                error_sd = sigma_e_hat,
                int_sd = sigma_int_hat ) )
  }


## Function 3: Post-hoc analysis for pairwise comparisons
# Arguments:  
# value: 
multi_lab_PostHoc <- function( Multi.lab.fit = NULL,  mult.comp.method = 'none',
                               data.summary = NULL ) {
  summary.tibble <- difflsmeans( Multi.lab.fit , test.effs = "strain")

  pairs <- tibble( pairs = rownames(summary.tibble) ) %>%
      separate(pairs, c('strain.1', "strain.2"), sep = "([-])") %>%
      mutate(strain.1 = str_remove_all(strain.1, "[ strain]"), 
             strain.2 = str_remove_all(strain.2, "[ strain]") ) %>%
    mutate(reordered = (strain.1 < strain.2)) 
  
  
  summary.tibble <- as.tibble(summary.tibble) %>%
    mutate(strain.1 = pairs$strain.1, strain.2 = pairs$strain.2 ) %>%
    mutate(pvalue.adjusted = p.adjust(`Pr(>|t|)`,method = mult.comp.method) ) 
  summary.tibble[pairs$reordered,c('strain.1', 'strain.2', 'lower', 'upper' ) ] <-
    summary.tibble[pairs$reordered,c('strain.2', 'strain.1', 'upper' , 'lower' ) ]
  summary.tibble[pairs$reordered,c('Estimate', 't value', 'lower', 'upper' ) ] <-
    -1 * summary.tibble[pairs$reordered,c('Estimate', 't value', 'upper', 'lower' ) ]
  
  #   mutate(strain.first = NA, strain.second = NA, reordered.strain = F)
  # 
  # summary.tibble$strain.first <- apply(summary.tibble %>% select(strain.1, strain.2),
  #                         MARGIN = 1, FUN = function(x) sort(x)[1])
  # summary.tibble$strain.second <- apply(summary.tibble %>% select(strain.1, strain.2),
  #                                      MARGIN = 1, FUN = function(x) sort(x)[2])
  # summary.tibble <- summary.tibble %>%
  #   mutate( reordered.strain = (strain.first != strain.1) ) %>%
  #   mutate(  Estimate =  Estimate * (-1)^reordered.strain, 
  #           `t value` = `t value` * (-1)^reordered.strain, 
  #           strain.1 = strain.first, strain.2 = strain.second ) %>%
  summary.tibble <-  summary.tibble %>% 
    dplyr::select(strain.1, strain.2, Estimate, `Std. Error`, df, `t value`, `Pr(>|t|)`, pvalue.adjusted )
  
  
  if(!is.null(data.summary)){
    data.summary <- data.summary %>% ungroup()
    summary.tibble <- data.summary %>% 
      mutate(strain.1 = strain, sd.1 = sdd, mean.1 = meann, n.1 = n ) %>% 
      dplyr::select(-c('strain', 'sdd', 'meann', 'n')) %>%
      right_join( summary.tibble )
    
    
    summary.tibble <- data.summary  %>% 
      mutate(strain.2 = strain, sd.2 = sdd, mean.2 = meann, n.2 = n  ) %>%  
      dplyr::select(-c('strain', 'sdd', 'meann', 'n')) %>%
      right_join( summary.tibble ) %>%
      mutate(Diff = mean.1 - mean.2) %>%
      rename( 'Diff Multi'="Estimate", "Std. Error REML"="Std. Error", "df REML"="df", "t-value REML"="t value",
              "pval REML"="Pr(>|t|)","Mult. comps pvalue REML"="pvalue.adjusted") %>%
      mutate(s.pooled = sqrt(((n.1-1)*sd.1^2 + (n.2-1)*sd.2^2)/(n.1+n.2-2)) , 
             `t-value common` = Diff / (s.pooled*sqrt(1/n.1 + 1/n.2 ))) %>%
      mutate(`pval common` = 2*pt(abs(`t-value common`), df = n.1+n.2-2, lower.tail = F ))
    
    summary.tibble <- summary.tibble %>%
      dplyr::select( treatment,lab,sex, strain.1, strain.2, mean.1, mean.2, sd.1, sd.2, n.1, n.2, , Diff, `Diff Multi`, s.pooled, `t-value common`,`pval common`,
                     `Std. Error REML`, `df REML`,`t-value REML`, `pval REML`,`Mult. comps pvalue REML`)
    
  }
  
  
  return(summary.tibble)
}

## Function 4: Ground truth based on Post-hoc analysis
# Arguments:  
# value: 
multi_lab_alti_truth <- function( Multi.lab.fit.summary = NULL ){
  unadjusted.comps <- Multi.lab.fit.summary # %>% #filter( `Pr(>|t|)` < 0.05 ) %>%
    # mutate(strain.first = NA, strain.second = NA, reordered.strain = F)
  if(nrow(unadjusted.comps)>0){
    unadjusted.comps$strain.first <- apply(unadjusted.comps %>% dplyr::select(strain.1, strain.2),
                                         MARGIN = 1, FUN = function(x) sort(x)[1])
    unadjusted.comps$strain.second <- apply(unadjusted.comps %>% dplyr::select(strain.1, strain.2),
                                          MARGIN = 1, FUN = function(x) sort(x)[2])
    # unadjusted.comps <- unadjusted.comps %>%
    #   dplyr::select(strain.1, strain.2, sex, Diff, `Std. Error`, df, `t value`, `Pr(>|t|)`, pvalue.adjusted )
    
    unadjusted.comps <- unadjusted.comps %>% 
      dplyr::select(strain.1, strain.2, sex,  Diff, `pval REML`) %>%
      rename('Diff Multi'='Diff')
  }else{
    unadjusted.comps <- NULL # 'All strains are similar'
  }
  
  if(all( Multi.lab.fit.summary$`Pr(>|t|)` == Multi.lab.fit.summary$pvalue.adjusted ) ){
    return(list(Unadjusted = unadjusted.comps))
  }
  
  # adjusted.comps <- Multi.lab.fit.summary %>% filter( pvalue.adjusted < 0.05 ) # %>%  mutate(strain.greater = strain.1, strain.smaller = strain.2 )
  if(nrow(adjusted.comps)>0){
    # adjusted.comps[adjusted.comps$Estimate < 0, c('strain.greater','strain.smaller')] <- 
    # adjusted.comps[adjusted.comps$Estimate < 0, c('strain.2','strain.1')]
  adjusted.comps <- adjusted.comps %>% 
    dplyr::select(strain.1, strain.2, sex, Estimate, pvalue.adjusted) # %>% mutate(Estimate = abs(Estimate) )
  }else{
    adjusted.comps <- NULL # 'All strains are similar after the adjustment for multiple comparisons'
  }
  
  return(list(Adjusted = adjusted.comps, Unadjusted = unadjusted.comps))

}



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
