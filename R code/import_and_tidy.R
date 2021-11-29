library('dplyr')
library("tidyverse")
library(readr)
library(readxl)

my.path <- '~/Library/Mobile Documents/com~apple~CloudDocs/Research proposal files/TAU_JAX_trial'

(files_names <- list.files(path = paste0(my.path, '/Data - Original files/'),
                           pattern="*.xlsx"))


for ( contemp_filename in files_names[8:12]){
  print( contemp_filename )
  cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
  print(paste(excel_sheets(path = cln ), collapse = "', '"))
  # contemp_xl <- read_excel( path = cln )
}

contemp_filename <- "GxL_JAX_TST_females.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
tst_FM_JAX <- read_excel( path = cln, sheet = "Data", skip = 3) 
contemp_filename <- "GxL_JAX_TST_males.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
tst_FM_JAX <- rbind( tst_FM_JAX, read_excel( path = cln, sheet = "Data", skip = 3) ) %>%
  group_by(strain, sex, treatment,id) %>%
  summarise('tst_6min_percent' = sum(s...17[Order<7])/(60*6),
            'tst_7min_percent' = sum(s...17[Order<8])/(60*7) ) 
tst_FM_JAX <- as_tibble(tst_FM_JAX) %>%
  filter(!is.na(strain)) %>%
  mutate( lab = 'JAX',
          strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",  "C3H/HeJ"= "C3H/HeJ",
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male')) %>%
  select( lab, strain, sex, treatment, tst_6min_percent, tst_7min_percent)


contemp_filename <- "GxL_JAX_OF_females.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
oft_FM_JAX <- read_excel( path = cln, sheet = "Data") 
contemp_filename <- "GxL_JAX_OF_males.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
oft_FM_JAX <- rbind( oft_FM_JAX, read_excel( path = cln, sheet = "Data") ) %>%
  rename('strain' = 'STRAIN', 'sex' = 'SEX', 'treatment' = 'TREATMENT',
         'sample' = 'SAMPLE', 'dist' = 'TOTAL DISTANCE (cm)',
         'centertime' = 'CENTER TIME LEGACY (s)') 

oft10_FM_JAX <- oft_FM_JAX %>%
  filter(!is.na(strain), sample < 3 ) %>%
  mutate( lab = 'JAX', strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",   "C3H/HeJ"= "C3H/HeJ",
                                        'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                                        'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                                        'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                                        'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male')) %>%
  group_by(strain, sex, treatment,ID) %>%
  summarise(OFTlarge_dist_10m_sec = sum(dist), # /(10*60),
            OFTlarge_centertime_10m_percent = sum(centertime)/(10*60),
            lab = 'JAX') %>%
  select( lab, strain, sex, treatment, 
          OFTlarge_centertime_10m_percent, 
          OFTlarge_dist_10m_sec)


oft20_FM_JAX <- oft_FM_JAX %>%
  filter(!is.na(strain)) %>%
  mutate( lab = 'JAX', strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",   "C3H/HeJ"= "C3H/HeJ",
                                        'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                                        'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                                        'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                                        'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male')) %>%
  group_by(strain, sex, treatment,ID) %>%
  summarise(OFTlarge_dist_20m_sec = sum(dist), #/(20*60)
            OFTlarge_centertime_20m_percent = sum(centertime)/(20*60),
            lab = 'JAX') %>%
  select( lab, strain, sex, treatment, 
          OFTlarge_centertime_20m_percent,
          OFTlarge_dist_20m_sec)


contemp_filename <- "GxL_JAX_Grip_males.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
grip_M_JAX <- read_excel( path = cln, sheet = "Data") %>%
  mutate(grip.avg = `average of forepaws`, lab = 'JAX' ) %>%
  select( lab, treatment, sex, strain, grip.avg) %>%
  filter(!is.na(strain)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",   "C3H/HeJ"= "C3H/HeJ",
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))

contemp_filename <- "GxL_JAX_Grip_females.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
grip_F_JAX <- read_excel( path = cln, sheet = "Data") %>%
  mutate(grip.avg = `average of forepaws`, lab = 'JAX' ) %>%
  select( lab, treatment, sex, strain, grip.avg) %>%
  filter(!is.na(strain)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",   "C3H/HeJ"= "C3H/HeJ",
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))

contemp_filename <- "GxL_JAX_BW_females.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
bw_F_JAX <- read_excel( path = cln, sheet = "Data") %>%  mutate(`Body Weight` = NA, lab = 'JAX' )
bw_F_JAX$`Body Weight` <- rowMeans(bw_F_JAX[,5:10], na.rm = T)
bw_F_JAX <- bw_F_JAX %>%
  filter(!is.na(strain)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",   "C3H/HeJ"= "C3H/HeJ",
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))

contemp_filename <- "GxL_JAX_BW_males.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
bw_M_JAX <- read_excel( path = cln, sheet = "Data") %>%  mutate(`Body Weight` = NA , lab = 'JAX')
bw_M_JAX$`Body Weight` <- rowMeans(bw_M_JAX[,5:10], na.rm = T)
bw_M_JAX <- bw_M_JAX %>%
  filter(!is.na(strain)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR', "BTBR T<+> Itpr3<tf>/J"="BTBR",  "C3H/HeJ"= "C3H/HeJ",
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))


contemp_filename <- "tail suspension for 6 minutes.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c("male orgenized", "female orgenized")
tst6_F_TAUM <- read_excel( path = cln, sheet = "female orgenized",skip = 3 )
tst6_M_TAUM <- read_excel( path = cln, sheet = "male orgenized", skip = 3)
names(tst6_F_TAUM) <-  c('cage', 'treatment', 'sex', 'strain', 'tst_6min_mean', 'tst_6min_sec')
names(tst6_M_TAUM)[1:6] <-  c('cage', 'treatment', 'sex', 'strain', 'tst_6min_mean', 'tst_6min_sec')
tst6_F_TAUM <- tst6_F_TAUM %>% filter(!is.na(strain)) %>%
  mutate( lab = 'TAUM', tst_6min_percent = tst_6min_sec/(60*6)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 
                           'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J', 'CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = 'Control',
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
tst6_M_TAUM <- tst6_M_TAUM[,1:6] %>% 
 mutate( lab = 'TAUM', 
         tst_6min_percent = tst_6min_sec/(60*6)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))


contemp_filename <- "tail suspension for 7 min (002).xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c("male orgenized", "female orgenized")
tst7_F_TAUM <- read_excel( path = cln, sheet = "female orgenized" ,skip = 3 )
tst7_M_TAUM <- read_excel( path = cln, sheet = "male orgenized" ,skip = 3 )
names(tst7_F_TAUM) <-  c('cage', 'treatment', 'sex', 'strain', 'tst_7min_sec')
names(tst7_M_TAUM) <-  c('cage', 'treatment', 'sex', 'strain', 'tst_7min_sec')
tst7_F_TAUM <- tst7_F_TAUM %>% filter(!is.na(strain)) %>%
mutate( lab = 'TAUM', treatment = 'Control', 
                                       tst_7min_percent = tst_7min_sec/(60*7)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J','DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
tst7_M_TAUM <- tst7_M_TAUM %>% filter(!is.na(strain)) %>%
mutate( lab = 'TAUM', tst_7min_percent = tst_7min_sec/(60*7)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))

contemp_filename <- "Grip test final (2).xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c('Sheet1', 'Sheet2', 'Sheet3', 'Sheet4', 'Sheet5', 'Sheet6', 'Sheet8', 'Sheet9', 'all orgenized')
grip_FM_TAUM <- read_excel( path = cln, sheet = 'all orgenized' , skip = 2) %>% mutate(grip.avg = NA) %>%
  filter( !is.na(`weight (gr)`))
grip_FM_TAUM$grip.avg <- rowMeans(grip_FM_TAUM[,c('trial 1 (gf)','trial 2 (gf)','trial 3 (gf)')], na.rm = T ) 
names(grip_FM_TAUM)[1:5] <- c('cage',  'strain', 'sex','treatment', 'Body Weight')
grip_FM_TAUM <- grip_FM_TAUM[,c(1:5,11)] %>%
  filter(!is.na(strain)) %>%
  mutate(lab = 'TAUM') %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                           'CBA ' = 'CBA/J', 'CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))

contemp_filename <- "Mice TAUL 04.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c('Sheet1', 'Sheet2', 'Sheet3')
all_FM_TAUL <- read_excel( path = cln, sheet = 'Sheet1' )
all_F_TAUL <- all_FM_TAUL %>% filter(sex=='F') %>% 
  rename( 'ID'='Mouse' , 'strain'='Strain', 'cage'='Cage ID', 'treatment' = 'Treatment',
          'grip.avg' = 'Grip AVG', 
          'OFTlarge_dist_10m_sec' = 'OFT1 Dist 10m', 
          'OFTlarge_dist_20m_sec' = 'OFT1 Dist 20m',
          'OFTlarge_centertime_10m_percent' ='OFT1 %C 10m' , 
          'OFTlarge_centertime_20m_percent' = 'OFT1 %C 20m',
          'tst_6min_percent'='TST 6m', 'tst_7min_percent'='TST m') %>%
  filter(!is.na(strain)) %>%
  mutate(OFTsmall_dist_10m_sec = NA, 
         OFTsmall_dist_20m_sec = NA ,
         OFTsmall_centertime_10m_percent = NA ,
         OFTsmall_centertime_20m_percent = NA ) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'),
          lab = 'TAUL')

all_M_TAUL <- all_FM_TAUL %>% filter(sex=='M') %>%
  rename( 'ID'='Mouse' , 'strain'='Strain', 'cage'='Cage ID', 'treatment' = 'Treatment',
          'grip.avg' = 'Grip AVG', 
          'OFTsmall_dist_10m_sec' = 'OFT1 Dist 10m', 
          'OFTsmall_centertime_10m_percent' ='OFT1 %C 10m' ,
          'OFTsmall_dist_20m_sec' = 'OFT1 Dist 20m',
          'OFTsmall_centertime_20m_percent' = 'OFT1 %C 20m',
          'OFTlarge_dist_10m_sec' = 'OFT2 Dist 10m', 
          'OFTlarge_centertime_10m_percent' ='OFT2 %C 10m' ,
          'OFTlarge_dist_20m_sec' = 'OFT2 Dist 20m',
          'OFTlarge_centertime_20m_percent' = 'OFT2 %C 20m',
          'tst_6min_percent'='TST 6m', 
          'tst_7min_percent'='TST m') %>%
  filter(!is.na(strain)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'),
          lab = 'TAUL')
all_M_TAUL <- all_M_TAUL %>% 
  filter(!is.na(strain)) %>%
  mutate(OFTsmall_dist_10m_sec = as.numeric(OFTsmall_dist_10m_sec),
            OFTsmall_centertime_10m_percent = as.numeric(OFTsmall_centertime_10m_percent),
            OFTsmall_dist_20m_sec = as.numeric(OFTsmall_dist_20m_sec),
            OFTsmall_centertime_20m_percent = as.numeric(OFTsmall_centertime_20m_percent),
            OFTlarge_dist_10m_sec = as.numeric(OFTlarge_dist_10m_sec),
            OFTlarge_centertime_10m_percent = as.numeric(OFTlarge_centertime_10m_percent),
            OFTlarge_dist_20m_sec = as.numeric(OFTlarge_dist_20m_sec),
            OFTlarge_centertime_20m_percent = as.numeric(OFTlarge_centertime_20m_percent), 
            tst_6min_percent = as.numeric(tst_6min_percent),
            tst_7min_percent = as.numeric(tst_7min_percent))
all_F_TAUL <- all_F_TAUL %>% 
  filter(!is.na(strain)) %>%
  mutate(OFTsmall_dist_10m_sec = as.numeric(OFTsmall_dist_10m_sec),
            OFTsmall_centertime_10m_percent = as.numeric(OFTsmall_centertime_10m_percent),
            OFTsmall_dist_20m_sec = as.numeric(OFTsmall_dist_20m_sec),
            OFTsmall_centertime_20m_percent = as.numeric(OFTsmall_centertime_20m_percent),
            OFTlarge_dist_10m_sec = as.numeric(OFTlarge_dist_10m_sec),
            OFTlarge_centertime_10m_percent = as.numeric(OFTlarge_centertime_10m_percent),
            OFTlarge_dist_20m_sec = as.numeric(OFTlarge_dist_20m_sec),
            OFTlarge_centertime_20m_percent = as.numeric(OFTlarge_centertime_20m_percent), 
            tst_6min_percent = as.numeric(tst_6min_percent),
            tst_7min_percent = as.numeric(tst_7min_percent))

temp <- (all_F_TAUL$OFTlarge_dist_20m_sec-all_F_TAUL$OFTlarge_dist_10m_sec)
hist(temp, breaks = 20, main = 'TAUL: Histogram of \n dist 20m - dist 10m')



contemp_filename <- "open field 10 min.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c('Female orgenized', '27X27 Male orgenized', '50X50 Male orgenized')
oft10_Flarge_TAUM <- read_excel( path = cln, sheet = 'Female orgenized',skip=3 )
oft10_Mlarge_TAUM <- read_excel( path = cln, sheet = '50X50 Male orgenized' ,skip=3 )
oft10_Msmall_TAUM <- read_excel( path = cln, sheet = '27X27 Male orgenized' ,skip=3 )
names(oft10_Flarge_TAUM) <- c('result','trial','arena', 'cage','sex','strain','OFTlarge_dist_10m_sec',
                              'velocity','activity','center-point Frequency','OFTlarge_centertime_10m_percent' )
names(oft10_Mlarge_TAUM) <- c('result','trial','arena', 'cage','treatment','sex','strain','OFTlarge_dist_10m_sec',
                              'velocity','activity','center-point Frequency','OFTlarge_centertime_10m_percent' )
names(oft10_Msmall_TAUM) <- c('result','trial','arena', 'cage','treatment','sex','strain','OFTsmall_dist_10m_sec',
                              'velocity','activity','center-point Frequency','OFTsmall_centertime_10m_percent' )
oft10_Flarge_TAUM <- oft10_Flarge_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( OFTlarge_centertime_10m_percent = OFTlarge_centertime_10m_percent / (10*60), 
                                                    treatmen='Control', lab = 'TAUM' ) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = 'Control',
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
oft10_Mlarge_TAUM <- oft10_Mlarge_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( OFTlarge_centertime_10m_percent = OFTlarge_centertime_10m_percent / (10*60), 
                                                    lab = 'TAUM') %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
oft10_Msmall_TAUM <- oft10_Msmall_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( OFTsmall_centertime_10m_percent = OFTsmall_centertime_10m_percent / (10*60), 
                                                    lab = 'TAUM') %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J', 'CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))


contemp_filename <- "open field 20 min NK.xlsx"
cln <- paste0(c(my.path, '/Data - Original files/', contemp_filename), collapse = '')
# cln_sheets <- c('Female Orgenized', 'Male 27X27 Orgenized', 'Male 50X50 Orgenized')
oft20_Flarge_TAUM <- read_excel( path = cln, sheet = 'Female Orgenized', skip = 3 )
oft20_Mlarge_TAUM <- read_excel( path = cln, sheet = 'Male 50X50 Orgenized', skip = 3 )
oft20_Msmall_TAUM <- read_excel( path = cln, sheet = 'Male 27X27 Orgenized', skip = 3 )
names(oft20_Flarge_TAUM)[c(1:7)] <-  c('trial','arena', 'cage', 'sex','strain','OFTlarge_dist_20m_sec','OFTlarge_centertime_20m_percent')
names(oft20_Mlarge_TAUM)[c(1:8)] <-  c('trial','arena', 'cage', 'treatment', 'sex','strain','OFTlarge_dist_20m_sec','OFTlarge_centertime_20m_percent')
names(oft20_Msmall_TAUM)[c(1:8)] <-  c('trial','arena', 'cage', 'treatment', 'sex','strain','OFTsmall_dist_20m_sec','OFTsmall_centertime_20m_percent')
oft20_Flarge_TAUM <- oft20_Flarge_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( lab = 'TAUM', treatment = 'Control', OFTlarge_centertime_20m_percent = OFTlarge_centertime_20m_percent / (20*60) ) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
oft20_Mlarge_TAUM <- oft20_Mlarge_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( lab = 'TAUM' , OFTlarge_centertime_20m_percent = OFTlarge_centertime_20m_percent/ (20*60)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J', 
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))
oft20_Msmall_TAUM <- oft20_Msmall_TAUM %>%
  filter(!is.na(strain)) %>%
  mutate( lab = 'TAUM' , OFTsmall_centertime_20m_percent = OFTsmall_centertime_20m_percent/ (20*60)) %>%
  mutate( strain = recode( na_if(strain, ''), .missing = 'NA', 'BTBR'='BTBR',
                           'BALB' = 'BALB/cJ','BALB ' = 'BALB/cJ','Balbc' = 'BALB/cJ', 'c57' = 'C57BL/6J',
                           'C57' = 'C57BL/6J','c57bl' = 'C57BL/6J','C57bl' = 'C57BL/6J',
                           'CBA ' = 'CBA/J','CBA' = 'CBA/J',
                           'DBA' = 'DBA/2J' ,'SWR ' = 'SWR/J', 'SWR'='SWR/J'),
          treatment = recode(na_if(treatment,''), .missing = 'Control',  'C' = 'Control','control' = 'Control',
                             'drug' = 'Fluoxetine', 'Fl' = 'Fluoxetine','fluoxetine' = 'Fluoxetine', 'water' = 'Control'),
          sex = recode(sex, 'f'='Female', 'm'='Male', 'F'='Female', 'M'='Male', 'female'='Female', 'male'='Male'))



################# Combine & save

################ OFT data

oft20_Msmall_TAUM <- oft20_Msmall_TAUM %>% 
  mutate(OFTlarge_dist_20m_sec = NA, OFTlarge_centertime_20m_percent = NA,
         OFTlarge_dist_10m_sec = NA, OFTlarge_centertime_10m_percent = NA,
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)


oft20_Mlarge_TAUM <- oft20_Mlarge_TAUM  %>%
  mutate(OFTlarge_dist_10m_sec = NA, OFTlarge_centertime_10m_percent = NA,
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA,
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)


oft20_Flarge_TAUM <- oft20_Flarge_TAUM  %>%
  mutate(OFTlarge_dist_10m_sec = NA, OFTlarge_centertime_10m_percent = NA,
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA,
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

oft10_Msmall_TAUM  <- oft10_Msmall_TAUM %>%
  mutate(OFTlarge_dist_20m_sec = NA, OFTlarge_centertime_20m_percent = NA,
         OFTlarge_dist_10m_sec = NA, OFTlarge_centertime_10m_percent = NA,
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

oft10_Mlarge_TAUM <- oft10_Mlarge_TAUM %>%
  mutate(OFTlarge_dist_20m_sec = NA, OFTlarge_centertime_20m_percent = NA,
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA,
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

oft10_Flarge_TAUM <- oft10_Flarge_TAUM %>%
  mutate(OFTlarge_dist_20m_sec = NA, OFTlarge_centertime_20m_percent = NA,
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA,
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA ) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

oft10_FM_JAX <- oft10_FM_JAX %>% 
  mutate(OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA,
         OFTlarge_dist_20m_sec = NA,  OFTlarge_centertime_20m_percent = NA, 
         OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

oft20_FM_JAX <- oft20_FM_JAX %>% 
  mutate(OFTsmall_dist_10m_sec = NA, OFTsmall_centertime_10m_percent = NA,
         OFTlarge_dist_10m_sec = NA,  OFTlarge_centertime_10m_percent = NA, 
         OFTsmall_dist_20m_sec = NA, OFTsmall_centertime_20m_percent = NA) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

OFT_F_TAUL <-  all_F_TAUL %>%
  mutate(OFTlarge_centertime_20m_percent = OFTlarge_centertime_20m_percent,
         OFTsmall_centertime_20m_percent = OFTsmall_centertime_20m_percent,
         OFTlarge_centertime_10m_percent = OFTlarge_centertime_10m_percent,
         OFTsmall_centertime_10m_percent = OFTsmall_centertime_10m_percent) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

OFT_M_TAUL <- all_M_TAUL %>%
  mutate(OFTlarge_centertime_20m_percent = OFTlarge_centertime_20m_percent,
         OFTsmall_centertime_20m_percent = OFTsmall_centertime_20m_percent,
         OFTlarge_centertime_10m_percent = OFTlarge_centertime_10m_percent,
         OFTsmall_centertime_10m_percent = OFTsmall_centertime_10m_percent) %>%
  select(lab, treatment, sex, strain, 
         OFTlarge_dist_20m_sec, OFTlarge_centertime_20m_percent,
         OFTlarge_dist_10m_sec, OFTlarge_centertime_10m_percent,
         OFTsmall_dist_10m_sec, OFTsmall_centertime_10m_percent,
         OFTsmall_dist_20m_sec, OFTsmall_centertime_20m_percent)

OFT_combined_data <- bind_rows(oft20_Msmall_TAUM, oft20_Mlarge_TAUM, OFT_F_TAUL, OFT_M_TAUL, oft10_FM_JAX, oft20_FM_JAX, oft10_Flarge_TAUM, oft10_Mlarge_TAUM, oft10_Msmall_TAUM, oft20_Flarge_TAUM ) %>%
  as_tibble(.)
OFT_combined_data <- OFT_combined_data %>% pivot_longer(cols = -c(lab, treatment, sex, strain))  %>% filter(!is.na(value))  %>%
  filter( strain !='C3H/HeJ') %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))

ggplot(OFT_combined_data) +
  geom_boxplot(aes(x=name, y = value, color = treatment))+
  facet_grid(sex~lab)+
  theme(axis.text = element_text(angle=90))+
  scale_y_log10()
write.csv(OFT_combined_data, file = 'OFT_combined_data.csv')


OFT_combined_data %>% 
  filter(name %in% c('OFTsmall_dist_10m_sec','OFTsmall_dist_20m_sec','OFTlarge_dist_10m_sec','OFTlarge_dist_20m_sec'),
         lab=='TAUL', treatment=='Control') %>%
  ggplot(.) +
  geom_histogram(aes(value))+
  facet_grid(name~sex)+
  theme(axis.text = element_text(angle=90))+
  ggtitle('Distance Travelled, TAUL')

################ TST data
tst7_M_TAUM <- tst7_M_TAUM %>% mutate(tst_6min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst7_F_TAUM <- tst7_F_TAUM %>% mutate(tst_6min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst6_M_TAUM <- tst6_M_TAUM %>% mutate(tst_7min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst6_F_TAUM <- tst6_F_TAUM %>% mutate(tst_7min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst_FM_JAX <-tst_FM_JAX %>%  
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst_F_TAUL <- all_F_TAUL %>% # mutate(tst_6min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)
tst_M_TAUL <- all_M_TAUL %>% # mutate(tst_6min_percent = NA ) %>%
  select(lab, treatment, sex, strain, tst_6min_percent, tst_7min_percent)

TST_combined_data <- bind_rows(tst_M_TAUL, tst_F_TAUL, tst_FM_JAX, tst6_F_TAUM, tst6_M_TAUM, tst7_F_TAUM, tst7_M_TAUM ) %>%
  as_tibble(.) %>%
  filter(!is.na(strain),!is.na(sex))  %>%
  filter( strain !='C3H/HeJ') %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))

write.csv(TST_combined_data, file = 'TST_combined_data.csv')

################ bodyweight data
bw_M_JAX <- bw_M_JAX %>% select(lab, treatment, sex, strain, `Body Weight`)
bw_F_JAX <- bw_F_JAX %>% select(lab, treatment, sex, strain, `Body Weight`)
bw_F_TAUL <- all_F_TAUL %>% select(lab, treatment, sex, strain, `Body Weight`)
bw_M_TAUL <- all_M_TAUL %>% select(lab, treatment, sex, strain, `Body Weight`)
bw_FM_TAUM <- grip_FM_TAUM %>% select(lab, treatment, sex, strain, `Body Weight`)

bw_combined_data <- bind_rows(bw_FM_TAUM, bw_M_TAUL, bw_F_TAUL, bw_F_JAX, bw_M_JAX) %>% as_tibble(.)  %>%
  filter( strain !='C3H/HeJ') %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))

write.csv(bw_combined_data, file = 'bw_combined_data.csv')

################ grip data
grip_FM_TAUM <- grip_FM_TAUM %>% select(lab, treatment, sex, strain, grip.avg)
grip_F_JAX <- grip_F_JAX %>% select(lab, treatment, sex, strain, grip.avg)
grip_M_JAX <- grip_M_JAX %>% select(lab, treatment, sex, strain, grip.avg)
grip_F_TAUL <- all_F_TAUL %>% select(lab, treatment, sex, strain, grip.avg)
grip_M_TAUL <- all_M_TAUL %>% select(lab, treatment, sex, strain, grip.avg)

grip_combined_data <- bind_rows(grip_FM_TAUM, grip_M_TAUL, grip_F_TAUL, grip_F_JAX, grip_M_JAX) %>% as_tibble(.)  %>%
  filter( strain !='C3H/HeJ',
          !is.na(grip.avg)) %>%
  mutate(strain = factor(strain, levels = c("BALB/cJ","BTBR","C57BL/6J","DBA/2J","SWR/J", "CBA/J",' ', '  ')))

write.csv(grip_combined_data, file = 'grip_combined_data.csv')
