library(tidyverse)
library(dplyr)

# setwd('~/Dropbox/A cross laboratory investigation of timing endophenotypes in mouse behavior/MPD pairs/full sheets June 19')
strains.formal.names <- sort( c('BALB/cJ' , 'BTBR', 'C3H/HeJ'  , 'C57BL/6J'  , 'CBA'  , 'DBA/2J'  , 'SWR/J'))

####################################################################################################
############################################ Wiltshire2 ############################################
####################################################################################################
dd <- read_csv('~/Library/Mobile Documents/com~apple~CloudDocs/Research proposal files/TAU_JAX_trial/MPD pairs/full sheets June 19/Wiltshire2.csv')
setwd('~/Library/Mobile Documents/com~apple~CloudDocs/Research proposal files/TAU_JAX_trial/MPD pairs/MPD data - arranged IJ')

dd$OFT_thig_cont <- as.numeric(dd$OFT_thig_cont)
# as.numeric(dd$OFT_thig_fluox) -> dd$OFT_thig_fluox
# as.numeric(dd$OFT_dist_cont) -> dd$OFT_dist_cont
# as.numeric(dd$OFT_dist_fluox) -> dd$OFT_dist_fluox
# as.numeric(dd$TST_fluox) -> dd$TST_fluox
# as.numeric(dd$TST_cont) -> dd$TST_cont

dd$OFT_thig_cont <- (100 - dd$OFT_thig_cont)/100
dd$OFT_thig_fluox<- (100 - dd$OFT_thig_fluox)/100

d1 <- dd[,c('strain' , 'OFT_thig_cont' , 'OFT_dist_cont','TST_cont')]
colnames(d1)[-1] <- c('OFT_thig' , 'OFT_dist','TST')
d1$group <- 'control'

d2 <- dd[,c('strain' , 'OFT_thig_fluox' , 'OFT_dist_fluox','TST_fluox')]
colnames(d2)[-1] <- c('OFT_thig' , 'OFT_dist','TST')
d2$group <- 'fluoxetine'

dd <- rbind(d1,d2) %>%
  mutate(treatment = group, 
         strain = recode(strain , 'BALB/cJ' = 'BALB/cJ', 'BTBR T+ tf/J' = 'BTBR', 'C3H/HeJ' = 'C3H/HeJ',
                         'C57BL/6J' = 'C57BL/6J', 'CBA/J' = 'CBA', 'DBA/2J' = 'DBA/2J', 'SWR/J' = 'SWR/J'))

write.csv(x = dd[,c('strain' , 'OFT_thig', 'group')]  %>% filter(!is.na(OFT_thig)) ,
          file = 'Lab_Wiltshire2_PercentCenter_Males.csv')

write.csv(x = dd[,c('strain' , 'OFT_dist', 'group')]  %>% filter(!is.na(OFT_dist)) ,
          file = 'Lab_Wiltshire2_DistanceTraveled_Males.csv')

write.csv(x = dd[,c('strain' , 'TST', 'group')]  %>% filter(!is.na(TST)) ,
          file = 'Lab_Wiltshire2_TST_Males.csv')

#####################################################################################################
#################################### Crowley2 body weight Males #####################################
#####################################################################################################

dd <- read.csv('mpd39401.csv' ,stringsAsFactors = F )
dd <- dd %>% filter(label == 'starting') %>% mutate(group = 'control')
dd <- dd[,c('strain' , 'value', 'group')] %>%
  mutate(strain = recode(strain, 
                         'BALB/cByJ' = 'BALB/cJ', 'BTBR T<+> Itpr3<tf>/J' = 'BTBR', 
                         'C3H/HeJ' = 'C3H/HeJ', 'C57BL/6J' = 'C57BL/6J',
                         'CBA/J' = 'CBA' , 'DBA/2J' = 'DBA/2J'))

colnames(dd)[2] <- 'MALES_body_weight'



write.csv(x = dd , file = 'Lab_Crowley2_body_weight_Males.csv')

#####################################################################################################
################ Tarantino2 females open field, distance traveled and percent center ################
#####################################################################################################

dd <- read.csv('mpd50601.csv' ,stringsAsFactors = F )
dd <- dd %>% filter(label == 'control') %>% mutate(group = 'control')
dd <- dd[,c('strain' , 'value', 'group')]
colnames(dd)[2] <- 'FEMALES_body_weight'
write.csv(x = dd , file = 'Lab_Tarantino2_DistanceTraveled_Females.csv')


dd <- read.csv('mpd50626.csv' ,stringsAsFactors = F )
dd <- dd %>% filter(label == 'control')  %>% mutate(group = 'control')

dd <- dd[,c('strain' , 'value', 'group')] %>%
  mutate(strain = recode(strain, 'BALB/cJ' = 'BALB/cByJ','BTBR' = 'BTBR T<+> Itpr3<tf>/J',
                         'C3H/HeJ' = 'C3H/HeJ','C57BL/6J' = 'C57BL/6J','CBA' = 'CBA/J',
                         'DBA/2J' = 'DBA/2J','SWR/J' = 'SWR/J'),
         value = value / 100 )
colnames(dd)[2] <- 'FEMALES_percentage_center_time'
write.csv(x = dd , file = 'Lab_Tarantino2_PercentCenter_Females.csv')


#####################################################################################################
############################## Tordoff3 body weight, Males and Females ##############################
#####################################################################################################

dd <- read.csv('mpd10305.csv' ,stringsAsFactors = F )
dd <- dd %>% filter(label == 'start')
dd <- dd %>% mutate( MALES_body_weight = replace( value , sex=='f' , NA ) ) %>%
  mutate( FEMALES_body_weight = replace( value , sex=='m' , NA ) , group = 'control')
dd <- dd[ , c('strain' ,'FEMALES_body_weight' , 'MALES_body_weight', 'group' ) ] %>%
  mutate(strain = recode(strain, 'BALB/cJ' = 'BALB/cByJ',
                         'BTBR' = 'BTBR T<+> Itpr3<tf>/J',
                         'C3H/HeJ' = 'C3H/HeJ', 'C57BL/6J' = 'C57BL/6J', 
                         'CBA' = 'CBA/J', 'DBA/2J' = 'DBA/2J', 'SWR/J' = 'SWR/J'))


write.csv(x = dd[ , c('strain' ,'FEMALES_body_weight', 'group' ) ] , file = 'Lab_Tordoff3_body_weight_Females.csv')
write.csv(x = dd[ , c('strain' ,'MALES_body_weight', 'group' ) ] , file = 'Lab_Tordoff3_body_weight_Males.csv')


#####################################################################################################
########################### Crabbe4 grip strength in males and in females ###########################
#####################################################################################################

dd <- read.csv('mpd19106.csv' ,stringsAsFactors = F )
dd <- dd %>% filter(label == 'baseline')
dd <- dd %>% mutate( MALES_grip_strength = replace( value , sex=='f' , NA ) ) %>%
  mutate( FEMALES_grip_strength = replace( value , sex=='m' , NA ), group = 'control')
dd <- dd[ , c('strain' ,'FEMALES_grip_strength' , 'MALES_grip_strength', 'group' ) ] %>%
  mutate(strain = recode( strain, 'BALB/cJ' = 'BALB/cByJ', 'BTBR' = 'BTBR T<+> Itpr3<tf>/J',
                          'C3H/HeJ' = 'C3H/HeJ', 'C57BL/6J' = 'C57BL/6J', 'CBA' = 'C58/J',
                          'DBA/2J' = 'DBA/2J', 'SWR/J' = 'SWR/J'))

# 
# strains.edit.names <- sort( unique(dd$strain))[c(4:7,9,11,21)]
# paste(strains.formal.names, strains.edit.names, sep = "' = '") %>% 
#   paste(., collapse = "', '")
# 

write.csv(x = dd[ , c('strain' ,'MALES_grip_strength', 'group' ) ] ,
          file = 'Lab_Crabbe4_grip_strength_Males.csv')

write.csv(x = dd[ , c('strain' ,'FEMALES_grip_strength', 'group' ) ] ,
          file = 'Lab_Crabbe4_grip_strength_Females.csv')


