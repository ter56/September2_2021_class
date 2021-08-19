library(arm)
library(lme4)
library(Hmisc)
library(plyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GAPIT3)
library(tidyverse)
library(readxl)
library(purrrlyr)

# Let's start with the example we saw in the slides 5 loci, equal effect sizes.
small_loci_num = 5
P1_sl = rep(2,5)
P2_sl = rep(0,5)

numOffspring_sl = 200
numLoci_sl = 5

offspring_sl = data.frame()
for (i in 1:numOffspring_sl){
  for (j in 1:numLoci_sl){
    temp = sum(temp,sample(c(P1_sl[j],P2_sl[j]),1))
  }
  offspring_sl = rbind(offspring_sl, data.frame (ProgNum = i, value = temp))
  temp = 0
}
offspring_sl %>% ggplot(aes(x =value)) +geom_histogram(bins = 6)+labs(title = '5 loci segregating, equal effect')

# Now lets see what happens with effect sizes that differ
# lets create a vector of marker allele substitution effects
MarkerEffects_sl = c(1,6,2,1,1)
# Same thing as above - play with changing the effect sizes and see what happens
numOffspring_sl = 300
offspring_sl_Ef = data.frame()
numLoci_sl = 5

for (i in 1:numOffspring_sl){
  for (j in 1:numLoci_sl){
    temp = sum(temp,sample(c(P1_sl[j],P2_sl[j]),1)*MarkerEffects_sl[j])
  }
  offspring_sl_Ef = rbind(offspring_sl_Ef, data.frame (ProgNum = i, value = temp))
  temp = 0
}
offspring_sl_Ef %>% ggplot(aes(x =value)) +geom_histogram(bins = 12)+labs(title = '5 loci segregating, different marker effects')

# So now lets extend things: 
# A vector of 2's and 0's (fully inbred) is the same as a vector of 1's and 0's and is a little easier to think about
# Lets say now we have 200 loci, with equal effect sizes, but now there are some sites the parents have in common 
# (ie the sites are not segregating), and P2 no loger has only the negative allele. 
numLoci = 200
# 
# P1 = rep(1,numLoci)
# P2 = rep(0,numLoci)

P1 = sample(c(0,1), numLoci, replace = T,prob = c(0.5,0.5))
P2 = sample(c(0,1), numLoci, replace = T,prob = c(0.5,0.5))

numOffspring = 500
offspring = data.frame() 
temp=0

# mutate(x1 = sample(P1,P2))

for (i in 1:numOffspring){
  for (j in 1:numLoci){
    temp = sum(temp,sample(c(P1[j],P2[j]),1))
  }
  offspring = rbind(offspring, data.frame (ProgNum = i, value = temp))
  temp = 0
}
offspring %>% ggplot(aes(x = value))+geom_density()+geom_vline(xintercept = c(sum(P1),sum(P2))) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = c(sum(P1),sum(P2)), label = c('P1','P2')),
             aes(x =x, y = y, label = label))

P1s = c(sample(c(0,1), numLoci, replace = T), 20)
P2s = c(sample(c(0,1), numLoci, replace = T), 0)

offsprings = data.frame() 
temp=0
for (i in 1:numOffspring){
  for (j in 1:numLoci+1){
    temp = sum(temp,sample(c(P1s[j],P2s[j]),1))
  }
  offsprings = rbind(offsprings, data.frame (ProgNum = i, value = temp))
  temp = 0
}
sinlgeLargeEffect = offsprings %>% ggplot(aes(x = value))+geom_density()+geom_vline(xintercept = c(sum(P1s),sum(P2s))) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = c(sum(P1s),sum(P2s)), label = c('P1','P2')),
             aes(x =x, y = y, label = label))


P1ranef = sample(c(0,1), numLoci, replace = T)
P2ranef = sample(c(0,1), numLoci, replace = T)
MarkerEf = rnorm(numLoci, mean = 1)

offspringranef = data.frame() 
temp=0
for (i in 1:numOffspring){
  for (j in 1:numLoci){
    temp = sum(temp,sum((sample(c(P1s[j],P2s[j]),1)))*MarkerEf[j])
  }
  offspringranef = rbind(offspringranef, data.frame (ProgNum = i, value = temp))
  temp = 0
}
ranefPop = offspringranef %>% ggplot(aes(x = value))+geom_density()+geom_vline(xintercept = c(sum(P1ranef*MarkerEf),sum(P2ranef*MarkerEf))) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = c(sum(P1ranef*MarkerEf),sum(P2ranef*MarkerEf)), label = c('P1','P2')),
             aes(x =x, y = y, label = label))

Norm+sinlgeLargeEffect+ranefPop +plot_layout(ncol = 1)



#########################################################################################
setwd('WMB_Data/')
DHs = rbind(read_excel("DHtp1_all.xlsx")%>%mutate(TP = 'TP1',PM_date =5),
            read_excel("DHtp2_all.xlsx")%>%mutate(TP = 'TP2',PM_date =19),
            read_excel("DHtp3_all.xlsx")%>%mutate(TP = 'TP3',PM_date =47),
            read_excel("DHtp4_all.xlsx")%>%mutate(TP = 'TP4',PM_date =96),
            read_excel("DHtp5_all.xlsx")%>%mutate(TP = 'TP5',PM_date =152)) %>%
  mutate(GE =(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
         GI = 10*GE*(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+2*Day2Germ+3*Day3Germ),
         GI = ifelse(is.nan(GI),0,GI),
         GE5 = (Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
         GI5 = 10*GE5*(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ*2+Day3Germ*3+Day4Germ*4+Day5Germ*5),
         GE4 =(Day1Germ+Day2Germ+Day3Germ+Day4Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
         GI4 = GE4 *10*(Day1Germ+Day2Germ+Day3Germ+Day4Germ)/(Day1Germ+Day2Germ*2+Day3Germ*3+Day4Germ*4),
         Location = ifelse(substr(PLOT,1,2)=='11','Snyder','Caldwell'),
         rep = as.factor(replication),
         taxa = mapvalues(Entry, from = c('Check 1','Check 2','Check 3','Check 4','Check 5','Check 6'),
                          to = c('Flavia', 'Scala','DH130910','SY Tepee','Wintmalt','Charles')),
         taxa =  mapvalues(taxa, from = c("Check 1-Flavia","Check 2-Scala","Check 3-DH130910","Check 3-SY Tepee",
                                          "Check 4-SY Tepee","Check 5-Wintmalt","Check 6-Charles"),
                           to = c('Flavia', 'Scala','DH130910','DH130910','SY Tepee','Wintmalt','Charles')))
setwd(rprojroot::find_rstudio_root_file())
load('SMB_Data/PHS_BLUEs_GGS1920.RData')

Model_get_blups_H2 = function(trait.lm, trait){
  # This function fits trait~Year+(1|taxa)+Location+rep using lmer()
  # the `trait` is a string names of the trait measured ('betaGlucan' for example)
  # This function will output the cullis H2, as the data is unbalanced
  # and will return the model blups, (randon effect of entry plus the fixed effect of intercept)
  model.blups=as.data.frame(ranef(trait.lm))[,3:4] %>% mutate(condval =condval+as.numeric(fixef(trait.lm)[1]),
                                                              trait = trait) %>%
    rename(taxa = grp, BLUP = condval )
  # H2 calculation according to Cullis et al 2006
  ses<- se.ranef(trait.lm)$'taxa' #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(trait.lm, comp="Variance")$'taxa'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated
  H2<- round(mean(Reliability),3) #This is equivalent to broad-sense heritability on the line-mean (or family-mean, if your individuals are non-inbred families) basis
  print(H2)
  return(list(BLUP = model.blups, H2 = H2))
}
get_BLUES = function(dataframe, trait, FirstTaxa){
  temp.lm = broom::tidy(lm(dataframe[[trait]] ~ taxa +Location+rep, data = dataframe))
  Intercept = as.numeric(temp.lm[which(temp.lm$term=='(Intercept)'),2])
  temp.blues = temp.lm %>% filter(substr(term,1,4) =='taxa') %>% 
    add_row(term = FirstTaxa, estimate = 0) %>% # add back in the first entry as it drops it and puts effect as the intercept. 
    mutate(BLUE = estimate+Intercept, trait = trait,
           term = gsub(pattern = 'taxa', replacement = '', x = term)) %>%
    rename(taxa = term)
  return(temp.blues)
}

DH_blups_H2 = data.frame(trait= NA, H2 = NA)
DH_blups = data.frame()
TP_list = c('TP1','TP2','TP3','TP4','TP5')
FirstTaxaList = c('BS611-2','BS611-1','BS611-1','BS611-1','BS611-1')
DH_blues = data.frame()
counter = 1
for (i in TP_list) {
  print(i)
  Model = Model_get_blups_H2(lmer(GI ~ (1|taxa)+Location+rep,
                                  DHs %>% filter(PLOT %nin%c('SN','SR') & TP ==i)),
                             trait = paste0('GI_',i))
  DH_blups = rbind(DH_blups,Model$BLUP)
  DH_blups_H2 = DH_blups_H2 %>% add_row(trait = paste0('GI_',i), H2 = Model$H2)
  
  ModelGE = Model_get_blups_H2(lmer(GE ~ (1|taxa) +Location+rep,
                                    DHs %>% filter(PLOT %nin%c('SN','SR') & TP ==i)),
                               trait = paste0('GE_',i))
  DH_blups = rbind(DH_blups,ModelGE$BLUP)
  DH_blups_H2 = DH_blups_H2 %>% add_row(trait = paste0('GE_',i), H2 = ModelGE$H2)
  
  TempGE = get_BLUES(DHs %>% filter(PLOT %nin%c('SN','SR') & TP ==i),
                     trait = 'GE',FirstTaxa = FirstTaxaList[counter]) %>%
    mutate(trait = paste0(trait,'_',i))
  TempGI = get_BLUES(DHs %>% filter(PLOT %nin%c('SN','SR') & TP ==i),
                     trait = 'GI',FirstTaxa = FirstTaxaList[counter])%>%
    mutate(trait = paste0(trait,'_',i))
  DH_blues  = rbind(DH_blues, TempGE, TempGI)
  counter = counter + 1
}

DH_blues = DH_blues %>% separate(remove = F, col = trait,into = c('Trait','TP')) %>% 
  arrange(TP, Trait, taxa)%>%
  mutate(Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9',
                                                       'DH1','Fla','SY ','Sca','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910',
                                   'SY Tepee/DH130910','Wintmalt/DH130910',
                                   'Parent','Parent','Parent','Parent','Parent')))
DH_blups = DH_blups %>% separate(remove = F, col = trait,into = c('Trait','TP'))%>% 
  arrange(TP, Trait, taxa)%>%
  mutate(Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9',
                                                       'DH1','Fla','SY ','Sca','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910',
                                   'SY Tepee/DH130910','DH130910/Wintmalt',
                                   'Parent','Parent','Parent','Parent','Parent')))


DH_blues %>% filter(TP == 'TP5' & Trait == 'GI' & Family !='Parent') %>% filter(taxa != 'Charles') %>% 
  ggplot()+geom_density(aes(x = BLUE)) +facet_wrap(facets = vars(Family)) +
  geom_vline(data = DH_blues %>% 
               filter(TP == 'TP5' & Trait == 'GI' & Family =='Parent')%>% 
               mutate(Family = paste0(taxa,'/DH130910')) %>% filter(taxa != 'DH130910') %>%
               rbind(
                 DH_blues %>% 
                   filter(TP == 'TP5' & Trait == 'GI' & taxa == 'DH130910') %>%
                   select(!Family) %>% mutate(count = 4) %>% uncount(count) %>%
                   mutate(Family = c('Flavia/DH130910','Scala/DH130910',
                                     'SY Tepee/DH130910','Wintmalt/DH130910'))),
             aes(xintercept = BLUE)
  )+
  geom_text(data = DH_blues %>% 
              filter(TP == 'TP5' & Trait == 'GI' & Family =='Parent')%>% 
              mutate(Family = paste0(taxa,'/DH130910')) %>% filter(taxa != 'DH130910') %>%
              rbind(
                DH_blues %>% 
                  filter(TP == 'TP5' & Trait == 'GI' & taxa == 'DH130910') %>%
                  select(!Family) %>% mutate(count = 4) %>% uncount(count) %>%
                  mutate(Family = c('Flavia/DH130910','Scala/DH130910',
                                    'SY Tepee/DH130910','Wintmalt/DH130910'))),
            aes(x = BLUE, y = .5, label = taxa)
  ) +
  xlab('Germination Rate at TP5') +ylab('Density') + theme_bw()+labs(title = 'Germination Rate in winter\nmalting barley biparental crosses')

PHS.blues %>% ggplot(aes(x = PHS)) +geom_density() +theme_bw() +labs(title = 'PHS in barley Nam-style Cross ')


hieght = read_excel('SMB_Data/GGS_20 fieldbook_DS.xlsx',sheet = 'Helfer')
hieght %>% ggplot(aes(x = Ht)) +geom_density() +labs(title = 'Plant height in barley Nam-style Cross')+theme_bw() 

sample(c('Wr','S'),size = 500, replace=T ) %>% data.frame() %>%
  ggplot(aes( x = .))+geom_bar() +labs(title = 'Mendels Peas')

load('SMB_Data/GGS2020_BLUE_summary_allTP.RData')
all_BLUE %>% select(!c(GE5,GI5scale)) %>% rename(GI3 = GI3scale) %>%
  pivot_longer(cols = c(GE3, GI3))%>% mutate(facet = paste0(name,TP))%>%
  ggplot(aes(value))+facet_wrap(facets = vars(facet), scales = 'free', nrow = 2) +geom_density() +
  theme_bw()
