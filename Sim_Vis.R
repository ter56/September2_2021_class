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
library(rrBLUP)

# Let's start with the example we saw in the slides: 5 loci, equal effect sizes.
P1_sl = rep(2,5) #P1 haplotype always carrying the 2
P2_sl = rep(0,5) #P2 haplotype always carrying the 0

numOffspring_sl = 200 #number of offspring to generate
numLoci_sl = 5  # our 5 loci!
Parents_sl = matrix(c(P1_sl,P2_sl), nrow = 5, ncol = 2)
offspring_sl = matrix(NA, nrow = 5, ncol = 200)
for (row in 1:nrow(Parents_sl)){
  ithLocus = sample(c(Parents_sl[row,]),200, replace = T)
  offspring_sl[row,] = ithLocus
}

markerEffects_SL = rep(1,5)
offspring_sl_BV = data.frame(BV = t(t(markerEffects_SL)%*%offspring_sl ))
Parents_sl_BV = data.frame(BV = t(t(markerEffects_SL)%*%Parents_sl ))
# If you'd like to look at what this is doing as a for loop here it is essentially:
# temp = 0
# offspring_sl_BV = data.frame()
# for (i in 1:numOffspring_sl){
#   for (j in 1:numLoci_sl){
#     temp = sum(temp,sample(c(P1_sl[j],P2_sl[j]),1))
#   }
#   offspring_sl = rbind(offspring_sl, data.frame (ProgNum = i, value = temp))
#   temp = 0
# }
# Plot the results
offspring_sl_BV %>% ggplot(aes(x =BV)) +geom_histogram(bins = 6)+labs(title = '5 loci segregating, equal effect')
# This is the same as a binomial distribution with n = 5, k = number positive loci, p=q=0.5

# Now lets see what happens with effect sizes that differ
# lets create a vector of marker allele substitution effects
DifMarkerEffects_sl = c(1,60,2,1,1)
# Same thing as above - play with changing the effect sizes and see what happens
offspring_SLD_BV = data.frame(BV = t(t(DifMarkerEffects_sl)%*%offspring_sl ))
Parents_SLD_BV = data.frame(BV = t(t(DifMarkerEffects_sl)%*%Parents_sl ))
# Plot the results
offspring_SLD_BV %>% ggplot(aes(x =BV)) +geom_histogram(bins = 6)+labs(title = '5 loci segregating, differnt effect size')

###### MORE MARKERS! ######
# So now lets extend things: 
# A vector of 2's and 0's (fully inbred) is the same as a vector of 1's and 0's (just divide by 2) and is a little easier to think about.
# Lets say now we have 200 loci, but now there are some sites the parents have in common 
# (ie the sites are not segregating), and P2 no loger has only the negative allele. 
numLoci = 200
# Lets create parents with more loci
P1 = c(sample(c(0,1), numLoci-1, replace = T,prob = c(0.5,0.5)),1)
P2 = c(sample(c(0,1), numLoci-1, replace = T,prob = c(0.5,0.5)),0)

#put those parents into a matrix
ParentMatrix = matrix(c(P1,P2), nrow = numLoci, ncol = 2)
numOffspring = 500 #decide on the number of offspinrg to simulate
offspring = matrix(NA, nrow = numLoci, ncol = numOffspring) #and make a matrix with that size.
# We will populate the offspirng matrix with the offspring each locus at a time (as things are unlinked) using this 'for' loop:
for (row in 1:nrow(ParentMatrix)){
  ithLocus = sample(c(ParentMatrix[row,]),numOffspring, replace = T)
  offspring[row,] = ithLocus
}
# which makes a matrix with the offsping in columns, and the loci in rows!
# now let's use this matrix to examine offsping values: 
# Let's assume equal effect sizes, 1 being positive and 0 being the negative. 
# That is for each offsping we are going to find the sum of the marker effects, assign that 
# as the individuals breeding value and then plot them!
# All loci have equal effects, so lets assing a vector of '1's as the effect size (this is that vector of a's that we saw on the slide)
Equaleffect_markers = rep(1,numLoci)
Parent_BV = data.frame(BV = t(Equaleffect_markers %*% ParentMatrix)) #
Offspring_BV =  data.frame(BV = t(Equaleffect_markers %*% offspring))
Offspring_BV%>% ggplot(aes(x = BV))+geom_density()+geom_vline(xintercept = Parent_BV$BV) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = Parent_BV$BV, label = c('P1','P2')),
             aes(x =x, y = y, label = label))

# Now lets do the same thing, but include a Single Large Effect (SLE) locus at the end of the parnetal haplotypes that 
# the parents are segregating for - see above where P1 and P2 are inititalized
SLE_markerEffect = c(rep(1,numLoci-1),20)
Parent_SLE_BV = data.frame(BV = t(t(SLE_markerEffect) %*% ParentMatrix))
Offspring_SLE_BV = data.frame(BV = t(t(SLE_markerEffect) %*% offspring))

Offspring_SLE_BV %>% ggplot(aes(x = BV))+geom_density()+geom_vline(xintercept = Parent_SLE_BV$BV) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = Parent_SLE_BV$BV, label = c('P1','P2')),
             aes(x =x, y = y, label = label))




# Now what happens when we assign a random effect size to all the markers (ie marker effects are not constant) 
# in this the effects are drawn from a normal distribution with mean and sd of 1. 
# The parental haplotype now inditcates presence/absence!
Random_MarkerEffect = rnorm(numLoci, mean = 1, sd= 1)
Parent_Randef_BV = data.frame(BV = t(t(Random_MarkerEffect) %*% ParentMatrix))
Offspring_Randef_BV = data.frame(BV = t(t(Random_MarkerEffect) %*% offspring))

Offspring_Randef_BV %>% ggplot(aes(x = BV))+geom_density()+geom_vline(xintercept = Parent_Randef_BV$BV) +
  geom_label(data = data.frame(y = c(0.05,0.05), x = Parent_Randef_BV$BV, label = c('P1','P2')),
             aes(x =x, y = y, label = label))


###  Back to the presenation! #####

x = sample(c(1,0),100,replace = T)
y = rnorm(mean = x*6, sd=1, n = 100)
data.frame(x = x,y=y) %>% ggplot(aes(x =x, y =y))+geom_point() +xlab('Marker Score')+
  ylab('Breeding Value') +geom_smooth(formula = y~x, se = F, method = lm)




############## Let use the random vs fixed effects models and see how they respond
# Lets also write a function to give us the BV, marker matrix, and etc:
GenerateBiparental = function(NumLoci, PopulationSize, MarkerEffects){
  ParentMatrix = matrix(sample(c(0,1),NumLoci*2,replace = T), nrow = NumLoci, ncol = 2)
  Offspring = matrix(NA, nrow = NumLoci, ncol = PopulationSize) 
  ExtraOffspring = matrix(NA, nrow = NumLoci, ncol = PopulationSize)#and make a matrix with that size.
  # We will populate the offspirng matrix with the offspring each locus at a time (as things are unlinked) using this 'for' loop:
  for (row in 1:nrow(ParentMatrix)){
    ithLocus = sample(c(ParentMatrix[row,]),PopulationSize, replace = T)
    Offspring[row,] = ithLocus
    ithLocus2 = sample(c(ParentMatrix[row,]),PopulationSize, replace = T)
    ExtraOffspring[row,] = ithLocus2
  }
  Parent_BV = data.frame(BV = t(MarkerEffects %*% ParentMatrix)) #
  Offspring_BV =  data.frame(BV = t(MarkerEffects %*% Offspring))
  ExtraOffspring_BV =  data.frame(BV = t(MarkerEffects %*% ExtraOffspring))
  
  Plot = Offspring_BV %>% ggplot(aes(x = BV))+geom_density()+geom_vline(xintercept = Parent_BV$BV) +
    geom_label(data = data.frame(y = c(0.05,0.05), x = Parent_BV$BV, label = c('P1','P2')),
               aes(x =x, y = y, label = label))
  
  return(list(ParentMatrix =ParentMatrix,
              Parent_BV = Parent_BV,
              Offspring_BV = Offspring_BV,
              MarkerEffects = MarkerEffects,
              Offspring = Offspring,
              Plot = Plot,
              ExtraOffspring = ExtraOffspring,
              ExtraOffspring_BV = ExtraOffspring_BV))
  
}


library(rrBLUP)
# 5 loci with numerous individuals works well
t(offspring_sl) %>% data.frame() %>% cbind(., offspring_sl_BV$BV) %>%
  rename(BV = 'offspring_sl_BV$BV') %>%
  lm(BV ~ X1+X2+X3+X4+X5, data = .)

# now 200 markers with 500 individuals - equal effect markers remember
t(offspring) %>% data.frame() %>% cbind(., Offspring_BV$BV) %>%
  rename(BV = 'Offspring_BV$BV') %>%
  lm(BV ~ . , data = .)
# works well everything had an effect size of 1

# How about about normally distributed markers?
randeffect.lm = t(offspring) %>% data.frame() %>% cbind(., Offspring_Randef_BV$BV) %>%
  rename(BV = 'Offspring_Randef_BV$BV') %>%
  lm(BV ~ . , data = .)
# lets check correltion of marker effects.
Random_MarkerEffect %>% cbind(data.frame(randeffect.lm$coefficients)[-1,]) %>% cor(., use = 'complete.obs')
# That works well to - mono morphic markers are dropped from the equation, so really we only have about 100 markers 
# that we are estimating. ie N > P by 5x

#now lets get more realistic
P4000_N500 = GenerateBiparental(NumLoci = 1000,PopulationSize = 400, MarkerEffects = rnorm(1000,mean = 0, sd = 1))
P4000_N500$Plot
P4k_N500.lm = t(P4000_N500$Offspring) %>% data.frame() %>%cbind(., P4000_N500$Offspring_BV) %>%
  lm(BV~.,data = .)
P4000_N500$MarkerEffects %>% cbind(data.frame(P4k_N500.lm$coefficients)[-1,]) %>% cor(., use = 'complete.obs')
# Low or No correlation between estimated effects and actual known effect!
MarkerCoef4k_n500= as.data.frame(P4k_N500.lm$coefficients) #let get out the coef
MarkerCoef4k_n500[is.na(MarkerCoef4k_n500)] <- 0 #replace nas with 0s so that things can be calcuated  - notice there are 499 of them 
cor(P4000_N500$ExtraOffspring_BV$BV,#TrueValues
    t(matrix(MarkerCoef4k_n500$`P4k_N500.lm$coefficients`)[-1,]  %*% P4000_N500$ExtraOffspring))#What we calcuate
#NO Correlation! everything is over fit and underwhelming!

# Now what if we used a ridge regression process to regularize and make the predictors less sensitive?
# rrblup is a good framework to work through 
Ueffects = mixed.solve(y = P4000_N500$Offspring_BV$BV, Z = t(P4000_N500$Offspring))
#Now lets check to see if these marker effects along with the extra offspring simulated and there true breeding values are correlated
cor(t(as.matrix(P4000_N500$ExtraOffspring)) %*% as.matrix(Ueffects$u),P4000_N500$ExtraOffspring_BV$BV)
 
# Lets look at the Markers now and the estimated vs actual effects:
data.frame(Estimated = as.matrix(Ueffects$u), Known =P4000_N500$MarkerEffects) %>%
  ggplot(aes(x = Known, y = Estimated))+geom_point()+ #Markers that are monomorphic get estimate effect sizes of zero hence those at zero
  geom_abline(slope = 1, intercept = 0)







####################### Actual GP #########
setwd('WMB_Data/')
load('WinterGD.RData')
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
WinterGD[1:4,1:4]
dim(WinterGD)
# This data is set up with all this location and so need to be summarized into sinlge trait values per line

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


## Basis:
Phenotype = DH_blues %>% filter(TP == 'TP2' & Trait == 'GI')%>% arrange(taxa) %>% filter(taxa %in% WinterGD$taxa)
myGD = WinterGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa) 
dim(Phenotype);dim(myGD)
if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
  print('taxa lists are not correctly aligning')
}
# 80/20 training and testing split
TrueandPredicted = data.frame()
datasplit = 0.8
Training_size = round(datasplit*dim(Phenotype)[1],0)

set.seed(1)
Phenotype = as.vector(Phenotype$BLUE) #correct format for RRblup
num_entries = length(Phenotype)
myGD = myGD[,-1]-1 #correct format for RRblup
trainingSet = sort(sample(1:num_entries,Training_size))
testingSet = setdiff(1:num_entries,trainingSet)

y_train = Phenotype[trainingSet] 
y_test = Phenotype[testingSet]

marker_train = myGD[trainingSet,]
marker_test = myGD[testingSet,]

trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)

PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)
print(cor(y_test,PredictedPheno))


# Now wrap it into a function
FoldCVGPv2 = function(Phenotype, myGD, numfolds,datasplit, trait_name){
  # Phenotype is full vector of phenotypes as a dataframe with column names of 'taxa' 
  # and 'BLUE'
  # myGD is a n taxa x N markers dataframe with 0,1,2 coding and a 'taxa' column as the first column in it
  # numfolds is the number of folds to cv (ie 5 for five-fold cv)
  # datasplit is the percentage as a decimal for the training/testing split.
  # Trait name is a string of the desired name output.
  
  Phenotype = Phenotype %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  
  TrueandPredicted = data.frame()
  
  Training_size = round(datasplit*num_entries,0)
  start_time <- Sys.time() 
  set.seed(1)
  Phenotype = as.vector(Phenotype$BLUE)
  num_entries = length(Phenotype)
  myGD = myGD[,-1]-1
  for (i in 1:numfolds){
    
    trainingSet = sort(sample(1:num_entries,Training_size))
    testingSet = setdiff(1:num_entries,trainingSet)
    
    y_train = Phenotype[trainingSet] 
    y_test = Phenotype[testingSet]
    
    marker_train = myGD[trainingSet,]
    marker_test = myGD[testingSet,]
    
    trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
    
    PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+as.numeric(trained.Model$beta)
    print(cor(y_test,PredictedPheno))
    
    TrueandPredicted = rbind(TrueandPredicted, 
                             data.frame(TruePheno = y_test,
                                        PredPheno = PredictedPheno,
                                        fold = i,
                                        trait = trait_name,
                                        correlation = cor(y_test,PredictedPheno)))
  }
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}

TP2_GI_Prediction = 
  FoldCVGPv2(Phenotype = DH_blues %>% filter(TP == 'TP2' & Trait == 'GI')%>% arrange(taxa),
             myGD = WinterGD,
             numfolds = 5,
             datasplit = .8,
             trait_name = 'TP2_GI')

TP2_GI_Prediction %>% ggplot(aes(x = TruePheno, y = PredPheno))+
  geom_point()+geom_abline(slope = 1, intercept = 0, color = 'red')+
  geom_smooth(method = 'lm')

####Association analysis example using GAPIT #############
load('SMB_Data/PHS_BLUEs_GGS1920.RData')
load(file = 'SMB_Data/myGD_LDpruned_w_KASP.RData')
load(file = 'SMB_Data/myGM_LDpruned_w_KASP.RData')

DF_OperationsV2 = function(dataframe){
  dataframe = dataframe[order(dataframe$P.value),]
  dataframe$logPercentileQQplot = -log10(c(1:length(dataframe$SNP))/length(dataframe$SNP))
  dataframe$rank = c(1:length(dataframe$SNP))
  dataframe$FDRPval = length(dataframe$SNP)*dataframe$P.value/dataframe$rank
  dataframe = dataframe[order(dataframe$Chromosome, as.numeric(dataframe$Position)),]
  dataframe$log10PVal = -log10(dataframe$P.value)
  dataframe$ordinal = c(1:length(dataframe$SNP))
  return(dataframe)
} #Operations to be performed on GAPIToutput$GWAS so that homeMadeManhattanLDPruned works
TableOutput = function(GWAS.sum.dataframe, n = 20){
  GWAS.sum.dataframe[order(GWAS.sum.dataframe$rank),c('SNP','Chromosome','Position','P.value','maf','FDRPval')][1:n,] %>% 
    mutate(maf = as.numeric(as.character(maf))) %>%
    kable(digits = c(0,1,0,8,2,6), align = 'lcccrr')}


PHS.MLM.gwas = PHS.blues %>% GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                                   Geno.View.output=F, model="MLM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
PHS.MLMM.gwas = PHS.blues %>% GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                                   Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)

SNPperChr = table(PHS.MLM.gwas$GWAS$Chromosome)
MiddleOfChr = SNPperChr/2
breaks = c(); breaks[1]= MiddleOfChr[1];
ChromLines = c();ChromLines[1] = SNPperChr[1]
for(i in 2:length(SNPperChr)){
  breaks[i] = sum(SNPperChr[1:i-1])+MiddleOfChr[i]
  ChromLines[i] = sum(SNPperChr[1:i])
}


DF_OperationsV2(PHS.MLM.gwas$GWAS) %>% filter(maf >0.07) %>%
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'All MTA, MAF>0.07') 

DF_OperationsV2(PHS.MLMM.gwas$GWAS) %>% filter(maf >0.07) %>%
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  theme_bw()+labs(title = 'All MTA, MAF>0.07')
  

##################Random Plots made for the presentation #######################################################################

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


# This is a test

#for louis 

# For David 

