# Script for making tables and Figures (except figures 7 and S9) in:
  
# Kardan, O., Stier, A. J., Cardenas-Iniguez, C., Schertz, K. E., Pruin, J. C., Deng, Y.,
# Chamberlain, T., Meredith, W. J.,  Zhang, X., Bowman, J. E., Lakhtakia, T., Tindel, L.,
# Avery, E. W., Yoo, K., Lin, Q., Chun, M. M., Berman, M. G., & Rosenberg, M. D.
# (accepted manuscript). Differences in the functional brain architecture of sustained attention
# and working memory in youth and adults. PLOS Biology 

# Data for this paper are at NDA (https://nda.nih.gov/) Study 1849 DOI: 10.15154/1528288

# questions? email Omid Kardan omidk@med.umich.edu

setwd("abcd/n-back_block_level_conn_analysis/")

library(tidyverse)
library(RColorBrewer)
library(stargazer)
library(lme4)
library(sjPlot)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(psych)
library(ppcor)

df.block <- read_csv("ABCD_nbkRun1and2_blocks_beh_brain_data.csv")
df.block_FDfilter <- filter(df.block, FDav < .2, FDmax < 2, nbkqc_eprm ==0 ,laptop_nback ==1, boxswitch ==0) #exclude runs with high motion etc.

# supplementary less strict about motion
# change above filter to FDav < .5, FDmax < 5
#

df.block_FDfilter$Sex <- as.factor(df.block_FDfilter$Sex)
df.block_FDfilter$N <- as.factor(df.block_FDfilter$N)
df.block_FDfilter$Type <- as.factor(df.block_FDfilter$Type)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)
df.block_FDfilter$Site <- as.factor(df.block_FDfilter$Site)

summary(df.block_FDfilter$N) # NaN

df.block_FDfilter <- filter(df.block_FDfilter, N != "NaN", blockFD != "NaN") %>% droplevels() # 
(SitesN <- df.block_FDfilter %>% distinct(subs, .keep_all = TRUE) %>% group_by(Site) %>% count() )
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
length(unique(df.block_FDfilter$subs)) #1545 subs and 19216 blokcs

#### race/ethnicity representations 
ethnicity <- read_csv('/acspsw03.csv')  # download file acspsw03 from NDA before running this line
gg<-df.block_FDfilter %>% group_by(subs) %>% summarise(subjectkey = first(subs))
eths <- merge(ethnicity,gg,by = 'subjectkey')
#eths$race_ethnicity <- as.factor(eths$race_ethnicity)
hist(eths$race_ethnicity)
sum(eths$race_ethnicity==1,na.rm = TRUE)/1545 # white
sum(eths$race_ethnicity==2,na.rm = TRUE)/1545 # black
sum(eths$race_ethnicity==3,na.rm = TRUE)/1545 # hispanic
sum(eths$race_ethnicity==4,na.rm = TRUE)/1545 # asian
sum(eths$race_ethnicity==5,na.rm = TRUE)/1545 # other
#####

#### sex and age
df.block_FDfilter %>% group_by(subs) %>% summarise(Sex = first(Sex)) %>% group_by(Sex) %>% count()
# Male = 694 Female = 851
ages<-df.block_FDfilter %>% group_by(subs) %>% summarise(Age = mean(Age))

mean(ages$Age)/12 # 10.03 years
####

filter(df.block_FDfilter, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores

df.block_FDfilter$Type <- factor(df.block_FDfilter$Type, levels = c("NeutFace", "NegFace", "PosFace", "Place"))
df.block_FDfilter$Sex <- factor(df.block_FDfilter$Sex,  levels = c("M", "F"))
df.block_FDfilter$Sex<-recode(df.block_FDfilter$Sex, `F` = "Female", M = "Male")
df.block_FDfilter$blockOrdered <- factor(df.block_FDfilter$block, ordered = T)

## for corrected r values in the supplementary analysis non-participation weights
nps <- read_csv("/ABCD_nonparticipationweights_HPG2.csv")

##########################################################################
#################BETWEEN-SUBJECT DATA AND PLOTS#######################
######################################################################
df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp) , brainscNIHcpm = mean(brainscNIHcpm))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN") 

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

##########
# mixed effects regressions across subs 

lmsa <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory + (1|Site), data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory + (1|Site), data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","0-back Accuracy"))
gg<-pcor(df.Sublevel.back2[,c("AccZ","WorkingMemory","Age","Motion")])
gg$estimate
gg$p.value


lmsa <- lmer(AccZ ~ Age+Sex+ Motion +SustainedAttention + (1|Site), data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lmer(AccZ ~ Age+Sex+ Motion +SustainedAttention + (1|Site), data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","0-back Accuracy"))
gg<-pcor(df.Sublevel.back0[,c("AccZ","SustainedAttention","Age","Motion")])
gg$estimate
gg$p.value


lmsa <- lmer(AccZ ~ Age+Sex+ Motion +CognitiveCompositeCPM + (1|Site), data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lmer(AccZ ~ Age+Sex+ Motion +CognitiveCompositeCPM + (1|Site), data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsb,lmsa, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("0-back Accuracy","2-back Accuracy"))
gg<-pcor(df.Sublevel.back2[,c("AccZ","CognitiveCompositeCPM","Age","Motion")])
gg$estimate
gg$p.value

##### test difference between betas
#Function to get difference between bootstrapped coefficients cc and wm
boot_func <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledCC <- lm(AccZ ~ Age + Sex + Motion + CognitiveCompositeCPM, data=resampled_data)
  lmsa_resampledWM <- lm(AccZ ~ Age + Sex + Motion + WorkingMemory, data=resampled_data)
  
  unname(lmsa_resampledCC$coefficients["CognitiveCompositeCPM"]) - unname(lmsa_resampledWM$coefficients["WorkingMemory"])
}

#Function to get difference between bootstrapped coefficients cc and SA
boot_func2 <- function(data){
  resampled_data = data %>% 
    sample_frac(.5, replace=TRUE)
  
  lmsa_resampledCC <- lm(AccZ ~ Age + Sex + Motion + CognitiveCompositeCPM, data=resampled_data)
  lmsa_resampledWM <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data)
  
  unname(lmsa_resampledCC$coefficients["CognitiveCompositeCPM"]) - unname(lmsa_resampledWM$coefficients["SustainedAttention"])
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0))
#histogram of these permuted correlations
hist(boot_coefs)
(pval = 1- sum(boot_coefs>0)/length(boot_coefs))

#################### permutations and np corrected r values 
cor.test(df.Sublevel.back2$AccZ, df.Sublevel.back2$WorkingMemory, method = "spearman" )
cor.test(df.Sublevel.back2$AccZ, df.Sublevel.back2$WorkingMemory, method = "pearson" )
# corr with non-participation weights
np_df.Sublevel.back2 <- merge(df.Sublevel.back2,nps,by.y = 'src_subject_id',by.x= 'subs')
xi <- np_df.Sublevel.back2$WorkingMemory
yi <- np_df.Sublevel.back2$AccZ
cor.test(xi, yi )
wi <- np_df.Sublevel.back2$wt_NR_HPGXacs
muwx = sum(wi*xi)/sum(wi)
muwy = sum(wi*yi)/sum(wi)
sig2wx = sum(wi*(xi-muwx)^2)/sum(wi)
sig2wy = sum(wi*(xi-muwy)^2)/sum(wi) 
covw = sum(wi*(xi-muwx)*(yi-muwy))/sum(wi)
r_corrected = covw/(sqrt(sig2wx)*sqrt(sig2wy))

#Smaller dataframe just to keep what we need
perm_Sub2back <- df.Sublevel.back2 %>% select(Acc, brainscRosenbergSA)
#Function to get permuted correlation (could be written more generally but right now just for these two columns and this one data frame)
perm_cors_func <- function(data){
  #Getting a shuffled acc column
  perm_Sub2back <- transform( perm_Sub2back, Acc2 = sample(Acc) )
  #Getting the correlation between WM and new shuffled column
  this_cor <- cor(perm_Sub2back$brainscRosenbergSA, perm_Sub2back$Acc2)
  #returning this correlation as the output of the function
  this_cor
}
#replicating the function 1000 times and saving the output as a vector
perm_cors <- replicate(1000, perm_cors_func(perm_Sub2back))
#histogram of these permuted correlations
hist(perm_cors)

############## scatter plots
p1 <- ggplot(data=df.Sublevel.back2) + 
  aes(x = WorkingMemory, y = Acc*100) + 
  geom_point(alpha=.15, size = 2,color = "#00AFBB") +
  geom_smooth(method = "lm", aes(x = WorkingMemory, y = Acc*100), size=1, color = "#00AFBB") +
  #  annotate("text", label = "r = .135, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = "2-back Accuracy (%)", x = " ", color="N-back", title = " ") 
p1
Fig4_tab4 <- data.frame(df.Sublevel.back2$WorkingMemory,df.Sublevel.back2$Acc*100)
write_csv(Fig4_tab4,"/figure files/Fig4_bottomR.csv")

cor.test(df.Sublevel.back2$Acc, df.Sublevel.back2$SustainedAttention )

p2 <- ggplot(data=df.Sublevel.back2) + 
  aes(x = SustainedAttention, y = Acc*100) + 
  geom_point(alpha=.15, size = 2, color = "#E7B800") +
  geom_smooth(method = "lm", aes(x = SustainedAttention, y = Acc*100), size=1, color = "#E7B800") +
  #  annotate("text", label = "r = -.040, N.S.", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = " ", x = " ", color="N-back", title = " ") 
p2
Fig4_tab3 <- data.frame(df.Sublevel.back2$SustainedAttention,df.Sublevel.back2$Acc*100)
write_csv(Fig4_tab3,"/figure files/Fig4_bottomL.csv")

#ggsave("SA_2back.png")
cor.test(df.Sublevel.back2$Acc, df.Sublevel.back2$CognitiveCompositeCPM)

np_df.Sublevel.back2 <- merge(df.Sublevel.back2,nps,by.y = 'src_subject_id',by.x= 'subs')
xi <- np_df.Sublevel.back2$CognitiveCompositeCPM
yi <- np_df.Sublevel.back2$Acc
cor.test(xi, yi )
wi <- np_df.Sublevel.back2$wt_n_NR_HPG
muwx = sum(wi*xi)/sum(wi)
muwy = sum(wi*yi)/sum(wi)
sig2wx = sum(wi*(xi-muwx)^2)/sum(wi)
sig2wy = sum(wi*(xi-muwy)^2)/sum(wi) 
covw = sum(wi*(xi-muwx)*(yi-muwy))/sum(wi)
r_corrected = covw/(sqrt(sig2wx)*sqrt(sig2wy))

p3 <- ggplot(data=df.Sublevel.back2) + 
  aes(x = CognitiveCompositeCPM, y = Acc*100) + 
  geom_point(alpha=.15, size = 2,color = "#FC4E07") +
  geom_smooth(method = "lm", aes(x = CognitiveCompositeCPM, y = Acc*100), size=1, color = "#FC4E07") +
  #  annotate("text", label = "r = .316, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = " ", x = " ", color="N-back", title = " ") 
p3
Fig5_tab2 <- data.frame(df.Sublevel.back2$CognitiveCompositeCPM,df.Sublevel.back2$Acc*100)
write_csv(Fig5_tab2,"/figure files/Fig5_bottomR.csv")


cor.test(df.Sublevel.back0$Acc, df.Sublevel.back0$WorkingMemory)

np_df.Sublevel.back0 <- merge(df.Sublevel.back0,nps,by.y = 'src_subject_id',by.x= 'subs')
xi <- np_df.Sublevel.back0$WorkingMemory
yi <- np_df.Sublevel.back0$AccZ
cor.test(xi, yi )
wi <- np_df.Sublevel.back0$wt_NR_HPGXacs
muwx = sum(wi*xi)/sum(wi)
muwy = sum(wi*yi)/sum(wi)
sig2wx = sum(wi*(xi-muwx)^2)/sum(wi)
sig2wy = sum(wi*(xi-muwy)^2)/sum(wi) 
covw = sum(wi*(xi-muwx)*(yi-muwy))/sum(wi)
r_corrected = covw/(sqrt(sig2wx)*sqrt(sig2wy))

p4 <- ggplot(data=df.Sublevel.back0) + 
  aes(x = WorkingMemory, y = Acc*100) + 
  geom_point(alpha=.15, size = 2,color = "#00AFBB") +
  geom_smooth(method = "lm", aes(x = WorkingMemory, y = Acc*100), size=1, color = "#00AFBB") +
  #  annotate("text", label = "r = .152, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = "0-back Accuracy (%)", x = "WM Network Strength", color="N-back", title = " ") 
p4
Fig4_tab2 <- data.frame(df.Sublevel.back0$WorkingMemory,df.Sublevel.back0$Acc*100)
write_csv(Fig4_tab2,"/figure files/Fig4_topR.csv")


cor.test(df.Sublevel.back0$AccZ, df.Sublevel.back0$SustainedAttention)

np_df.Sublevel.back0 <- merge(df.Sublevel.back0,nps,by.y = 'src_subject_id',by.x= 'subs')
xi <- np_df.Sublevel.back0$SustainedAttention
yi <- np_df.Sublevel.back0$AccZ
cor.test(xi, yi )
wi <- np_df.Sublevel.back0$wt_NR_HPGXacs
muwx = sum(wi*xi)/sum(wi)
muwy = sum(wi*yi)/sum(wi)
sig2wx = sum(wi*(xi-muwx)^2)/sum(wi)
sig2wy = sum(wi*(xi-muwy)^2)/sum(wi) 
covw = sum(wi*(xi-muwx)*(yi-muwy))/sum(wi)
r_corrected = covw/(sqrt(sig2wx)*sqrt(sig2wy))

p5 <- ggplot(data=df.Sublevel.back0) + 
  aes(x = SustainedAttention, y = Acc*100) + 
  geom_point(alpha=.25, size = 2, color = "#E7B800") +
  geom_smooth(method = "lm", aes(x = SustainedAttention, y = Acc*100), size=1, color = "#E7B800") +
  #  annotate("text", label = "r = .193, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = " ", x = "SA Network Strength", color="N-back", title = " ") 
p5
Fig4_tab1 <- data.frame(df.Sublevel.back0$SustainedAttention,df.Sublevel.back0$Acc*100)
write_csv(Fig4_tab1,"/figure files/Fig4_topL.csv")
#ggsave("SA_0back.png")

cor.test(df.Sublevel.back0$Acc, df.Sublevel.back0$CognitiveCompositeCPM)
p6 <- ggplot(data=df.Sublevel.back0) + 
  aes(x = CognitiveCompositeCPM, y = Acc*100) + 
  geom_point(alpha=.15,size = 2, color = "#FC4E07") +
  geom_smooth(method = "lm", aes(x = CognitiveCompositeCPM, y = Acc*100), size=1, color = "#FC4E07") +
  #  annotate("text", label = "r = .225, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = " ", x = "Cog. Comp. Network Strength", color="N-back", title = " ") 
p6
Fig5_tab1 <- data.frame(df.Sublevel.back0$CognitiveCompositeCPM,df.Sublevel.back0$Acc*100)
write_csv(Fig5_tab1,"/figure files/Fig5_topR.csv")

ggarrange(p1,p2,p3,p4,p5,p6,ncol = 3,nrow=2)
#ggsave("Six_scatterplots.png")


######### composite of nih toolbox performance
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
df.back02s <- df.block_FDfilter
df.back02s <- filter(df.back02s,  !is.na(avnih5),Acc !="NaN", avnih5 != "NaN") 

df.Sublevel.back02 <- df.back02s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), 
            FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp),brainscNIHcpm = mean(brainscNIHcpm)) %>% ungroup()

df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(composite = (avnih5-mean(avnih5, na.rm=T))/sd(avnih5, na.rm=T))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))


p7 <- ggplot(data=df.Sublevel.back02) + 
  aes(x = CognitiveCompositeCPM, y = composite) + 
  geom_point(alpha=.15,size = 2, color = "#FC4E07") +
  geom_smooth(method = "lm", aes(x = CognitiveCompositeCPM, y = composite), size=1, color = "#FC4E07") +
  annotate("text", label = "r = .295, p < .001", x = 0, y = 3, size=7) +
  theme_minimal(base_size = 20) + coord_cartesian(xlim = c(-5,5), ylim = c(-4,4))+
  labs(y = "Cognitive Composite Score ", x = "Cognitive Composite Network Strength", color="N-back", title = " ") 
p7
FigS4_dat <- data.frame(df.Sublevel.back02$CognitiveCompositeCPM,df.Sublevel.back02$composite)
write_csv(FigS4_dat,"/figure files/FigS4_dat.csv")

# site-wise cog comp cpm performances
#SitesN <- as.character(unique(df.Sublevel.back02$Site))
SitesN <- c(2:16,18,20,21)
cor.test(df.Sublevel.back02$avnih5, df.Sublevel.back02$CognitiveCompositeCPM, method = "spearman")
cor.test(df.Sublevel.back0$avnih5, df.Sublevel.back0$AccZ, method = "spearman")
cor.test(df.Sublevel.back2$avnih5, df.Sublevel.back2$AccZ, method = "spearman")

library(magrittr)

siteCor <- function(data, site){
  siteCor <- data %>% filter(Site == site) %$% cor(avnih5, CognitiveCompositeCPM)
  siteCor
}
cors = c(1:18)
for(i in 1:18){
  j = SitesN[i]
  cors[i] <- siteCor(df.Sublevel.back02, j)
}

sitecors <- as.data.frame(cors)


sitePerm <- function(data, site){
  sitedat <- data %>% filter(Site == site) 
  temp <- transform(sitedat, nih5 = sample(avnih5) )
  null_cor <- cor(temp$nih5, temp$CognitiveCompositeCPM)
  null_cor
}
nullcorgen<- function(){
  null_cors = c(rep(0,times=18))
  for(i in 1:18){
    j = SitesN[i]
    null_cors[i] <- sitePerm(df.Sublevel.back02, j)
  }
  null_cors
}

perm_cors <- replicate(1000, nullcorgen())

col_names = c("S02","S03","S04","S05","S06","S07","S08","S09", "S10", "S11","S12","S13","S14","S15","S16","S18","S20","S21")
perm_cors.df <- as.data.frame(t(perm_cors))
colnames(perm_cors.df) <- col_names
perm_cors.df.long <- pivot_longer(perm_cors.df, everything())

sitecors$name <- col_names
write_csv(sitecors,"/figure files/FigS5_model.csv")
sitecounts <-df.Sublevel.back02 %>% group_by(Site) %>% count()
sitecounts$sites <- col_names
write_csv(perm_cors.df.long,"/figure files/FigS5_null.csv")
p05upper <- perm_cors.df.long %>% group_by(name) %>% summarize(Q=quantile(value,.95))

p14 <- ggplot(data=perm_cors.df.long) + aes(x=name, y = value) + geom_violin(fill = 'gray', color='gray') +
  stat_mean(data=sitecors, inherit.aes = FALSE, aes(x = name, y =cors), color = "blue",shape=19, size=2)+
  #  geom_text(data=sitecounts, aes(x = sites, y = -.5, label = n))+
  stat_mean(data=p05upper, inherit.aes = FALSE, aes(x = name, y =Q), color = "red", shape="-", size=10)+
  theme_minimal(base_size = 20)+
  theme(axis.ticks.x = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
  labs(x = "Sites ", y = "Correlation value in each site", title = " ") 
p14

################

#William's and Stiger's tests
library(psych)
r.test(1545, r12 = .135, r34 = NULL, r23 = -.11, r13 = -.040, r14 = NULL, r24 = NULL, 
       n2 = NULL,pooled=TRUE, twotailed = TRUE)

r.test(1545, r12 = .193, r34 = NULL, r23 = .16, r13 = .152, r14 = NULL, r24 = NULL, 
       n2 = NULL,pooled=TRUE, twotailed = TRUE)

r.test(1545, r12 = .193, r34 = NULL, r23 = .62, r13 = -.04, r14 = NULL, r24 = NULL, 
       n2 = NULL,pooled=TRUE, twotailed = TRUE)

r.test(1545, r12 = .135, r34 = NULL, r23 = .62, r13 = .152, r14 = NULL, r24 = NULL, 
       n2 = NULL,pooled=TRUE, twotailed = TRUE)
# corplots

togetherforCor <- merge(df.Sublevel.back0[,c("subs","WorkingMemory","SustainedAttention", "CognitiveCompositeCPM")], df.Sublevel.back2[,c("subs","WorkingMemory","SustainedAttention", "CognitiveCompositeCPM")], by = "subs")
togetherforCor <- togetherforCor %>% rename(WorkingMemory_0back = WorkingMemory.x,
                                            WorkingMemory_2back = WorkingMemory.y,
                                            SustainedAttention_0back = SustainedAttention.x,
                                            SustainedAttention_2back = SustainedAttention.y,
                                            CognitiveCompositeCPM_0back = CognitiveCompositeCPM.x,
                                            CognitiveCompositeCPM_2back = CognitiveCompositeCPM.y) %>% dplyr::select(-subs)
(corM_both <- cor(togetherforCor[,c(1,4,2,5,3,6)]))
corrplot(corM_both, method="color", type="lower", tl.pos = "n", diag=FALSE)
corrplot(corM_both, method="number",number.cex = 1.5, type="lower", tl.col = "black",tl.pos = "n", col="black", bg=NULL, cl.pos = "n", diag=FALSE, tl.srt = 45, add=TRUE)

ggsave("corrmat_nolabel.png")

corrplot(corM_both, method="color", type="lower", tl.pos = "n")
corrplot(corM_both, method="number",number.cex = 2, type="lower", tl.col = "black", col="black", bg=NULL, tl.srt = 45, add=TRUE)

# scatterplot matrix
togetherforCor <- merge(df.Sublevel.back0[,c("subs","WorkingMemory","SustainedAttention")], df.Sublevel.back2[,c("subs","WorkingMemory","SustainedAttention")], by = "subs")
togetherforCor <- togetherforCor %>% rename(WorkingMemory_0back = WorkingMemory.x,
                                            WorkingMemory_2back = WorkingMemory.y,
                                            SustainedAttention_0back = SustainedAttention.x,
                                            SustainedAttention_2back = SustainedAttention.y,
) %>% dplyr::select(-subs)
pairs(togetherforCor)
#write_csv(togetherforCor,'togetherforCor.csv')   ## -> data for figure 3


#######################################RecMem
# collapse across both 0-back and 2-back blocks to test rec mem across the entire n-back
df.back02s <- df.block_FDfilter
df.back02s <- filter(df.back02s, Acc !="NaN") 

df.Sublevel.back02 <- df.back02s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIHcpm = mean(brainscNIHcpm) ,brainscAveryTask = mean(brainscAveryTask),
            FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
lms1 <- lmer(recmem ~ Age+Sex+ Motion +AccZ+WorkingMemory+SustainedAttention +(1|Site), data=df.Sublevel.back02)
summary(lms1)
tab_model(lms1,show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("Recognition Memory d'"))

# scatterplots of recmem
library(ppcor)
df.Sublevel.back02 <- filter(df.Sublevel.back02, recmem !="NaN") 
cor.test(df.Sublevel.back02$recmem, df.Sublevel.back02$WorkingMemory, method = "spearman")

np_df.Sublevel.back02 <- merge(df.Sublevel.back02,nps,by.y = 'src_subject_id',by.x= 'subs')
xi <- np_df.Sublevel.back02$recmem
yi <- np_df.Sublevel.back02$WorkingMemory
cor.test(xi, yi )
wi <- np_df.Sublevel.back02$wt_NR_HPGXacs
muwx = sum(wi*xi)/sum(wi)
muwy = sum(wi*yi)/sum(wi)
sig2wx = sum(wi*(xi-muwx)^2)/sum(wi)
sig2wy = sum(wi*(xi-muwy)^2)/sum(wi) 
covw = sum(wi*(xi-muwx)*(yi-muwy))/sum(wi)
r_corrected = covw/(sqrt(sig2wx)*sqrt(sig2wy))


gg<-pcor(df.Sublevel.back02[,c("recmem","WorkingMemory","Motion")])
gg$estimate
gg$p.value
p1 <- ggplot(data=df.Sublevel.back02) + 
  aes(x = WorkingMemory, y = recmem) + 
  geom_point(alpha=.15,size =4, color = "#00AFBB") +
  geom_smooth(method = "lm", aes(x = WorkingMemory, y = recmem), size=1, color = "#00AFBB") +
  #  annotate("text", label = "r = .123, p < .001", x = 0, y = -1, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5))+
  labs(y = "Recognition Memory d'", x = "Working Memory Network Strength", title = "") 
p1
FigS3_tab2 <- data.frame(df.Sublevel.back02$WorkingMemory,df.Sublevel.back02$recmem)
write_csv(FigS3_tab2,"/figure files/FigS3_right.csv")


cor.test(df.Sublevel.back02$recmem, df.Sublevel.back02$SustainedAttention, method = "spearman")
gg<-pcor(df.Sublevel.back02[,c("recmem","SustainedAttention","Motion")])
gg$estimate
gg$p.value
p2 <- ggplot(data=df.Sublevel.back02) + 
  aes(x = SustainedAttention, y = recmem) + 
  geom_point(alpha=.15, size =4,color = "#E7B800") +
  geom_smooth(method = "lm", aes(x = SustainedAttention, y = recmem), size=1, color = "#E7B800") +
  # annotate("text", label = "r = .007, N.S.", x = 0, y = -1, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5))+
  labs(y = "Recognition Memory d'", x = "Sustained Attention Network Strength",  title = "") 
p2
FigS3_tab1 <- data.frame(df.Sublevel.back02$SustainedAttention,df.Sublevel.back02$recmem)
write_csv(FigS3_tab1,"/figure files/FigS3_left.csv")



cor.test(df.Sublevel.back02$recmem, df.Sublevel.back02$CognitiveCompositeCPM)


p3 <- ggplot(data=df.Sublevel.back02) + 
  aes(x = CognitiveCompositeCPM, y = recmem) + 
  geom_point(alpha=.15, size =2,color = "#FC4E07") +
  geom_smooth(method = "lm", aes(x = CognitiveCompositeCPM, y = recmem), size=1, color = "#FC4E07") +
  annotate("text", label = "r = .241, p < .001", x = 0, y = -1, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5))+
  labs(y = "Recognition Memory d'", x = "Cognitive Composite Network Strength", title = "") 
p3

ggarrange(p1,p2,p3,ncol = 1,nrow=3)
ggsave("recmem_3scatterplots.png")


############################################WITHIN-SUBJECT#####################################

# z-scoring variables so that beta coefficients become comparable
df.back0 <- filter(df.block_FDfilter, N == '0')
df.back0 <- filter(df.back0,Acc !="NaN") 
df.back0 <- df.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.back0 <- df.back0 %>% mutate(CognitiveComposite = (brainscNIH-mean(brainscNIH))/sd(brainscNIH))
df.back0 <- df.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.back0 <- df.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.back0 <- df.back0 %>% mutate(BlockMotion = (blockFD-mean(blockFD))/sd(blockFD))
df.back0 <- df.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.back0 <- df.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
length(unique(df.back0$subs))


df.back2 <- filter(df.block_FDfilter, N == '2')
df.back2 <- filter(df.back2,Acc !="NaN") 
df.back2 <- df.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.back2 <- df.back2 %>% mutate(CognitiveComposite = (brainscNIH-mean(brainscNIH))/sd(brainscNIH))
df.back2 <- df.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.back2 <- df.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.back2 <- df.back2 %>% mutate(BlockMotion = (blockFD-mean(blockFD))/sd(blockFD))
df.back2 <- df.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.back2 <- df.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


###### regressions within-sub predicting acc

lm11 <- lmer(AccZ ~  BlockMotion+Type +Run+CognitiveCompositeCPM+(1|subs), data=df.back0)
lm12 <- lmer(AccZ ~  BlockMotion+Type +Run+CognitiveCompositeCPM+(1|subs), data=df.back2)
tab_model(lm11,lm12,  show.stat = TRUE, string.stat = "t-Statistic",
          show.se = FALSE, string.se = "SE",dv.labels =c( "0-back Accuracy","2-back Accuracy"),digits = 2, digits.p = 3)


lm11 <- lmer(AccZ ~  BlockMotion+Type+Run+WorkingMemory+SustainedAttention +(1|subs), data=df.back0)
lm12 <- lmer(AccZ ~  BlockMotion+Type+Run+WorkingMemory+SustainedAttention +(1|subs), data=df.back2)
tab_model(lm11,lm12,  show.stat = TRUE, string.stat = "t-Statistic",
          show.se = FALSE, string.se = "SE",dv.labels =c( "0-back Accuracy","2-back Accuracy"),digits = 2, digits.p = 3)


#mean-centering within subjects to remove individual differences for supplementary mediation analyses with stim types
df.back0 <- df.back0 %>% group_by(subs) %>% 
  mutate(ACCz_MC = (AccZ-mean(AccZ))) %>% 
  ungroup()

df.back2 <- df.back2 %>% group_by(subs) %>% 
  mutate(ACCz_MC = (AccZ-mean(AccZ))) %>% 
  ungroup()

df.back0 <- df.back0 %>% group_by(subs) %>% 
  mutate(CognitiveComposite_MC = (CognitiveComposite-mean(CognitiveComposite))) %>% 
  ungroup()

df.back2 <- df.back2 %>% group_by(subs) %>% 
  mutate(CognitiveComposite_MC = (CognitiveComposite-mean(CognitiveComposite))) %>% 
  ungroup()
df.back0 <- df.back0 %>% group_by(subs) %>% 
  mutate(WorkingMemory_MC = (WorkingMemory-mean(WorkingMemory))) %>% 
  ungroup()

df.back2 <- df.back2 %>% group_by(subs) %>% 
  mutate(WorkingMemory_MC = (WorkingMemory-mean(WorkingMemory))) %>% 
  ungroup()
df.back0 <- df.back0 %>% group_by(subs) %>% 
  mutate(SustainedAttention_MC = (SustainedAttention-mean(SustainedAttention))) %>% 
  ungroup()

df.back2 <- df.back2 %>% group_by(subs) %>% 
  mutate(SustainedAttention_MC = (SustainedAttention-mean(SustainedAttention))) %>% 
  ungroup()

## performance by en and stimulus type


df.back0typs1 <- df.back0 %>% group_by(subs,Type) %>% summarise(AccTp = mean(ACCz_MC))
summary(aov(AccTp~Type , data = df.back0typs1))
TukeyHSD(aov(AccTp~Type,data=df.back0typs1))

plot_type01 <- ggplot(data=df.back0typs1) + aes(x = Type, y = AccTp) + geom_violin(fill="grey70", color="grey70")+ stat_summary(fun.y=median, geom="point", shape=95, size=8, color="red") + geom_hline(aes(yintercept=0),color = "grey40") +theme_minimal(base_size = 15) + labs(y = "0-back Accuracy (z)", x = "0-back Task", color="N-back") + scale_y_continuous(breaks=c(seq(-5,4,1)), limits=c(-5,4)) + theme(panel.grid.minor=element_blank())
write_csv(df.back0typs1,"/figure files/FigS1_topleft.csv")

TukeyHSD(aov(WMTp~Type,data=df.back0typs2))

df.back0typs3 <- df.back0 %>% group_by(subs,Type) %>% summarise(SATp = mean(SustainedAttention_MC))
plot_type03 <- ggplot(data=df.back0typs3) + aes(x = Type, y = SATp) + geom_violin(fill="#E7B800", color="#E7B800")+ stat_summary(fun.y=median, geom="point", shape=95, size=8, color="red") + geom_hline(aes(yintercept=0),color = "grey40")+ theme_minimal(base_size = 15) + labs(y = "Sustained Attention Network Strength", x = "0-back Task", color="N-back") + scale_y_continuous(breaks=c(seq(-5,4,1)), limits=c(-5,4)) + theme(panel.grid.minor=element_blank())
write_csv(df.back0typs3,"/figure files/FigS1_bottomleft.csv")

summary(aov(SATp~Type,data=df.back0typs3))
TukeyHSD(aov(SATp~Type,data=df.back0typs3))

df.back2typs4 <- df.back2 %>% group_by(subs,Type) %>% summarise(AccTp = mean(ACCz_MC))
plot_type21 <- ggplot(data=df.back2typs4) + aes(x = Type, y = AccTp) + geom_violin(fill="grey70", color="grey70")+ stat_summary(fun.y=median, geom="point", shape=95, size=8, color="red") + geom_hline(aes(yintercept=0),color = "grey40")+ theme_minimal(base_size = 15) + labs(y = "2-back Accuracy (z)", x = "2-back Task", color="N-back") + scale_y_continuous(breaks=c(seq(-5,4,1)), limits=c(-5,4)) + theme(panel.grid.minor=element_blank())
write_csv(df.back2typs4,"/figure files/FigS1_topright.csv")


summary(aov(AccTp~Type,data=df.back2typs4))
TukeyHSD(aov(AccTp~Type,data=df.back2typs4))

df.back2typs5 <- df.back2 %>% group_by(subs,Type) %>% summarise(WMTp = mean(WorkingMemory_MC))
plot_type22 <- ggplot(data=df.back2typs5) + aes(x = Type, y = WMTp) + geom_violin(fill="#00AFBB", color="#00AFBB")+ stat_summary(fun.y=median, geom="point", shape=95, size=8, color="red") + geom_hline(aes(yintercept=0),color = "grey40")+ theme_minimal(base_size = 15) + labs(y = "Working Memory Network Strength", x = "2-back Task", color="N-back") + scale_y_continuous(breaks=c(seq(-5,4,1)), limits=c(-5,4)) + theme(panel.grid.minor=element_blank())
write_csv(df.back2typs5,"/figure files/FigS1_bottomright.csv")

TukeyHSD(aov(WMTp~Type,data=df.back2typs5))


ggarrange(plot_type01,plot_type21,plot_type03,plot_type22,ncol = 2,nrow=2)
#ggsave("types_4boxplots.png")

## meidiations
# WM and type in 2-back
df.back2typs5 <- df.back2 %>% group_by(subs,Type) %>% summarise(WMTp = mean(WorkingMemory_MC),ACCz_MCTp = mean(ACCz_MC),motionTp = mean(blockFD))
df.back2med1<-df.back2typs5
df.back2med1<-df.back2med1 %>% filter(Type==c("NeutFace","Place"))%>% droplevels()
med.fit <- lm(WMTp ~ Type+motionTp, data = df.back2med1)
summary(med.fit)
out.fit <- lm(ACCz_MCTp ~ WMTp+Type+motionTp,   data = df.back2med1)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "Type", mediator = "WMTp",
                     control.value = 0, treat.value = 1,sims = 500 )
summary(med23.out)

# SA and type in 0-back
df.back0med1<-df.back0
df.back0med1<-df.back0med1 %>% filter(Type==c("NeutFace","Place"))%>% droplevels()
med.fit <- lm(SustainedAttention_MC ~ Type+BlockMotion, data = df.back0med1)
summary(med.fit)
out.fit <- lm(ACCz_MC ~ SustainedAttention_MC+Type+BlockMotion,   data = df.back0med1)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "Type", mediator = "SustainedAttention_MC",
                     control.value = 0, treat.value = 1, sims = 100)
summary(med23.out)

# SA and valence in 0-back
df.back0med2<-df.back0
df.back0med2<-df.back0med2 %>% filter(Type==c("PosFace","NegFace"))%>% droplevels()
med.fit <- lm(SustainedAttention_MC ~ Type+BlockMotion, data = df.back0med2)
summary(med.fit)
out.fit <- lm(ACCz_MC ~ SustainedAttention_MC+Type+BlockMotion,   data = df.back0med2)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "Type", mediator = "SustainedAttention_MC",
                     control.value = 0, treat.value = 1, sims = 100)
summary(med23.out)

# WM and type in 2-back
df.back2med1<-df.back2typs5
df.back2med1<-df.back2med1 %>% filter(Type==c("NeutFace","Place"))%>% droplevels()
med.fit <- lm(WMTp ~ Type+BlockMotion, data = df.back2med1)
summary(med.fit)
out.fit <- lm(ACCz_MC ~ WMTp+Type+BlockMotion,   data = df.back2med1)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "Type", mediator = "WMTp",
                     control.value = 0, treat.value = 1, sims = 100)
summary(med23.out)


# WM and valence in 2-back
df.back2med2<-df.back2
df.back2med2<-df.back2med2 %>% filter(Type==c("PosFace","NegFace"))%>% droplevels()
med.fit <- lm(WorkingMemory_MC ~ Type, data = df.back2med2)
summary(med.fit)
out.fit <- lm(ACCz_MC ~ WorkingMemory_MC+Type+BlockMotion,   data = df.back2med2)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "Type", mediator = "WorkingMemory_MC",
                     control.value = 0, treat.value = 1, sims = 100)
summary(med23.out)
# supplementary: partial r^2 values
lmr2 <- lmer(AccZ ~ Age+Sex+ Motion + WorkingMemory +SustainedAttention+(1|Site), data=df.Sublevel.back2)
(partR2(lmr2,partvars = c("Age","Sex","Motion","WorkingMemory","SustainedAttention"), 
        R2_type = "marginal",max_level = 1, nboot = 100))

lmr2 <- lmer(AccZ ~ Age+Sex+ Motion + WorkingMemory +SustainedAttention+(1|Site), data=df.Sublevel.back0)
(partR2(lmr2,partvars = c("Age","Sex","Motion","WorkingMemory","SustainedAttention"), 
        R2_type = "marginal",max_level = 1, nboot = 100))

#####
########### comparison to canonical functional networks
setwd("E:/Omid/abcd/n-back_block_level_conn_analysis/manuscript/Plos_Bio/Revisions/figure files")
df.block1 <- read_csv("ABCD_nbkRun1and2_blocks_beh_brain_data.csv")
df.block2 <- read_csv("ABCD_nbkRun1and2_blocks_beh_brain_data_otherNetworks.csv")

df.block <- cbind.data.frame(df.block1,df.block2[,-1])
df.block_FDfilter <- filter(df.block, FDav < .2, FDmax < 2, nbkqc_eprm ==0 ,laptop_nback ==1, boxswitch ==0) #exclude runs with high motion etc.

df.block_FDfilter$Sex <- as.factor(df.block_FDfilter$Sex)
df.block_FDfilter$N <- as.factor(df.block_FDfilter$N)
df.block_FDfilter$Type <- as.factor(df.block_FDfilter$Type)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)
df.block_FDfilter$Site <- as.factor(df.block_FDfilter$Site)

df.block_FDfilter <- filter(df.block_FDfilter, N != "NaN", blockFD != "NaN") %>% droplevels() # 
(SitesN <- df.block_FDfilter %>% distinct(subs, .keep_all = TRUE) %>% group_by(Site) %>% count() )
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
length(unique(df.block_FDfilter$subs)) #1545 subs and 19216 blokcs
filter(df.block_FDfilter, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores
df.block_FDfilter$Type <- factor(df.block_FDfilter$Type, levels = c("NeutFace", "NegFace", "PosFace", "Place"))
df.block_FDfilter$Sex <- factor(df.block_FDfilter$Sex,  levels = c("M", "F"))
df.block_FDfilter$Sex<-recode(df.block_FDfilter$Sex, `F` = "Female", M = "Male")
df.block_FDfilter$blockOrdered <- factor(df.block_FDfilter$block, ordered = T)

df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

colnames <- colnames(df.back0s)
colnames_canons <- colnames[34:41]

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), brainscNIHcpm = mean(brainscNIHcpm), across(brainsc_medialprefrontal:brainsc_visualassociation, mean))

# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(across(brainsc_medialprefrontal:brainsc_visualassociation, ~ (.x - mean(.x))/sd(.x)))


df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN")

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm), across(brainsc_medialprefrontal:brainsc_visualassociation, mean))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(across(brainsc_medialprefrontal:brainsc_visualassociation, ~ (.x - mean(.x))/sd(.x)))


betas_0back <- data.frame()
for (i in 1:8){
  model_formula <- as.formula(paste0("AccZ ~ Age+Sex+ Motion + (1|Site) + ",colnames_canons[i]))
  
  model_fit <- lmer(data = df.Sublevel.back0,
                    formula = model_formula)
  betas_0back <- rbind(betas_0back, c(model_fit@beta[5], summary(model_fit)$coef[5, 2, drop = TRUE]))
}

betas_2back <- data.frame()
for (i in 1:8){
  model_formula <- as.formula(paste0("AccZ ~ Age+Sex+ Motion + (1|Site) + ",colnames_canons[i]))
  
  model_fit <- lmer(data = df.Sublevel.back2,
                    formula = model_formula)
  betas_2back <- rbind(betas_2back, c(model_fit@beta[5], summary(model_fit)$coef[5, 2, drop = TRUE]))
}


lmsaa <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory + (1|Site), data=df.Sublevel.back2)
lmsbb <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory + (1|Site), data=df.Sublevel.back0)


lmsa <- lmer(AccZ ~ Age+Sex+ Motion +SustainedAttention + (1|Site), data=df.Sublevel.back2)
lmsb <- lmer(AccZ ~ Age+Sex+ Motion +SustainedAttention + (1|Site), data=df.Sublevel.back0)

# z-test of beta coefficents based on Clogg, C. C., Petkova, E., & Haritou, A. (1995). Statistical methods for comparing regression coefficients between models. American Journal of Sociology, 100(5), 1261-1293.
k = 8 # set network 1 through 8
se1 = summary(lmsbb)$coef[5, 2, drop = TRUE]   
se2 = betas_0back[k,2]
beta1 = lmsbb@beta[5]
beta2 = betas_0back[k,1]
Z0wm = (beta1 - beta2)/(sqrt(se1^2 + se2^2))

se1 = summary(lmsaa)$coef[5, 2, drop = TRUE]  # 2-back and WM
se2 = betas_2back[k,2]
beta1 = lmsaa@beta[5]
beta2 = betas_2back[k,1]
(Z2wm = (beta1 - beta2)/(sqrt(se1^2 + se2^2)))

se1 = summary(lmsb)$coef[5, 2, drop = TRUE] # 0-back and SA
se2 = betas_0back[k,2]
beta1 = lmsb@beta[5]
beta2 = betas_2back[1,1]
(Z0sa = (beta1 - beta2)/(sqrt(se1^2 + se2^2)))

se1 = summary(lmsa)$coef[5, 2, drop = TRUE]   
se2 = betas_2back[k,2]
beta1 = lmsa@beta[5]
beta2 = betas_2back[k,1]
Z2sa = (beta1 - beta2)/(sqrt(se1^2 + se2^2))

2*pnorm(q=2.389, lower.tail=FALSE)  # smallest wm-2 difference for network 7 still significant p =.017
2*pnorm(q=4.738, lower.tail=FALSE)  # smallest Sa-0 difference for network 7 still significant p =.017


########### comparison to 0-back cpm trained in HCP dataset
df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA),
            brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav),
            FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp),
            Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm),
            brainscadul0bk = mean(brainscHCP0BKcpm))

# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(adul0bkCPM = (brainscadul0bk-mean(brainscadul0bk))/sd(brainscadul0bk))


df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN")

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA),
            brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav),
            FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp),
            Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm),
            brainscadul0bk = mean(brainscHCP0BKcpm))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(adul0bkCPM = (brainscadul0bk-mean(brainscadul0bk))/sd(brainscadul0bk))

lmsaa <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory + (1|Site), data=df.Sublevel.back2)
lmsbb <- lmer(AccZ ~ Age+Sex+ Motion +adul0bkCPM + (1|Site), data=df.Sublevel.back2)

tab_model(lmsaa,lmsbb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","2-back Accuracy"))


lmsa <- lmer(AccZ ~ Age+Sex+ Motion +SustainedAttention + (1|Site), data=df.Sublevel.back0)
lmsb <- lmer(AccZ ~ Age+Sex+ Motion +adul0bkCPM + (1|Site), data=df.Sublevel.back0)

tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("0-back Accuracy","0-back Accuracy"))

gg<-pcor(df.Sublevel.back2[,c("AccZ","WorkingMemory","Age","Motion")])
gg$estimate
gg$p.value

gg<-pcor(df.Sublevel.back0[,c("AccZ","adul0bkCPM","Age","Motion")])
gg$estimate
gg$p.value

##### test difference between betas
#Function to get difference between bootstrapped coefficients 0bkcpm and 2bkcpm
boot_func <- function(data0,data2){
  resampled_data0 = data0 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA <- lmer(AccZ ~ Age + Sex + Motion + adul0bkCPM +(1|Site), data=resampled_data0)
  lmsa_resampledWM <- lmer(AccZ ~ Age + Sex + Motion + WorkingMemory + (1|Site), data=resampled_data2)
  
  unname(lmsa_resampledSA@beta[5]) - unname(lmsa_resampledWM@beta[5])
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func(df.Sublevel.back0,df.Sublevel.back2))
#histogram of these permuted correlations
hist(boot_coefs)
(pval = 1- sum(boot_coefs>0)/length(boot_coefs))



########### Random nets analysis

df.block <- read_csv("ABCD_nbkRun1and2_blocks_march21b_famreclapbox_RandNets.csv")
df.block_FDfilter <- filter(df.block, FDav < .2, FDmax < 2, nbkqc_eprm ==0 ,laptop_nback ==1, boxswitch ==0) #exclude runs with high motion etc.
df.block_FDfilter$Sex <- as.factor(df.block_FDfilter$Sex)
df.block_FDfilter$N <- as.factor(df.block_FDfilter$N)
df.block_FDfilter$Type <- as.factor(df.block_FDfilter$Type)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)
df.block_FDfilter$Site <- as.factor(df.block_FDfilter$Site)

summary(df.block_FDfilter$N) # NaN
df.block_FDfilter <- filter(df.block_FDfilter, N != "NaN", blockFD != "NaN") %>% droplevels() # 
(SitesN <- df.block_FDfilter %>% distinct(subs, .keep_all = TRUE) %>% group_by(Site) %>% count() )
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
length(unique(df.block_FDfilter$subs)) #1545 subs and 19216 blokcs


filter(df.block_FDfilter, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores

df.block_FDfilter$Type <- factor(df.block_FDfilter$Type, levels = c("NeutFace", "NegFace", "PosFace", "Place"))
df.block_FDfilter$Sex <- factor(df.block_FDfilter$Sex,  levels = c("M", "F"))
df.block_FDfilter$Sex<-recode(df.block_FDfilter$Sex, `F` = "Female", M = "Male")
df.block_FDfilter$blockOrdered <- factor(df.block_FDfilter$block, ordered = T)


df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

colnames <- colnames(df.back0s)
colnames_rand <- colnames[33:232]

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp) , brainscNIHcpm = mean(brainscNIHcpm), across(brainsc_Rand_1:brainsc_Rand_200, mean))

# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(across(brainsc_Rand_1:brainsc_Rand_200, ~ (.x - mean(.x))/sd(.x)))

df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN")

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm), across(brainsc_Rand_1:brainsc_Rand_200, mean))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(across(brainsc_Rand_1:brainsc_Rand_200, ~ (.x - mean(.x))/sd(.x)))

betas_0back <- data.frame()
for (i in 1:200){
  model_formula <- as.formula(paste0("AccZ ~ Age+Sex+ Motion + (1|Site) + ",colnames_rand[i]))
  
  model_fit <- lmer(data = df.Sublevel.back0,
                    formula = model_formula)
  betas_0back <- rbind(betas_0back, model_fit@beta[5])
}

betas_2back <- data.frame()
for (i in 1:200){
  model_formula <- as.formula(paste0("AccZ ~ Age+Sex+ Motion + (1|Site) + ",colnames_rand[i]))
  
  model_fit <- lmer(data = df.Sublevel.back2,
                    formula = model_formula)
  betas_2back <- rbind(betas_2back, model_fit@beta[5])
}

############# removing siblings analysis
# making sub-level dataframe and pick one sub from sharing famid, then expanding back to block-level

df.Sublevel <- df.block_FDfilter %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIHcpm = mean(brainscNIHcpm) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), recmem = first(allrecmem_dp), Site = first(Site), famid = first(famid) ) %>% ungroup()
df.Sublevel$famid <- as.factor(df.Sublevel$famid)

#no sibs in study (664)
df.Sublevel.noSibs <- df.Sublevel %>% filter(famid == "NaN")

#famid with more than 1 sib left
multSibleft <- df.Sublevel %>% group_by(famid) %>% count() %>% filter(n != 1) %>% droplevels()
df.Sublevel.multSibleft <- df.Sublevel %>% filter(famid %in% multSibleft$famid, famid != "NaN")

#only one sib 
df.Sublevel.onlySibleft <- df.Sublevel %>% filter(!(famid %in% multSibleft$famid))

#pulling randomly from the families with 2 siblings
df.Sublevel.multSibleft_choose1 <- df.Sublevel.multSibleft %>% distinct(famid, .keep_all = T)

#Binding back together
df.Sublevel.OneSubPerFam <- rbind(df.Sublevel.noSibs, df.Sublevel.onlySibleft, df.Sublevel.multSibleft_choose1)

df.Sublevel <- df.Sublevel %>% mutate(IncludedOnePerFamID = if_else(subs %in% df.Sublevel.OneSubPerFam$subs, 1, 0))

df.block_FDfilter <- df.block_FDfilter %>% mutate(IncludedOnePerFamID = if_else(subs %in% df.Sublevel.OneSubPerFam$subs, 1, 0))

#################BETWEEN-SUBJECT DATA AND PLOTS#######################
df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s, avnih5 != "NaN",Acc !="NaN", IncludedOnePerFamID == 1) 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIHcpm = mean(brainscNIHcpm) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site) )
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(avnih5 = (avnih5-mean(avnih5))/sd(avnih5))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s, avnih5 != "NaN",Acc !="NaN", IncludedOnePerFamID == 1) 

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIHcpm = mean(brainscNIHcpm) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site) )
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(avnih5 = (avnih5-mean(avnih5))/sd(avnih5))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

# mixed effects regressions across subs
lms1 <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory+SustainedAttention +(1|Site), data=df.Sublevel.back2)
summary(lms1)

lms2 <- lmer(AccZ ~ Age+Sex+ Motion +WorkingMemory+SustainedAttention + (1|Site), data=df.Sublevel.back0)
summary(lms2)

tab_model(lms2,lms1, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("0-back Accuracy","2-back Accuracy"))
cor.test(df.Sublevel.back2$Acc, df.Sublevel.back2$WorkingMemory)

lms1 <- lmer(AccZ ~ Age+Sex+ Motion +CognitiveCompositeCPM+SustainedAttention +(1|Site), data=df.Sublevel.back2)
summary(lms1)

lms2 <- lmer(AccZ ~ Age+Sex+ Motion +CognitiveCompositeCPM+SustainedAttention + (1|Site), data=df.Sublevel.back0)
summary(lms2)
tab_model(lms2,lms1, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("0-back Accuracy","2-back Accuracy"))


##########RecMem
# collapse across both 0-back and 2-back blocks to test rec mem across the entire n-back
df.back02s <- df.block_FDfilter
df.back02s <- filter(df.back02s, Acc !="NaN", IncludedOnePerFamID == 1) 

df.Sublevel.back02 <- df.back02s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIHcpm = mean(brainscNIHcpm) ,brainscAveryTask = mean(brainscAveryTask),
            FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back02 <- df.Sublevel.back02 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))
lms1 <- lmer(recmem ~ Age+Sex+ Motion +AccZ+WorkingMemory+SustainedAttention +(1|Site), data=df.Sublevel.back02)
summary(lms1)
tab_model(lms1,show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("Recognition Memory d'"))

#############
##################### HCP dataset  ####################

df.block <- read_csv("/HCP_nbkRunRLLR_blocks_beh_brain_data.csv")
demog <- read_csv("/HCP1200/unrestricted_okardan_5_12_2021_14_45_23.csv")
df.bl.withdemo <- merge(df.block,demog[,c('Subject','Gender','Age','QC_Issue')],by.x = "subs"
                        ,by.y = "Subject",all.x = TRUE)

df.Sublevel <- df.bl.withdemo %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), 
             brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs), motion = mean(blockFD), FDmax = mean(FDmaxs)
             )

df.block_FDfilter <- filter(df.bl.withdemo, FDavs < .2, FDmaxs < 2, is.na(QC_Issue)) #exclude runs with high motion etc.

df.block_FDfilter$Type <- as.factor(df.block_FDfilter$blockType)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)

summary(df.block_FDfilter$Ns) # NaN

df.block_FDfilter <- filter(df.block_FDfilter, Ns != "NaN") %>% droplevels() # 

length(unique(df.block_FDfilter$subs)) #881 subs and 12648 blokcs or 754 subs 10824 blocks

filter(df.block_FDfilter, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores
df.block_FDfilter$Gender <- factor(df.block_FDfilter$Gender,  levels = c("M", "F"))

#################BETWEEN-SUBJECT DATA
df.back0s <- filter(df.block_FDfilter, Ns == '0')
df.back0s <- filter(df.back0s,Accs !="NaN") 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs), motion = mean(blockFD), FDmax = mean(FDmaxs)
             )
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


df.back2s <- filter(df.block_FDfilter, Ns == '2')
#df.back2s <- filter(df.back2s, avnih5 != "NaN",Acc !="NaN") 
df.back2s <- filter(df.back2s,Accs !="NaN") 

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs),motion = mean(blockFD), FDmax = mean(FDmaxs))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

select(df.Sublevel.back2, Age ) %>% group_by(Age) %>% count()
select(df.Sublevel.back0, Age ) %>% group_by(Age) %>% count()
select(df.Sublevel.back0, Gender ) %>% group_by(Gender) %>% count()

# regressions across subs
lmsa <- lm(AccZ ~  Age + Gender+ Motion +cSA, data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lm(AccZ ~ Motion +cSA , data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","0-back Accuracy"))

lmsa <- lm(AccZ ~  Age + Gender+Motion +SustainedAttention, data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lm(AccZ ~ Age + Gender+Motion +SustainedAttention , data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsb,lmsa, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("0-back Accuracy","2-back Accuracy"))

##### test difference between betas
#Adding 2 back and 0 back together
df.Sublevel.back0$N <- 0
df.Sublevel.back2$N <- 2


colnames(df.Sublevel.back0) <-paste(colnames(df.Sublevel.back0),"0",sep="_")
colnames(df.Sublevel.back2) <-paste(colnames(df.Sublevel.back2),"2",sep="_")
df.Sublevel.all <- merge(df.Sublevel.back0, df.Sublevel.back2, by.x = "subs_0", by.y = "subs_2")

#Function to get difference between bootstrapped coefficients
boot_func <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampled2 <- lm(AccZ_2 ~ Age_2 + Gender_2 + Motion_2 + SustainedAttention_2, data=resampled_data)
  lmsa_resampled0 <- lm(AccZ_0 ~ Age_0 + Gender_0 + Motion_0 + SustainedAttention_0, data=resampled_data)
  
  unname(lmsa_resampled0$coefficients["SustainedAttention_0"]) - unname(lmsa_resampled2$coefficients["SustainedAttention_2"])
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func(df.Sublevel.all))
#histogram of these permuted correlations
hist(boot_coefs)
(pval = 1- sum(boot_coefs>0)/length(boot_coefs))

######
lmsa <- lm(AccZ ~  Motion +WorkingMemory, data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lm(AccZ ~ Motion +WorkingMemory , data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","0-back Accuracy"))

lmsa <- lm(AccZ ~  Motion +SustainedAttention, data=df.Sublevel.back2)
summary(lmsa)
lmsb <- lm(AccZ ~ Motion +SustainedAttention, data=df.Sublevel.back0)
summary(lmsb)
tab_model(lmsa,lmsb, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("2-back Accuracy","0-back Accuracy"))

cor.test(df.Sublevel.back2$SustainedAttention, df.Sublevel.back0$SustainedAttention)
cor.test(df.Sublevel.back2$AccZ, df.Sublevel.back0$AccZ)
cor.test(df.Sublevel.back2$SustainedAttention, df.Sublevel.back0$AccZ)
cor.test(df.Sublevel.back2$AccZ, df.Sublevel.back0$SustainedAttention)
r.test(754, r12 = .066, r34 = .168, r23 = .030, r13 = .588, r14 = .131, r24 = .407, 
       n2 = NULL,pooled=TRUE, twotailed = TRUE)

#plots

cor.test(df.Sublevel.back2$Acc, df.Sublevel.back2$SustainedAttention, method = "spearman")
p2 <- ggplot(data=df.Sublevel.back2) + 
  aes(x = SustainedAttention, y = Acc*100) + 
  geom_point(alpha=.15, size = 2, color = "#E7B800") +
  geom_smooth(method = "lm", aes(x = cSA, y = Acc*100), size=1, color = "#E7B800") +
  annotate("text", label = "r = .066, p = .070", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = "2-back Accuracy (%) ", x = "Sustained Attention Network Strength", color="N-back", title = " ") 
p2
FigS6_tab2 <- data.frame(df.Sublevel.back2$SustainedAttention,df.Sublevel.back2$Acc*100)
write_csv(FigS6_tab2,"/figure files/FigS6_right.csv")

cor.test(df.Sublevel.back0$Acc, df.Sublevel.back0$SustainedAttention, method = "pearson")
p5 <- ggplot(data=df.Sublevel.back0) + 
  aes(x = SustainedAttention, y = Acc*100) + 
  geom_point(alpha=.25, size = 2, color = "#E7B800") +
  geom_smooth(method = "lm", aes(x = cSA, y = Acc*100), size=1, color = "#E7B800") +
  annotate("text", label = "r = .168, p < .001", x = 0, y = 45, size=5) +
  theme_minimal(base_size = 15) + coord_cartesian(xlim = c(-5,5), ylim = c(0,100))+
  labs(y = "0-back Accuracy (%) ", x = "Sustained Attention Network Strength", color="N-back", title = " ") 
p5
FigS6_tab1 <- data.frame(df.Sublevel.back0$SustainedAttention,df.Sublevel.back0$Acc*100)
write_csv(FigS6_tab1,"E:/figure files/FigS6_left.csv")


ggarrange(p5,p2,ncol = 2,nrow=1)


###################WITHIN-SUBJECT

# z-scoring variables so that beta coefficients become comparable
df.back0 <- filter(df.block_FDfilter, Ns == '0')
df.back0 <- filter(df.back0,Accs !="NaN") 
df.back0 <- df.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.back0 <- df.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.back0 <- df.back0 %>% mutate(BlockMotion = (blockFD-mean(blockFD))/sd(blockFD))
df.back0 <- df.back0 %>% mutate(AccZ = (Accs-mean(Accs))/sd(Accs))
length(unique(df.back0$subs))


df.back2 <- filter(df.block_FDfilter, Ns == '2')
df.back2 <- filter(df.back2,Accs !="NaN") 
df.back2 <- df.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.back2 <- df.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.back2 <- df.back2 %>% mutate(BlockMotion = (blockFD-mean(blockFD))/sd(blockFD))
df.back2 <- df.back2 %>% mutate(AccZ = (Accs-mean(Accs))/sd(Accs))


###### regressions within-sub predicting acc

lm11 <- lmer(AccZ ~  Run+blockFD+Type+SustainedAttention+(1|subs), data=df.back2)
lm12 <- lmer(AccZ ~  Run+blockFD+Type+SustainedAttention+(1|subs), data=df.back0)
tab_model(lm12,lm11,  show.stat = TRUE, string.stat = "t-Statistic",
          show.se = FALSE, string.se = "SE",dv.labels =c( "0-back Accuracy","2-back Accuracy"),digits = 2, digits.p = 3)

####################### Fig 6 mixed HCP and BACD
#ABCD
df.block <- read_csv("ABCD_nbkRun1and2_blocks_march21b_famreclapbox_totfc.csv")
df.block$brainscSA <- df.block$SA_z  # within-individual z-scored network strength values
df.block_FDfilter <- filter(df.block, FDav < .2, FDmax < 2, nbkqc_eprm ==0 ,laptop_nback ==1, boxswitch ==0) #exclude runs with high motion etc.

df.block_FDfilter$Sex <- as.factor(df.block_FDfilter$Sex)
df.block_FDfilter$N <- as.factor(df.block_FDfilter$N)
df.block_FDfilter$Type <- as.factor(df.block_FDfilter$Type)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)
df.block_FDfilter$Site <- as.factor(df.block_FDfilter$Site)

summary(df.block_FDfilter$N) # NaN

df.block_FDfilter <- filter(df.block_FDfilter, N != "NaN", blockFD != "NaN") %>% droplevels() # 
(SitesN <- df.block_FDfilter %>% distinct(subs, .keep_all = TRUE) %>% group_by(Site) %>% count() )
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
length(unique(df.block_FDfilter$subs)) #1545 subs and 19216 blokcs

brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores

df.block_FDfilter$Type <- factor(df.block_FDfilter$Type, levels = c("NeutFace", "NegFace", "PosFace", "Place"))
df.block_FDfilter$Sex <- factor(df.block_FDfilter$Sex,  levels = c("M", "F"))
df.block_FDfilter$Sex<-recode(df.block_FDfilter$Sex, `F` = "Female", M = "Male")
df.block_FDfilter$blockOrdered <- factor(df.block_FDfilter$block, ordered = T)

df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp) , brainscNIHcpm = mean(brainscNIHcpm))

# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN") 

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

# HCP
df.blockhcp <- read_csv("HCP_nbkRunRLLR_blocks_jun21_totfc.csv")
df.blockhcp$brainscSA <- df.blockhcp$SA_z  # within-individual z-scored network strength values

demog <- read_csv("/HCP1200/unrestricted_okardan_5_12_2021_14_45_23.csv")
df.bl.withdemo <- merge(df.blockhcp,demog[,c('Subject','Gender','Age','QC_Issue')],by.x = "subs"
                        ,by.y = "Subject",all.x = TRUE)
df.block_FDfilterhcp <- filter(df.bl.withdemo, FDavs < .2, FDmaxs < 2, is.na(QC_Issue)) #exclude runs with high motion etc.
df.block_FDfilterhcp$Type <- as.factor(df.block_FDfilterhcp$blockType)
df.block_FDfilterhcp$Run <- as.factor(df.block_FDfilterhcp$Run)

summary(df.block_FDfilterhcp$Ns) # NaN

df.block_FDfilterhcp <- filter(df.block_FDfilterhcp, Ns != "NaN") %>% droplevels() # 
length(unique(df.block_FDfilterhcp$subs)) #881 subs and 12648 blokcs or 754 subs 10824 blocks

filter(df.block_FDfilterhcp, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilterhcp, brainscSA == "NaN") # check for NaN SA network scores
df.block_FDfilterhcp$Gender <- factor(df.block_FDfilterhcp$Gender,  levels = c("M", "F"))

df.back0shcp <- filter(df.block_FDfilterhcp, Ns == '0')
df.back0shcp <- filter(df.back0shcp,Accs !="NaN") 

df.Sublevel.back0hcp <- df.back0shcp %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs), motion = mean(blockFD), FDmax = mean(FDmaxs))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
f.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

df.back2shcp <- filter(df.block_FDfilterhcp, Ns == '2')
df.back2shcp <- filter(df.back2shcp,Accs !="NaN") 

df.Sublevel.back2hcp <- df.back2shcp %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs),motion = mean(blockFD), FDmax = mean(FDmaxs))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

# mixed graph

totaldf <- rbind(as.data.frame(cbind(df.Sublevel.back0$Acc,df.Sublevel.back0$brainscSA))
                 , as.data.frame(cbind(df.Sublevel.back0hcp$Acc,df.Sublevel.back0hcp$brainscSA)))
M<- mean(totaldf$V2)
SD<- sd(totaldf$V2)
MM<- mean(totaldf$V1)
SDSD<- sd(totaldf$V1)
Myouth<- mean(df.Sublevel.back0$brainscSA)
Myouth2<- mean(df.Sublevel.back0$Acc)
Madult<- mean(df.Sublevel.back0hcp$brainscSA)
Madult2<- mean(df.Sublevel.back0hcp$Acc)

t.test(df.Sublevel.back0hcp$brainscSA, y = df.Sublevel.back0$brainscSA)
grad <- as.data.frame(cbind(c((Myouth-M)/SD,(Madult-M)/SD),c((Myouth2-0)/1,(Madult2-0)/1)))

#install.packages("ggExtra")
library(ggExtra)
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(dataset="1")
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(dataset="2")
p12b <- ggplot(NULL,aes(x = (brainscSA-M)/SD, y = (Acc)*100/1, color= dataset)) + 
  geom_point(alpha=.20, size = 4, data=df.Sublevel.back0) +
  geom_smooth(method = "lm", aes(x = (brainscSA-M)/SD, y = (Acc)*100/1), size=1, data=df.Sublevel.back0)+
  scale_color_manual(values = c("#E7B800"), aesthetics = "color") +
  theme_minimal(base_size = 21) + coord_cartesian(xlim = c(-4,3), ylim = c(0,103))+
  theme(legend.position = "none") +
  labs(y = "0-back Accuracy (%)", x = 'Sustained attention network strength (z)', color="", title = "") 

p12b

df.Sublevel.back0_abcd_temp <- df.Sublevel.back0 %>% dplyr::select(brainscSA, Acc, dataset)
df.Sublevel.back0hcp_temp <- df.Sublevel.back0hcp %>% dplyr::select(brainscSA, Acc, dataset)
df.Sublevel.back0_abcd_hcp <- rbind(df.Sublevel.back0_abcd_temp, df.Sublevel.back0hcp_temp)


p12 <- ggplot(df.Sublevel.back0_abcd_hcp,aes(x = (brainscSA-M)/SD, y = (Acc)*100/1, color= dataset, group=dataset)) + 
  geom_point(alpha=.20, size = 4) +
  geom_smooth(method = "lm")+
 geom_point(data = grad,aes(x = V1,
                             y = V2*100, size=1),shape = 3, size=3, color = "red", inherit.aes = FALSE)+
  scale_color_manual(values = c("#E7B800","#c6a2c6"), aesthetics = "color") +
  theme_minimal(base_size = 16) + coord_cartesian(xlim = c(-4,3), ylim = c(0,103))+
  theme(legend.position = "none") +
  labs(y = "0-back Accuracy (%)", x = 'Sustained attention network strength (z)', color="", title = "") 

p12
Fig6_dat <- data.frame((df.Sublevel.back0_abcd_hcp$brainscSA - M)/SD,df.Sublevel.back0_abcd_hcp$Acc*100)
write_csv(Fig6_dat,"/figure files/Fig6_dat.csv")
ggsave(file = "Fig6.png", ggMarginal(p12, groupColour = TRUE, groupFill = TRUE), width = 7, height = 7)
t.test(df.Sublevel.back0$brainscSA, y = df.Sublevel.back0hcp$brainscSA,pair=FALSE)
t.test(df.Sublevel.back0$Acc, y = df.Sublevel.back0hcp$Acc,pair=FALSE)

############################ lesioning SA networks analysis
setwd("/SA_lesions_analysis")

library(tidyverse)
library(RColorBrewer)
library(stargazer)
library(lme4)
library(sjPlot)
library(ggplot2)
library(ggpubr)
library(corrplot)

df.block <- read_csv("ABCD_nbkRun1and2_blocks_march21b_famreclapbox_Lesions.csv")
df.block_FDfilter <- filter(df.block, FDav < .2, FDmax < 2, nbkqc_eprm ==0 ,laptop_nback ==1, boxswitch ==0) #exclude runs with high motion etc.
df.block_FDfilter$Sex <- as.factor(df.block_FDfilter$Sex)
df.block_FDfilter$N <- as.factor(df.block_FDfilter$N)
df.block_FDfilter$Type <- as.factor(df.block_FDfilter$Type)
df.block_FDfilter$Run <- as.factor(df.block_FDfilter$Run)
df.block_FDfilter$Site <- as.factor(df.block_FDfilter$Site)

summary(df.block_FDfilter$N) # NaN
df.block_FDfilter <- filter(df.block_FDfilter, N != "NaN", blockFD != "NaN") %>% droplevels() # 
(SitesN <- df.block_FDfilter %>% distinct(subs, .keep_all = TRUE) %>% group_by(Site) %>% count() )
df.block_FDfilter <- filter(df.block_FDfilter, Site != "22") %>% droplevels() # a site with only 3 subs should be removed
ength(unique(df.block_FDfilter$subs)) #1545 subs and 19216 blokcs
filter(df.block_FDfilter, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilter, brainscSA == "NaN") # check for NaN SA network scores
df.block_FDfilter$Type <- factor(df.block_FDfilter$Type, levels = c("NeutFace", "NegFace", "PosFace", "Place"))
df.block_FDfilter$Sex <- factor(df.block_FDfilter$Sex,  levels = c("M", "F"))
df.block_FDfilter$Sex<-recode(df.block_FDfilter$Sex, `F` = "Female", M = "Male")
df.block_FDfilter$blockOrdered <- factor(df.block_FDfilter$block, ordered = T)
df.back0s <- filter(df.block_FDfilter, N == '0')
df.back0s <- filter(df.back0s,Acc !="NaN") 

df.Sublevel.back0 <- df.back0s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp) , brainscNIHcpm = mean(brainscNIHcpm), brainscSA_Prefrntal = mean(brainscSA_Prefrntal),brainscSA_Motor = mean(brainscSA_Motor),brainscSA_Insula = mean(brainscSA_Insula),brainscSA_Parietal = mean(brainscSA_Parietal),brainscSA_Temporal = mean(brainscSA_Temporal),brainscSA_Occipital = mean(brainscSA_Occipital),brainscSA_Limbic = mean(brainscSA_Limbic), brainscSA_Cerebellum = mean(brainscSA_Cerebellum), brainscSA_Subcortex = mean(brainscSA_Subcortex), brainscSA_Brainstem = mean(brainscSA_Brainstem))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_prefrontal = (brainscSA_Prefrntal-mean(brainscSA_Prefrntal))/sd(brainscSA_Prefrntal))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_motor = (brainscSA_Motor-mean(brainscSA_Motor))/sd(brainscSA_Motor))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_insula = (brainscSA_Insula-mean(brainscSA_Insula))/sd(brainscSA_Insula))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_parietal = (brainscSA_Parietal-mean(brainscSA_Parietal))/sd(brainscSA_Parietal))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_temporal = (brainscSA_Temporal-mean(brainscSA_Temporal))/sd(brainscSA_Temporal))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_occipital = (brainscSA_Occipital-mean(brainscSA_Occipital))/sd(brainscSA_Occipital))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_limbic = (brainscSA_Limbic-mean(brainscSA_Limbic))/sd(brainscSA_Limbic))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_cerebellum = (brainscSA_Cerebellum-mean(brainscSA_Cerebellum))/sd(brainscSA_Cerebellum))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_subcortex = (brainscSA_Subcortex-mean(brainscSA_Subcortex))/sd(brainscSA_Subcortex))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(SA_brainstem = (brainscSA_Brainstem-mean(brainscSA_Brainstem))/sd(brainscSA_Brainstem))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


df.back2s <- filter(df.block_FDfilter, N == '2')
df.back2s <- filter(df.back2s,Acc !="NaN") 

df.Sublevel.back2 <- df.back2s %>% 
  group_by(subs) %>% 
  summarise(Age = mean(Age), Sex = first(Sex), Acc = mean(Acc), RT = mean(RT), brainscSA = mean(brainscSA), brainscNIH = mean(brainscNIH) ,brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDav), FDmax = mean(FDmax), avnih5 = mean(avnih5), Motion = mean(blockFD),numResp = mean(numResp), Site = first(Site), recmem = mean(allrecmem_dp), brainscNIHcpm = mean(brainscNIHcpm), brainscSA_Prefrntal = mean(brainscSA_Prefrntal),brainscSA_Motor = mean(brainscSA_Motor),brainscSA_Insula = mean(brainscSA_Insula),brainscSA_Parietal = mean(brainscSA_Parietal),brainscSA_Temporal = mean(brainscSA_Temporal),brainscSA_Occipital = mean(brainscSA_Occipital),brainscSA_Limbic = mean(brainscSA_Limbic), brainscSA_Cerebellum = mean(brainscSA_Cerebellum), brainscSA_Subcortex = mean(brainscSA_Subcortex), brainscSA_Brainstem = mean(brainscSA_Brainstem))

# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_prefrontal = (brainscSA_Prefrntal-mean(brainscSA_Prefrntal))/sd(brainscSA_Prefrntal))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_motor = (brainscSA_Motor-mean(brainscSA_Motor))/sd(brainscSA_Motor))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_insula = (brainscSA_Insula-mean(brainscSA_Insula))/sd(brainscSA_Insula))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_parietal = (brainscSA_Parietal-mean(brainscSA_Parietal))/sd(brainscSA_Parietal))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_temporal = (brainscSA_Temporal-mean(brainscSA_Temporal))/sd(brainscSA_Temporal))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_occipital = (brainscSA_Occipital-mean(brainscSA_Occipital))/sd(brainscSA_Occipital))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_limbic = (brainscSA_Limbic-mean(brainscSA_Limbic))/sd(brainscSA_Limbic))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_cerebellum = (brainscSA_Cerebellum-mean(brainscSA_Cerebellum))/sd(brainscSA_Cerebellum))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_subcortex = (brainscSA_Subcortex-mean(brainscSA_Subcortex))/sd(brainscSA_Subcortex))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(SA_brainstem = (brainscSA_Brainstem-mean(brainscSA_Brainstem))/sd(brainscSA_Brainstem))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Motion = (Motion-mean(Motion))/sd(Motion))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(Age = (Age-mean(Age))/sd(Age))
df.Sublevel.back2 <- df.Sublevel.back2 %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

#HCP
df.blockhcp <- read_csv("HCP1200/HCP_nbkRunRLLR_blocks_jun21_Lesions.csv")
demog <- read_csv("HCP1200/unrestricted_okardan_5_12_2021_14_45_23.csv")
df.bl.withdemo <- merge(df.blockhcp,demog[,c('Subject','Gender','Age','QC_Issue')],by.x = "subs"
                        ,by.y = "Subject",all.x = TRUE)
df.block_FDfilterhcp <- filter(df.bl.withdemo, FDavs < .2, FDmaxs < 2, is.na(QC_Issue)) #exclude runs with high motion etc.
df.block_FDfilterhcp$Type <- as.factor(df.block_FDfilterhcp$blockType)
df.block_FDfilterhcp$Run <- as.factor(df.block_FDfilterhcp$Run)
summary(df.block_FDfilterhcp$Ns) # NaN
df.block_FDfilterhcp <- filter(df.block_FDfilterhcp, Ns != "NaN") %>% droplevels() # 
length(unique(df.block_FDfilterhcp$subs)) # 754 subs 10824 blocks
filter(df.block_FDfilterhcp, brainscSA == "NaN") %>% count()
brainNaNs <- filter(df.block_FDfilterhcp, brainscSA == "NaN") # check for NaN SA network scores
df.block_FDfilterhcp$Gender <- factor(df.block_FDfilterhcp$Gender,  levels = c("M", "F"))
df.back0shcp <- filter(df.block_FDfilterhcp, Ns == '0')
df.back0shcp <- filter(df.back0shcp,Accs !="NaN") 

df.Sublevel.back0hcp <- df.back0shcp %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs), motion = mean(blockFD), FDmax = mean(FDmaxs), 
             brainscSA_Prefrntal = mean(brainscSA_Prefrntal),brainscSA_Motor = mean(brainscSA_Motor),brainscSA_Insula = mean(brainscSA_Insula),brainscSA_Parietal = mean(brainscSA_Parietal),brainscSA_Temporal = mean(brainscSA_Temporal),brainscSA_Occipital = mean(brainscSA_Occipital),brainscSA_Limbic = mean(brainscSA_Limbic), brainscSA_Cerebellum = mean(brainscSA_Cerebellum), brainscSA_Subcortex = mean(brainscSA_Subcortex), brainscSA_Brainstem = mean(brainscSA_Brainstem))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_prefrontal = (brainscSA_Prefrntal-mean(brainscSA_Prefrntal))/sd(brainscSA_Prefrntal))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_motor = (brainscSA_Motor-mean(brainscSA_Motor))/sd(brainscSA_Motor))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_insula = (brainscSA_Insula-mean(brainscSA_Insula))/sd(brainscSA_Insula))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_parietal = (brainscSA_Parietal-mean(brainscSA_Parietal))/sd(brainscSA_Parietal))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_temporal = (brainscSA_Temporal-mean(brainscSA_Temporal))/sd(brainscSA_Temporal))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_occipital = (brainscSA_Occipital-mean(brainscSA_Occipital))/sd(brainscSA_Occipital))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_limbic = (brainscSA_Limbic-mean(brainscSA_Limbic))/sd(brainscSA_Limbic))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_cerebellum = (brainscSA_Cerebellum-mean(brainscSA_Cerebellum))/sd(brainscSA_Cerebellum))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_subcortex = (brainscSA_Subcortex-mean(brainscSA_Subcortex))/sd(brainscSA_Subcortex))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(SA_brainstem = (brainscSA_Brainstem-mean(brainscSA_Brainstem))/sd(brainscSA_Brainstem))
df.Sublevel.back0 <- df.Sublevel.back0 %>% mutate(CognitiveCompositeCPM = (brainscNIHcpm-mean(brainscNIHcpm))/sd(brainscNIHcpm))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
df.Sublevel.back0hcp <- df.Sublevel.back0hcp %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))


df.back2shcp <- filter(df.block_FDfilterhcp, Ns == '2')
df.back2shcp <- filter(df.back2shcp,Accs !="NaN") 
df.Sublevel.back2hcp <- df.back2shcp %>% 
  group_by(subs) %>% 
  summarise( Acc = mean(Accs), Age = first(Age), Gender = first(Gender), brainscSA = mean(brainscSA), brainscAveryTask = mean(brainscAveryTask),FDav = mean(FDavs),motion = mean(blockFD), FDmax = mean(FDmaxs),
             brainscSA_Prefrntal = mean(brainscSA_Prefrntal),brainscSA_Motor = mean(brainscSA_Motor),brainscSA_Insula = mean(brainscSA_Insula),brainscSA_Parietal = mean(brainscSA_Parietal),brainscSA_Temporal = mean(brainscSA_Temporal),brainscSA_Occipital = mean(brainscSA_Occipital),brainscSA_Limbic = mean(brainscSA_Limbic), brainscSA_Cerebellum = mean(brainscSA_Cerebellum), brainscSA_Subcortex = mean(brainscSA_Subcortex), brainscSA_Brainstem = mean(brainscSA_Brainstem))
# z-scoring variables so that beta coefficients become comparable
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SustainedAttention = (brainscSA-mean(brainscSA))/sd(brainscSA))

df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_prefrontal = (brainscSA_Prefrntal-mean(brainscSA_Prefrntal))/sd(brainscSA_Prefrntal))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_motor = (brainscSA_Motor-mean(brainscSA_Motor))/sd(brainscSA_Motor))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_insula = (brainscSA_Insula-mean(brainscSA_Insula))/sd(brainscSA_Insula))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_parietal = (brainscSA_Parietal-mean(brainscSA_Parietal))/sd(brainscSA_Parietal))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_temporal = (brainscSA_Temporal-mean(brainscSA_Temporal))/sd(brainscSA_Temporal))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_occipital = (brainscSA_Occipital-mean(brainscSA_Occipital))/sd(brainscSA_Occipital))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_limbic = (brainscSA_Limbic-mean(brainscSA_Limbic))/sd(brainscSA_Limbic))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_cerebellum = (brainscSA_Cerebellum-mean(brainscSA_Cerebellum))/sd(brainscSA_Cerebellum))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_subcortex = (brainscSA_Subcortex-mean(brainscSA_Subcortex))/sd(brainscSA_Subcortex))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(SA_brainstem = (brainscSA_Brainstem-mean(brainscSA_Brainstem))/sd(brainscSA_Brainstem))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(WorkingMemory = (brainscAveryTask-mean(brainscAveryTask))/sd(brainscAveryTask))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(Motion = (FDav-mean(FDav))/sd(FDav))
df.Sublevel.back2hcp <- df.Sublevel.back2hcp %>% mutate(AccZ = (Acc-mean(Acc))/sd(Acc))

#Function to get difference between bootstrapped coefficients in full vs lesioned SA net
## prefrontal
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- summary(lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1))
  lmsa_resampledSA1def <- summary(lm(AccZ ~ Age + Sex + Motion + SA_prefrontal, data=resampled_data1))
  
  lmsa_resampledSA2full <- summary(lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2))
  lmsa_resampledSA2def <- summary(lm(AccZ ~ Age + Gender + Motion + SA_prefrontal, data=resampled_data2))
  
  #deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_prefrontal"])
  
  #deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_prefrontal"])
  
  deltaABCD <-  unname(lmsa_resampledSA1full$r.squared) - unname(lmsa_resampledSA1def$r.squared)
  
  deltaHCP <-   unname(lmsa_resampledSA2full$r.squared) - unname(lmsa_resampledSA2def$r.squared)
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/nrow(boot_deltas))


# motor
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_motor, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_motor, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_motor"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_motor"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/nrow(boot_deltas))

# insula
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_insula, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_insula, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_insula"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_insula"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/nrow(boot_deltas))


# parietal
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_parietal, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_parietal, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_parietal"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_parietal"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/length(boot_coefs))

## temporal
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- summary(lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1))
  lmsa_resampledSA1def <- summary(lm(AccZ ~ Age + Sex + Motion + SA_temporal, data=resampled_data1))
  
  lmsa_resampledSA2full <- summary(lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2))
  lmsa_resampledSA2def <- summary(lm(AccZ ~ Age + Gender + Motion + SA_temporal, data=resampled_data2))
  
  
  #deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_temporal"])
  
  #deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_temporal"])
  
  deltaABCD <-  unname(lmsa_resampledSA1full$r.squared) - unname(lmsa_resampledSA1def$r.squared)
  
  deltaHCP <-   unname(lmsa_resampledSA2full$r.squared) - unname(lmsa_resampledSA2def$r.squared)
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/nrow(boot_deltas))


# occipital
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_occipital, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_occipital, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_occipital"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_occipital"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/length(boot_coefs))


# limbic
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_limbic, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_limbic, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_limbic"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_limbic"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/length(boot_coefs))


# cerebellum
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_cerebellum, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_cerebellum, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_cerebellum"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_cerebellum"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/length(boot_coefs))


## subcortex
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- summary(lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1))
  lmsa_resampledSA1def <- summary(lm(AccZ ~ Age + Sex + Motion + SA_subcortex, data=resampled_data1))
  
  lmsa_resampledSA2full <- summary(lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2))
  lmsa_resampledSA2def <- summary(lm(AccZ ~ Age + Gender + Motion + SA_subcortex, data=resampled_data2))
  
  
  
  #  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_subcortex"])
  
  #  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_subcortex"])
  
  deltaABCD <-  unname(lmsa_resampledSA1full$r.squared) - unname(lmsa_resampledSA1def$r.squared)
  
  deltaHCP <-   unname(lmsa_resampledSA2full$r.squared) - unname(lmsa_resampledSA2def$r.squared)
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/nrow(boot_deltas))


# brainstem
boot_func2 <- function(data1,data2){
  resampled_data1 = data1 %>% 
    sample_frac(1, replace=TRUE)
  
  resampled_data2 = data2 %>% 
    sample_frac(1, replace=TRUE)
  
  lmsa_resampledSA1full <- lm(AccZ ~ Age + Sex + Motion + SustainedAttention, data=resampled_data1)
  lmsa_resampledSA1def <- lm(AccZ ~ Age + Sex + Motion + SA_brainstem, data=resampled_data1)
  
  lmsa_resampledSA2full <- lm(AccZ ~ Age + Gender + Motion + SustainedAttention, data=resampled_data2)
  lmsa_resampledSA2def <- lm(AccZ ~ Age + Gender + Motion + SA_brainstem, data=resampled_data2)
  
  deltaABCD <-  unname(lmsa_resampledSA1full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA1def$coefficients["SA_brainstem"])
  
  deltaHCP <-   unname(lmsa_resampledSA2full$coefficients["SustainedAttention"]) - unname(lmsa_resampledSA2def$coefficients["SA_brainstem"])
  
  return(c(deltaABCD,deltaHCP))
}

#replicating the function 1000 times and saving the output as a vector
boot_coefs <- replicate(1000, boot_func2(df.Sublevel.back0,df.Sublevel.back0hcp))
#histogram of these permuted correlations
boot_deltas <- as.data.frame(t(boot_coefs))
p1 = ggplot(boot_deltas)+ geom_histogram(aes(x = V1),color = 'blue') + geom_histogram(aes(x = V2),color = 'red')
p1
(pval = sum((boot_deltas$V1 - boot_deltas$V2)>0)/length(boot_coefs))

