setwd("/home/msupple/zoonomia/risk_prediction/zoonomia_prediction_nov18/")
infile="/home/msupple/zoonomia/zoonomia_dataset_17nov2021.csv"

library(ape)
library(phylolm)
library(dplyr)
library(tidyr)
library(MASS)
library(car)
library(effects)
library(gridExtra)
library(vioplot)

#########################################################################
##### read in data ######################################################
#########################################################################
### read in tree
ztree=read.tree('/home/msupple/zoonomia/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
#plot.phylo(ztree)

varset=c(52,53,55:58,62,64:69)

zoon=read.table(infile, header=T, sep=",") 

#streamline Order2 for 5 main categories
zoon$Order2=zoon$Order
levels(zoon$Order2)=c(levels(zoon$Order2),"other")
zoon$Order2[zoon$Order2!="CARNIVORA" & zoon$Order2!="CETARTIODACTYLA" &
              zoon$Order!="CHIROPTERA" & zoon$Order!="PRIMATES" & zoon$Order!="RODENTIA"]="other"
zoon$Order2=droplevels(zoon$Order2)

### remove domestic
zoon_nodom=subset(zoon, wild_status_reseq!="domestic" | is.na(wild_status_reseq))
zoon_nodom$wild_status_reseq=droplevels(zoon_nodom$wild_status_reseq)

### rescale variables--z score for all, transform some
#mean_phylop--no tx
z_mean_phylop=(zoon_nodom$mean_phylop-mean(zoon_nodom$mean_phylop, na.rm=T))/sd(zoon_nodom$mean_phylop, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_mean_phylop)
#ppn_conserved--log
z_log_ppn_conserved=(log(zoon_nodom$ppn_conserved)-mean(log(zoon_nodom$ppn_conserved), na.rm=T))/sd(log(zoon_nodom$ppn_conserved), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_ppn_conserved)
#phylop_kurtosis--no tx
z_phylop_kurtosis=(zoon_nodom$phylop_kurtosis-mean(zoon_nodom$phylop_kurtosis, na.rm=T))/sd(zoon_nodom$phylop_kurtosis, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_phylop_kurtosis)
#psmc_min--log
z_log_psmc_min=(log(zoon_nodom$psmc_min)-mean(log(zoon_nodom$psmc_min), na.rm=T))/sd(log(zoon_nodom$psmc_min), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_psmc_min)
#harm_mean_wt2--log
z_log_harm_mean_wt2=(log(zoon_nodom$harm_mean_wt2)-mean(log(zoon_nodom$harm_mean_wt2), na.rm=T))/sd(log(zoon_nodom$harm_mean_wt2), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_harm_mean_wt2)
#gw_het_mean--no tx
z_gw_het_mean=(zoon_nodom$gw_het_mean-mean(zoon_nodom$gw_het_mean, na.rm=T))/sd(zoon_nodom$gw_het_mean, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_gw_het_mean)
#outbred_het_mode--log
z_log_outbred_het_mode=(log(zoon_nodom$outbred_het_mode+0.0000001)-mean(log(zoon_nodom$outbred_het_mode+0.0000001), na.rm=T))/sd(log(zoon_nodom$outbred_het_mode+0.0000001), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_outbred_het_mode)
#froh--sqrt
z_sqrt_froh=(sqrt(zoon_nodom$froh)-mean(sqrt(zoon_nodom$froh), na.rm=T))/sd(sqrt(zoon_nodom$froh), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_sqrt_froh)
#ppn_miss_loeuf--no tx
z_ppn_miss_loeuf=(zoon_nodom$ppn_miss_loeuf-mean(zoon_nodom$ppn_miss_loeuf, na.rm=T))/sd(zoon_nodom$ppn_miss_loeuf, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_ppn_miss_loeuf)
#mean_phylop_missense--no tx
z_mean_phylop_missense=(zoon_nodom$mean_phylop_missense-mean(zoon_nodom$mean_phylop_missense, na.rm=T))/sd(zoon_nodom$mean_phylop_missense, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_mean_phylop_missense)
#mean_phylop_high--log
z_log_mean_phylop_high=(log(zoon_nodom$mean_phylop_high)-mean(log(zoon_nodom$mean_phylop_high), na.rm=T))/sd(log(zoon_nodom$mean_phylop_high), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_mean_phylop_high)
#ppn_miss_conserved--no tx
z_ppn_miss_conserved=(zoon_nodom$ppn_miss_conserved-mean(zoon_nodom$ppn_miss_conserved, na.rm=T))/sd(zoon_nodom$ppn_miss_conserved, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_ppn_miss_conserved)
#ppn_hi_loeuf--log
z_log_ppn_hi_loeuf=(log(zoon_nodom$ppn_hi_loeuf)-mean(log(zoon_nodom$ppn_hi_loeuf), na.rm=T))/sd(log(zoon_nodom$ppn_hi_loeuf), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_ppn_hi_loeuf)
#het_miss_impc_L--log
z_log_het_miss_impc_L=(log(zoon_nodom$het_miss_impc_L)-mean(log(zoon_nodom$het_miss_impc_L), na.rm=T))/sd(log(zoon_nodom$het_miss_impc_L), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_het_miss_impc_L)
#het_miss_impc_V--log
z_log_het_miss_impc_V=(log(zoon_nodom$het_miss_impc_V)-mean(log(zoon_nodom$het_miss_impc_V), na.rm=T))/sd(log(zoon_nodom$het_miss_impc_V), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_het_miss_impc_V)
#het_lof_impc_L--log(x+0.01)
z_log_het_lof_impc_L=(log(zoon_nodom$het_lof_impc_L+0.01)-mean(log(zoon_nodom$het_lof_impc_L+0.01), na.rm=T))/sd(log(zoon_nodom$het_lof_impc_L+0.01), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_het_lof_impc_L)
#het_lof_impc_V--log(x+0.001)
z_log_het_lof_impc_V=(log(zoon_nodom$het_lof_impc_V+0.001)-mean(log(zoon_nodom$het_lof_impc_V+0.001), na.rm=T))/sd(log(zoon_nodom$het_lof_impc_V+0.001), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_het_lof_impc_V)
#hom_miss_impc_L--log
z_log_hom_miss_impc_L=(log(zoon_nodom$hom_miss_impc_L)-mean(log(zoon_nodom$hom_miss_impc_L), na.rm=T))/sd(log(zoon_nodom$hom_miss_impc_L), na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_log_hom_miss_impc_L)
#hom_miss_impc_V--no tx
z_hom_miss_impc_V=(zoon_nodom$hom_miss_impc_V-mean(zoon_nodom$hom_miss_impc_V, na.rm=T))/sd(zoon_nodom$hom_miss_impc_V, na.rm=T)
zoon_nodom=cbind(zoon_nodom,z_hom_miss_impc_V)

### remove DD species
zoon_nodd=subset(zoon_nodom, IUCN!="DD")
zoon_nodd$IUCN=droplevels(zoon_nodd$IUCN)

### clean up dataset
row.names(zoon_nodd)=zoon_nodd$Sp
zoon_nodd$IUCN=factor(zoon_nodd$IUCN,levels=c('LC','NT','VU','EN','CR'), ordered=T)
zoon_nodd$threatened=factor(zoon_nodd$threatened)

#remove "other" Orders
zoon_order=subset(zoon_nodd, Order2!="other")
zoon_order$Order2=droplevels(zoon_order$Order2)



#########################################################################
##### LOGISTIC REGRESSION / ODDS RATIO ##################################
#########################################################################
##### model selection for logistic regression ###########################
### covariates ##########################################################
## individually
covm=phyloglm(threatened~wild_status_reseq,data=zoon_nodd,phy=ztree)
summary(covm)
covm=phyloglm(threatened~diet_hoc,data=zoon_nodd,phy=ztree)
summary(covm)

### predictor variables #################################################
## individually
#visualize for rln with order
for (i in varset)
{
  do.call("vioplot", list(formula=as.formula(paste0(colnames(zoon_nodd)[i], " ~ Order")), data=zoon_nodd, col="grey90"))
  points(jitter(as.numeric(zoon_nodd$Order)), zoon_nodd[,i], col=zoon_nodd$IUCN, pch=20)
  #vioplot(mean_phylop ~ Order, data=zoon_nodd, col="grey90")
  #points(jitter(as.numeric(zoon_nodd$Order)), zoon_nodd$mean_phylop, col=zoon_nodd$IUCN, pch=20)
}

for (i in varset)
{
  #print(colnames(zoon_nodd)[i])
  predm=do.call("phyloglm", list(formula=as.formula(paste0("threatened ~ ", colnames(zoon_nodd)[i])),
                                 data=zoon_nodd,phy=ztree))
  #summary(predm)   #doesn't print out properly
  print(2*pnorm(-abs(predm$coefficients[2]/predm$sd[2])))
}

## examine correlations 
x=cor(zoon_nodd[varset],zoon_nodd[varset], use="na.or.complete", method="pearson") %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(abs(value))) %>%
  group_by(value) %>%
  filter(row_number()==1)
x[1:20,]


##### model odds ratio ##################################################
#fulllr=phyloglm(threatened ~ diet_hoc + wild_status_reseq
#                + z_log_harm_mean_wt2 + z_gw_het_mean + z_log_outbred_het_mode + z_sqrt_froh
#                + z_log_ppn_conserved + z_phylop_kurtosis + z_ppn_miss_conserved
#                + z_log_het_miss_impc_L + z_log_het_miss_impc_V + z_log_het_lof_impc_L  
#                + z_log_het_lof_impc_V + z_log_hom_miss_impc_L + z_hom_miss_impc_V, 
#                data=zoon_nodd,phy=ztree)
fulllr=phyloglm(threatened ~ diet_hoc + wild_status_reseq
              + z_log_harm_mean_wt2 + z_log_outbred_het_mode + z_sqrt_froh 
              + z_hom_miss_impc_V + z_log_ppn_conserved + z_log_het_lof_impc_L,
              data=zoon_nodd,phy=ztree)
summary(fulllr)



#########################################################################
##### ORDINAL REGRESSION / PROPORTIONAL ODDS ############################
#########################################################################
##### model selection for ordinal regression ############################
### covariates ##########################################################
## individually
covm=polr(IUCN~Order2,data=zoon_order, Hess=T)
Anova(covm)
covm=polr(IUCN~wild_status_reseq, data=zoon_order, Hess=T)
Anova(covm)
covm=polr(IUCN~diet_hoc, data=zoon_order, Hess=T)
Anova(covm)

### predictor variables #################################################
## individually
for (i in varset)
{
  #print(colnames(zoon_order)[i])
  predm=do.call("polr", list(formula=as.formula(paste0("IUCN ~ ", colnames(zoon_order)[i])),
                                 data=zoon_order, Hess=T))
  #summary(predm)   #doesn't print out properly
  print(Anova(predm))
}

## examine correlations 
x=cor(zoon_order[varset],zoon_order[varset], use="na.or.complete", method="pearson") %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(abs(value))) %>%
  group_by(value) %>%
  filter(row_number()==1)
x[1:20,]


##### model ordinal regression ##########################################
#fullor=polr(IUCN ~ Order2 + diet_hoc + wild_status_reseq
#            + z_log_harm_mean_wt2 + z_gw_het_mean + z_log_outbred_het_mode + z_sqrt_froh
#            + z_log_ppn_conserved + z_phylop_kurtosis + z_ppn_miss_conserved
#            + z_log_het_miss_impc_L + z_log_het_miss_impc_V + z_log_het_lof_impc_L  
#            + z_log_het_lof_impc_V + z_log_hom_miss_impc_L + z_hom_miss_impc_V,
#            data=zoon_order, Hess=T)
fullor=polr(IUCN ~ Order2 + diet_hoc + wild_status_reseq
            + z_log_harm_mean_wt2
            + z_log_ppn_conserved + z_log_het_miss_impc_V + z_hom_miss_impc_V,
            data=zoon_order, Hess=T)
Anova(fullor)
#summary(fullor)
#S(fullor)

##### predict on DD species #############################################
subset(zoon_nodom, IUCN=="DD")$Sp
predictDD=predict(fullor, subset(zoon_nodom, IUCN=="DD"), type="probs")
predictDD
predictDD[,1]                #LC
rowSums(predictDD[,2:5])     #>LC
predictDD=predict(fullor, subset(zoon_nodom, IUCN=="DD"))
predictDD


##### cross validation #################################################
nreps=100
cverrs=rep(NA,nreps)
cverrs_iucn=rep(NA,nreps)
samplesize=0.60*nrow(zoon_order)

#set.seed(19802658)
for (i in 1:nreps)
{
  index=sample(seq_len(nrow(zoon_order)), size=samplesize)
  datatrain=zoon_order[index,]
  datatest=zoon_order[-index,]

  ### build model on training
  model=polr(IUCN ~ Order2 + diet_hoc + wild_status_reseq
          + z_log_harm_mean_wt2
          + z_log_ppn_conserved + z_log_het_miss_impc_V + z_hom_miss_impc_V,
          data=datatrain, Hess=T)
  #summary(model)
  ### test model
  predictIUCN=predict(model, datatest)
  table(datatest$IUCN, predictIUCN)
  cverrs_iucn[i]=mean(as.character(datatest$IUCN)!=as.character(predictIUCN), na.rm=T)
  risk_actual=datatest$IUCN
  levels(risk_actual)=c(levels(risk_actual), "not","risk")
  risk_actual[risk_actual=="LC"]="not"
  risk_actual[risk_actual=="NT"]="risk"
  risk_actual[risk_actual=="VU"]="risk"
  risk_actual[risk_actual=="EN"]="risk"
  risk_actual[risk_actual=="CR"]="risk"
  risk_actual=droplevels(risk_actual)
  risk_predict=predictIUCN
  levels(risk_predict)=c(levels(risk_predict), "not","risk")
  risk_predict[risk_predict=="LC"]="not"
  risk_predict[risk_predict=="NT"]="risk"
  risk_predict[risk_predict=="VU"]="risk"
  risk_predict[risk_predict=="EN"]="risk"
  risk_predict[risk_predict=="CR"]="risk"
  risk_predict=droplevels(risk_predict)
  table(risk_actual,risk_predict)
  cverrs[i]=mean(as.character(risk_actual)!=as.character(risk_predict), na.rm=T)
}
mean(cverrs)
mean(cverrs_iucn)
