setwd("/home/msupple/zoonomia/risk_prediction/zoonomia_prediction_nov18/")

infile="/home/msupple/zoonomia/zoonomia_dataset_17nov2021.csv"

#PCA regression manual (rather than pls package, including order, and selecting PCs)
#https://rpubs.com/esobolewska/pcr-step-by-step

library(pls)
library(corrplot)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
#palette(c("blue","forestgreen","brown","orange","red","black")) #LC->CR
palette(c("#377EB8","#4DAF4A","#A65628","#FF7F00","#E41A1C","gray")) #LC->CR
iucncols=c("#377EB8","#4DAF4A","#A65628","#FF7F00","#E41A1C","gray")
#########################################################################
##### read in data ######################################################
#########################################################################
zoon=read.table(infile, header=T, sep=",") 
row.names(zoon)=zoon$Sp

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


### clean up dataset
row.names(zoon_nodom)=zoon_nodom$Sp
zoon_nodom$IUCN=factor(zoon_nodom$IUCN,levels=c('LC','NT','VU','EN','CR','DD'), ordered=T)

#remove "other" Orders
zoon_order=subset(zoon_nodom, Order2!="other")
zoon_order$Order2=droplevels(zoon_order$Order2)


##########################################################################
##### PCA regression of full dataset #####################################
##########################################################################
##### prep data ###################################################################
#get subset of relevant data
indata=zoon_nodom[,c(50,9,10,52,53,55:58,62,64:69)]
#remove missing
indata=indata[complete.cases(indata[4:16]),]

res=cor(indata[,4:16], method="pearson")
corrplot::corrplot(res, method= "color", order = "hclust", tl.pos = 'n')

##### PCA ######################################################################
pcaout=prcomp(indata[,4:16], center=T, scale=T)
summary(pcaout)

#loadings
loads=pcaout$rotation
write.csv(loads, file="pc_loadings.csv")

#default
plot(summary(pcaout)$importance[3,])  #5 PCs gets 84% cum prop
fviz_pca_ind(pcaout,axes = c(1, 3), geom=c("point"), col.ind=indata$IUCN)
fviz_pca_var(pcaout,axes = c(1, 4), repel=T)
fviz_pca_biplot(pcaout,axes = c(1, 3), geom=c("point"), col.ind=indata$IUCN)

#pretty plots
pcaout_clean=pcaout
rownames(pcaout_clean$rotation)=c("ppn_conserved","phylop_kurtosis","harm_mean_Ne","gw_het_mean",
                         "outbred_het_mode","froh","ppn_miss_conserved","het_miss_impc_L",
                         "het_miss_impc_V","het_lof_impc_L","het_lof_impc_V",
                         "hom_miss_impc_L","hom_miss_impc_V")


fviz_pca_biplot(pcaout_clean,axes = c(1, 3), geom=c("point"), repel=T, label="none",
                arrowsize=1,
                col.ind=NA, alpha.ind=0, col.var="black") +
  #col.ind=indata$IUCN, habillage=indata$IUCN, addEllipses = T,alpha.ind=0.1,
#  geom_point(color=as.numeric(indata$IUCN), alpha=0.4, size=3) +
  geom_point(color=as.numeric(indata$IUCN), alpha=1, size=2) +
  labs(title="", x="PC 1 (35%)", y="PC 3 (13%)", color="IUCN") +
  scale_shape_manual(values=c(20,20,20,20,20,20)) +
  scale_color_manual(name="IUCN", breaks=c("abc","def"), values=c("abc"="red", "def"="blue"))
#legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", title="IUCN", legend = levels(indata$IUCN), pch=20, cex=2,
       bty='n', col = alpha(iucncols, 1))

  


plot(pcaout$x[,1], pcaout$x[,3], pch=20, col=indata$IUCN, cex=1.2,
     xlab="PC 1 (34.5%)", ylab="PC 2 (13.3%)")
legend("bottomright", levels(indata$IUCN), pch=20, cex=1.4, col=iucncols)

oitpcs=cbind(indata$Order2,indata$IUCN,indata$threatened,as.data.frame(pcaout$x[,1:5]))
colnames(oitpcs)[1:3]=c("Order2","IUCN","threatened")
#plot(oitpcs$threatened,oitpcs$PC1)  #give binary response, not very instructive

##### regression #################################################################
### remove DD species
zoon_nodd=subset(oitpcs, IUCN!="DD")
zoon_nodd$IUCN=droplevels(zoon_nodd$IUCN)

### model
fullmod=lm(threatened ~ PC1 + PC2 + PC3 + PC4 + PC5, data=zoon_nodd)
summary(fullmod)


##### predict DD species #####################################
rownames(subset(oitpcs, IUCN=="DD"))
predictDD=predict.lm(fullmod, subset(oitpcs, IUCN=="DD"))
predictDD


##### cross validation #################################################
#fullmod=lm(threatened ~ PC1 + PC2 + PC3 + PC4 + PC5, data=zoon_nodd)
nreps=100
cverrs=rep(NA,nreps)
samplesize=0.60*nrow(zoon_nodd)

### partition training and validation sets
#set.seed(19802658)
for (i in 1:nreps)
{
  index=sample(seq_len(nrow(zoon_nodd)), size=samplesize)
  datatrain=zoon_nodd[index,]
  datatest=zoon_nodd[-index,]
  
  ### build model on training
  model=lm(threatened ~ PC1 + PC2 + PC3 + PC4 + PC5, data=datatrain)
  #summary(model)
  ### test model
  predictIUCN=predict(model, datatest)
  #convert from numeric to integer (<0.5->0; >=0.5->1 threat)
  risk_predict=as.factor(predictIUCN)
  levels(risk_predict)=c(levels(risk_predict), "not","risk")
  risk_predict[predictIUCN<0.5]="not"
  risk_predict[predictIUCN>=0.5]="risk"
  risk_predict=droplevels(risk_predict)
  #convert IUCN to not/risk
  risk_actual=datatest$IUCN
  levels(risk_actual)=c(levels(risk_actual), "not","risk")
  risk_actual[risk_actual=="LC"]="not"
  risk_actual[risk_actual=="NT"]="risk"
  risk_actual[risk_actual=="VU"]="risk"
  risk_actual[risk_actual=="EN"]="risk"
  risk_actual[risk_actual=="CR"]="risk"
  risk_actual=droplevels(risk_actual)
  
  table(risk_actual,risk_predict)
  cverrs[i]=mean(as.character(risk_actual)!=as.character(risk_predict), na.rm=T)
}
mean(cverrs)




##########################################################################
##### PCA regression of five orders ######################################
##########################################################################
##### prep data ###################################################################
#get subset of relevant data
indata=zoon_order[,c(50,9,10,52,53,55:58,62,64:69)]
#remove missing
indata=indata[complete.cases(indata),]

res=cor(indata[,4:16], method="pearson")
corrplot::corrplot(res, method= "color", order = "hclust", tl.pos = 'n')


##### PCA ######################################################################
pcaout=prcomp(indata[,4:16], center=T, scale=T)
summary(pcaout)

plot(summary(pcaout)$importance[3,])  #5 PCs gets 85% cum prop
fviz_pca_var(pcaout,axes = c(1, 3))
plot(pcaout$x[,1], pcaout$x[,3], pch=20, col=indata$IUCN, cex=1.2)

oitpcs=cbind(indata$Order2,indata$IUCN,indata$threatened,as.data.frame(pcaout$x[,1:5]))
colnames(oitpcs)[1:3]=c("Order2","IUCN","threatened")
#plot(oitpcs$threatened,oitpcs$PC1)  #give binary response, not very instructive


##### regression #################################################################
### remove DD species
zoon_nodd=subset(oitpcs, IUCN!="DD")
zoon_nodd$IUCN=droplevels(zoon_nodd$IUCN)

### model
fullmod=lm(threatened ~ Order2 + PC1 + PC2 + PC3 + PC4 + PC5, data=zoon_nodd)
summary(fullmod)


##### predict DD species #####################################
rownames(subset(oitpcs, IUCN=="DD"))
predictDD=predict.lm(fullmod, subset(oitpcs, IUCN=="DD"))
predictDD


##### cross validation #################################################
#fullmod=lm(threatened ~ Order2 + PC1 + PC2 + PC3 + PC4 + PC5, data=zoon_nodd)
nreps=100
cverrs=rep(NA,nreps)
samplesize=0.60*nrow(zoon_nodd)

### partition training and validation sets
#set.seed(19802658)
for (i in 1:nreps)
{
  index=sample(seq_len(nrow(zoon_nodd)), size=samplesize)
  datatrain=zoon_nodd[index,]
  datatest=zoon_nodd[-index,]

  ### build model on training
  model=lm(threatened ~ Order2 + PC1 + PC2 + PC3 + PC4 + PC5, data=datatrain)
  #summary(model)
  ### test model
  predictIUCN=predict(model, datatest)
  #convert from numeric to integer (<0.5->0; >=0.5->1 threat)
  risk_predict=as.factor(predictIUCN)
  levels(risk_predict)=c(levels(risk_predict), "not","risk")
  risk_predict[predictIUCN<0.5]="not"
  risk_predict[predictIUCN>=0.5]="risk"
  risk_predict=droplevels(risk_predict)
  #convert IUCN to not/risk
  risk_actual=datatest$IUCN
  levels(risk_actual)=c(levels(risk_actual), "not","risk")
  risk_actual[risk_actual=="LC"]="not"
  risk_actual[risk_actual=="NT"]="risk"
  risk_actual[risk_actual=="VU"]="risk"
  risk_actual[risk_actual=="EN"]="risk"
  risk_actual[risk_actual=="CR"]="risk"
  risk_actual=droplevels(risk_actual)

  table(risk_actual,risk_predict)
  cverrs[i]=mean(as.character(risk_actual)!=as.character(risk_predict), na.rm=T)
}
mean(cverrs)


###################################################################
##### phyloreg ###################################################
library(phylolm)
### read in tree
ztree=read.tree('/home/msupple/zoonomia/risk_prediction/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
#plot.phylo(ztree)
### model
fulllr=phyloglm(threatened ~ PC1 + PC2 + PC3 + PC4 + PC5,
                data=zoon_nodd,phy=ztree)
summary(fulllr)

pcthreat=phylolm(PC4 ~ threatened, data=zoon_nodd,phy=ztree)
summary(pcthreat)
