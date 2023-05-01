setwd("/home/msupple/zoonomia/psmc")

#metafile="zoonomia_240spp.tsv"
metafile="/home/msupple/zoonomia/zoonomia_biodiversity/biodiversity_metadata.csv"
  #when downloading from google doc, need to open and hit save
indir="/home/msupple/zoonomia/psmc/psmc_txt_out/"

# read in metadata
meta=read.table(metafile, header=T, sep=",")
#meta$sp=gsub(" ","_", meta$Species)
#meta$IUCN=meta$IUCN..new.



###############################################################################
# 5 larger order
#set IUCN colors
#palette(c("red","darkgrey","orange","blue","green","brown"))
iucn_cols=c("#377EB8","#4DAF4A","#A65628","#FF7F00","#E41A1C","gray") #LC, NT, VU, EN, CR, DD
#palette(c("#377EB8","#4DAF4A","#A65628","#FF7F00","#E41A1C")) #LC, NT, VU, EN, CR
palette(c("#E41A1C","gray","#FF7F00","#377EB8","#4DAF4A","#A65628"))

# BY ORDER
#pdf("psmc_order_IUCN.pdf")
png("psmc_order_IUCN.png", width=7.6, height=7, units="in", res=1000)
meta$Order2=meta$Order
levels(meta$Order2)=c(levels(meta$Order2),"OTHER")
meta$Order2[meta$Order2!="CARNIVORA" & meta$Order2!="CETARTIODACTYLA" &
      meta$Order!="CHIROPTERA" & meta$Order!="PRIMATES" & meta$Order!="RODENTIA"]="OTHER"
meta$Order2=droplevels(meta$Order2)
ymaxs=c(10,20,30,15,60,30)

orders=levels(meta$Order2)
#orders=c(levels(meta$Order)[2],levels(meta$Order)[7])
par(mfrow=c(3,2), mar=c(4,4,2,0.5))

# loop over IUCN categories
for (i in 1:length(orders))
{
  order=orders[i]
  print(order)
  
  #subset
  metasub=subset(meta, Order2==order)
  
  #set up plot
  plot(c(10000,10000000), c(0,ymaxs[i]),pch="",       
       log='x', ylab="",xlab="", xaxt="n",
       main=order)
  axis(1, tcl=-0.7,
       at=c(10^(4:10)), 
       labels=sapply(c(4:10),function(i) as.expression(bquote(10^.(i)))))
  axis(1, tcl=-0.4, labels=NA,
       at=c(2*10^(4:10),3*10^(4:10),4*10^(4:10),5*10^(4:10),
            6*10^(4:10),7*10^(4:10),8*10^(4:10),9*10^(4:10)))
#  mtext("Years", side=1, line=2.5)
#  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
#        side=2, line=2)
  if (order=="CETARTIODACTYLA")
  {legend("topright",
          legend=c("CR","EN","VU","NT","LC","DD"),
          #legend=c("Critically Endangered","Endangered","Vulnerable","Near threatened","Least Concern","Data deficient"),
          col=iucn_cols,
          #col=c("red","orange","brown","green","blue","darkgrey"),
          #legend=levels(metasub$IUCN),
          #col=as.numeric(as.factor(levels(metasub$IUCN))),
          pch=20)}
  #loop over individuals in subset
  for (j in 1:dim(metasub)[1])
  {
    print(as.character(metasub$species[j]))
    #print(as.character(metasub$Order[j]))
    #print(as.character(metasub$IUCN[j]))
    
    #read in psmc
    if(file.exists(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep="")))
    {
      psmc=read.table(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep=""))
      lines(psmc$V1,psmc$V2, col=metasub$IUCN[j])
    }
  }
}
  mtext("Years", side=1, line=2.5, at=3000)
  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
        side=2, line=31, at=65)
  
dev.off()





#########################################################################
# BY IUCN CATEGORY (6 panels) with captive/wild/dom
pdf("psmc_iucn_wild.pdf")
par(mfrow=c(3,2), mar=c(4,4,2,0.5))
#iucncats=levels(meta$IUCN)
iucncats=c("CR","EN","VU","NT","LC","DD")
#iucncats=c(levels(meta$IUCN)[4],levels(meta$IUCN)[6])

#set wild/cap colors
palette(c("blue","green","red","black"))
#reorder 
meta$wild_status_reseq=factor(meta$wild_status_reseq, level=c("wild","captive","domesticated","unknown"))

# loop over IUCN categories
for (i in 1:length(iucncats))
{
  iucncat=iucncats[i]
  print(iucncat)
  
  #subset
  metasub=subset(meta, IUCN==iucncat)
  
  #set up plot
  plot(c(10000,10000000), c(0,50),pch="",       
       log='x', ylab="",xlab="", xaxt="n",
       main=iucncat
  )
  axis(1, tcl=-0.7,
       at=c(10^(4:10)), 
       labels=sapply(c(4:10),function(i) as.expression(bquote(10^.(i)))))
  axis(1, tcl=-0.4, labels=NA,
       at=c(2*10^(4:10),3*10^(4:10),4*10^(4:10),5*10^(4:10),
            6*10^(4:10),7*10^(4:10),8*10^(4:10),9*10^(4:10)))
#  mtext("Years", side=1, line=2.5)
#  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
#        side=2, line=2)
  if (iucncat=="EN")
  {legend("topright", 
         legend=levels(meta$wild_status_reseq),
         col=c(1:4),
#         col=as.numeric(as.factor(levels(meta$wild_status_reseq))),
         #legend=levels(metasub$IUCN),
         #col=as.numeric(as.factor(levels(metasub$IUCN))),
         pch=20)}
  
  #loop over individuals in subset
  for (j in 1:dim(metasub)[1])
  {
    print(as.character(metasub$species[j]))
    #print(as.character(metasub$Order[j]))
    #print(as.character(metasub$IUCN[j]))
    
    #read in psmc
    if(file.exists(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep="")))
    {
      psmc=read.table(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep=""))
      lines(psmc$V1,psmc$V2, col=as.numeric(metasub$wild_status_reseq[j]))
    }
  }
}
mtext("Years", side=1, line=2.5, at=3000)
mtext(expression(paste("Effective Population Size (x10"^"4",")")),
      side=2, line=28.5, at=110)
dev.off()


#########################################################################
# INDIVIDUAL SAMPLES
pdf("psmc_indiv.pdf")
# loop over individuals
for (i in 1:length(meta$species))
{
  print(as.character(meta$species[i]))
  #print(as.character(meta$Order[i]))
  #print(as.character(meta$IUCN[i]))
  
#check if file exists
if(file.exists(paste(indir,"/",as.character(meta$genus_species[i]),".0.txt", sep="")))
{
  #read in psmc
  psmc=read.table(paste(indir,"/",as.character(meta$genus_species[i]),".0.txt", sep=""))
  psmc[1,1]=psmc[1,1]+0.000001 #fix log 0 problem
  
  #plot
  plot(psmc$V1,psmc$V2, 
       xlim=c(10000,10000000), 
       ylim=c(0,30),
       log='x', type="l",
       main=paste(as.character(meta$species[i]),as.character(meta$Common.Name[i]),
                  as.character(meta$Order[i]),as.character(meta$IUCN[i]), sep=", "),
       ylab="",xlab="", xaxt="n")
  axis(1, tcl=-0.7,
       at=c(10^(4:10)), 
       labels=sapply(c(4:10),function(i) as.expression(bquote(10^.(i)))))
  axis(1, tcl=-0.4, labels=NA,
       at=c(2*10^(4:10),3*10^(4:10),4*10^(4:10),5*10^(4:10),
            6*10^(4:10),7*10^(4:10),8*10^(4:10),9*10^(4:10)))
  mtext("Years", side=1, line=2.5)
  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
        side=2, line=2)
}
}
dev.off()






# BY IUCN CATEGORY
pdf("psmc_iucn.pdf")
iucncats=levels(meta$IUCN)
#iucncats=c(levels(meta$IUCN)[4],levels(meta$IUCN)[6])

# loop over IUCN categories
for (i in 1:length(iucncats))
{
  iucncat=iucncats[i]
  print(iucncat)
  
  #subset
  metasub=subset(meta, IUCN==iucncat)
  
  #set up plot
  plot(c(10000,10000000), c(0,30),pch="",       
       log='x', ylab="",xlab="", xaxt="n",
       main=iucncat
       )
  axis(1, tcl=-0.7,
       at=c(10^(4:10)), 
       labels=sapply(c(4:10),function(i) as.expression(bquote(10^.(i)))))
  axis(1, tcl=-0.4, labels=NA,
       at=c(2*10^(4:10),3*10^(4:10),4*10^(4:10),5*10^(4:10),
            6*10^(4:10),7*10^(4:10),8*10^(4:10),9*10^(4:10)))
  mtext("Years", side=1, line=2.5)
  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
        side=2, line=2)
  
  #loop over individuals in subset
  for (j in 1:dim(metasub)[1])
  {
    print(as.character(metasub$species[j]))
    #print(as.character(metasub$Order[j]))
    #print(as.character(metasub$IUCN[j]))
    
    #read in psmc
    if(file.exists(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep="")))
    {
    psmc=read.table(paste(indir,"/",as.character(metasub$genus_species[j]),".0.txt", sep=""))
    lines(psmc$V1,psmc$V2)
    }
  }
}
dev.off()




# CALC METRICS and plot all
palette(c("#E41A1C","gray","#FF7F00","#377EB8","#4DAF4A","#A65628"))
meta$psmc_min=NA
meta$harm_mean_wt2=NA
pdf("psmc_all.pdf")

  #set up plot
  plot(c(10000,10000000), c(0,30),pch="",       
       log='x', ylab="",xlab="", xaxt="n",
       main="PSMC")
  axis(1, tcl=-0.7,
       at=c(10^(4:10)), 
       labels=sapply(c(4:10),function(i) as.expression(bquote(10^.(i)))))
  axis(1, tcl=-0.4, labels=NA,
       at=c(2*10^(4:10),3*10^(4:10),4*10^(4:10),5*10^(4:10),
            6*10^(4:10),7*10^(4:10),8*10^(4:10),9*10^(4:10)))
  mtext("Years", side=1, line=2.5)
  mtext(expression(paste("Effective Population Size (x10"^"4",")")),
        side=2, line=2)
  
  #loop over individuals in subset
  for (j in 1:dim(meta)[1])
  {
    print(as.character(meta$species[j]))
    #print(as.character(meta$Order[j]))
    #print(as.character(meta$IUCN[j]))
    
    #read in psmc
    if(file.exists(paste(indir,"/",as.character(meta$genus_species[j]),".0.txt", sep="")))
    {
    psmc=read.table(paste(indir,"/",as.character(meta$genus_species[j]),".0.txt", sep=""))
    lines(psmc$V1,psmc$V2, col=meta$IUCN[j])
    
    #calc metrics 
    #psmc min
    meta$psmc_min[j]=min(subset(psmc, V1>10000 & V1<1000000)$V2)
    
    #harmonic mean
    temp_data=subset(psmc, V1>10000 & V1<1000000)[,1:2]
    meta$harm_mean[j]=length(temp_data$V2)/sum(1/temp_data$V2)
    
    #harmonic mean weighted by generation length of interval
    #method 1:  based on older interval
    t2=psmc$V1[-1]
    t2[length(psmc$V1)]=NA
    psmc=cbind(psmc,t2)
    temp_data=psmc[-(length(psmc$V1)),]
    temp_data=subset(temp_data, V1>10000)
    temp_data$wt1=temp_data$t2-temp_data$V1
#    meta$harm_mean_wt1[j]=sum(temp_data$wt1)/sum(temp_data$wt1/temp_data$V2)
    
    #method 2:  based on center of interval
    t1=psmc$V1[-(length(psmc$V1))]
    t1=append(t1,NA, after=0)
    psmc=cbind(psmc,t1)
    temp_data=psmc[-(length(psmc$V1)),]
    temp_data=temp_data[-1,]
    temp_data=subset(temp_data, V1>10000)
    temp_data$wt2=temp_data$t2-temp_data$t1
    meta$harm_mean_wt2[j]=sum(temp_data$wt2)/sum(temp_data$wt2/temp_data$V2)
   
    }
  }
  legend("topright", 
         legend=c("Critically Endangered","Endangered","Vulnerable","Near threatened","Least Concern","Data deficient"),
         col=c("red","orange","brown","green","blue","darkgrey"),
         #legend=levels(metasub$IUCN),
         #col=as.numeric(as.factor(levels(metasub$IUCN))),
         pch=20)

dev.off()

write.csv(meta, "zoonomia_psmc_metrics.csv")

