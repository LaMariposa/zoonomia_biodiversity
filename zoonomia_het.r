setwd("/home/msupple/Desktop/zoonomia_to_sort/heterozygosity/")

#beddir="/home/msupple/zoonomia/zoonomia_biodiversity/het_soh/SOH_files/bed_transfer/"
beddir="/home/msupple/zoonomia/zoonomia_biodiversity/het_soh/SOH_files/SoH/"

#get species file list
#bed_files=grep(".bed.gz", list.files(beddir), value=T)
bed_files=grep(".soh", list.files(beddir), value=T)
data=data.frame(species=substr(bed_files,1,12),
#data=data.frame(species=substr(bed_files,1,6),
                gw_het_mean=numeric(length(bed_files)),
                outbred_het_mean=numeric(length(bed_files)),
                outbred_het_mode=numeric(length(bed_files)))

#loop through species
#for (i in 1:length(bed_files))
for (i in 114:114)
  {
  print(bed_files[i]) 
  
  #read in file
  soh=read.table(paste(beddir,bed_files[i], sep=""))
  colnames(soh)=c("scf","start","end","het","cat")
  
  #genome-wide mean het
  soh$length=soh$end-soh$start
  data$gw_het_mean[i]=weighted.mean(soh$het,soh$length)
  
  #outbred stats subset
  soh_out=subset(soh, cat!=0 & het>0.0)
  #minhet=min(subset(soh, cat!=0)$het)
  #maxhet=max(subset(soh, cat!=0 & het>0.0)$het)
  #if(minhet==0) {print(minhet)}
  #if(maxhet>1) {print(maxhet)}
    
  #outbred mean het 
  data$outbred_het_mean[i]=weighted.mean(soh_out$het,soh_out$length)

  #mode
  histdata=hist(soh_out$het, 
             #   plot=F,
                breaks=500, main=bed_files[i], xlim=c(0,0.01))
  data$outbred_het_mode[i]=histdata$breaks[match(max(histdata$counts), histdata$counts)]
  abline(v=data$outbred_het_mode[i], col="red")
  abline(v=0.0005, col="blue")
  
  #plot ROH (red) and nonROH (black)
  minhetval=0
  maxhetval=2
  ax=pretty(minhetval:maxhetval, n=10000)
  #hist(soh$het, breaks=ax, main=bed_files[i], xlim=c(0,0.01))
  #abline(v=data$outbred_het_mode[i], col="red")
  histdata=hist(subset(soh, cat==0)$het, breaks=ax, main=bed_files[i], 
                #ylim=c(0,5000),
                col="red", xlim=c(0,0.01))
  histdata=hist(subset(soh, cat!=0)$het, breaks=ax, main=bed_files[i],
                add=T)
  abline(v=data$outbred_het_mode[i], col="red")
  abline(v=0.0005, col="blue")
}
write.csv(data,"zoonomia_het.csv")