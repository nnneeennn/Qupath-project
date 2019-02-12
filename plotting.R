#### Qpath output ####
library(ggplot2)
library(gridExtra)


#####Collagen1  #######

setwd("/Volumes/Macintosh HD/DataJoe/VESoutputcollagen1/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/VES_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcVES<-DACmc[complete.cases(DACmc[,32]),]

setwd("/Volumes/Macintosh HD/DataJoe/TMAoutputcollagen1/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/TMA26_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcTMA<-DACmc[complete.cases(DACmc[,32]),]

CombDACmc<-rbind(DACmcVES,DACmcTMA)

# combine the two measurements pf each sample to a mean
Aggdf<-aggregate(CombDACmc$Positive...of.stained.pixels,by=list(CombDACmc$sample_id,CombDACmc$class), FUN=mean)
names(Aggdf)<-c("sample_id","class","meanPosPixels")

Collagen1<-Aggdf

a<-Aggdf[Aggdf$class=="A","meanPosPixels"]
b1<-Aggdf[Aggdf$class=="B1","meanPosPixels"]
b2<-Aggdf[Aggdf$class=="B2","meanPosPixels"]
c<-Aggdf[Aggdf$class=="C","meanPosPixels"]


tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)
name <- Aggdf$class
name[grep("A", Aggdf$class)] <- "DEV"
name[grep("B1", Aggdf$class)] <- "LMS-LIKE"
name[grep("B2", Aggdf$class)] <- "ECM"
name[grep("C", Aggdf$class)] <- "PROLIF"
Aggdf<-cbind(Aggdf,name)
Aggdf$name <- factor(Aggdf$name,levels=c("DEV","LMS-LIKE","ECM","PROLIF"),ordered=T)
#############
plotnew<-ggplot(Aggdf,aes(name,meanPosPixels))+
  geom_boxplot()+
  theme_bw()+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("Collagen 1")+
  geom_point(aes(colour=name),size=5)+
  ylab("Positive % of stained pixels")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))
 
distanz<-(max(Aggdf$meanPosPixels)-min(Aggdf$meanPosPixels))/20

if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$meanPosPixels)+distanz/2+distanz,by=distanz*2,length.out = length(sigv)),
             #label=c(sprintf("p=%.4f",tscores.df$tscores[sigv])),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*")
             ,colour="red",size=10)
}
plotnew



#####Collagen 6 #######

setwd("/Volumes/Macintosh HD/DataJoe/VESoutputcollagen6/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/VES_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcVES<-DACmc[complete.cases(DACmc[,32]),]

setwd("/Volumes/Macintosh HD/DataJoe/TMAoutputcollagen6/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/TMA26_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcTMA<-DACmc[complete.cases(DACmc[,32]),]

CombDACmc<-rbind(DACmcVES,DACmcTMA)




Aggdf<-aggregate(CombDACmc$Positive...of.stained.pixels,by=list(CombDACmc$sample_id,CombDACmc$class), FUN=mean)
names(Aggdf)<-c("sample_id","class","meanPosPixels")

Collagen6<-Aggdf


a<-Aggdf[Aggdf$class=="A","meanPosPixels"]
b1<-Aggdf[Aggdf$class=="B1","meanPosPixels"]
b2<-Aggdf[Aggdf$class=="B2","meanPosPixels"]
c<-Aggdf[Aggdf$class=="C","meanPosPixels"]


tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)

name <- Aggdf$class
name[grep("A", Aggdf$class)] <- "DEV"
name[grep("B1", Aggdf$class)] <- "LMS-LIKE"
name[grep("B2", Aggdf$class)] <- "ECM"
name[grep("C", Aggdf$class)] <- "PROLIF"
Aggdf<-cbind(Aggdf,name)
Aggdf$name <- factor(Aggdf$name,levels=c("DEV","LMS-LIKE","ECM","PROLIF"),ordered=T)

################
plotnew<-ggplot(Aggdf,aes(name,meanPosPixels))+
  geom_boxplot()+
  theme_bw()+
  geom_point(aes(colour=name),size=5)+
  ylab("Positive % of stained pixels")+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("Collagen 6")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))


distanz<-(max(Aggdf$meanPosPixels)-min(Aggdf$meanPosPixels))/20

if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$meanPosPixels)+distanz/2+distanz,by=distanz*2,length.out = length(sigv)),
             #label=c(sprintf("p=%.4f",tscores.df$tscores[sigv])),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*"),colour="red",size=10)
}
plotnew

#####LEM  #######

setwd("/Volumes/Macintosh HD/DataJoe/VESoutputLEM/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/VES_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcVES<-DACmc[complete.cases(DACmc[,32]),]

setwd("/Volumes/Macintosh HD/DataJoe/TMAoutputLEM/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/TMA26_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcTMA<-DACmc[complete.cases(DACmc[,32]),]

CombDACmc<-rbind(DACmcVES,DACmcTMA)


Aggdf<-aggregate(CombDACmc$Positive...of.stained.pixels,by=list(CombDACmc$sample_id,CombDACmc$class), FUN=mean)
names(Aggdf)<-c("sample_id","class","meanPosPixels")

LEM<-Aggdf

a<-Aggdf[Aggdf$class=="A","meanPosPixels"]
b1<-Aggdf[Aggdf$class=="B1","meanPosPixels"]
b2<-Aggdf[Aggdf$class=="B2","meanPosPixels"]
c<-Aggdf[Aggdf$class=="C","meanPosPixels"]


tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)

name <- Aggdf$class
name[grep("A", Aggdf$class)] <- "DEV"
name[grep("B1", Aggdf$class)] <- "LMS-LIKE"
name[grep("B2", Aggdf$class)] <- "ECM"
name[grep("C", Aggdf$class)] <- "PROLIF"
Aggdf<-cbind(Aggdf,name)
Aggdf$name <- factor(Aggdf$name,levels=c("DEV","LMS-LIKE","ECM","PROLIF"),ordered=T)

##############
plotnew<-ggplot(Aggdf,aes(name,meanPosPixels))+
  geom_boxplot()+
  theme_bw()+
  geom_point(aes(colour=name),size=5)+
  ylab("Positive % of stained pixels")+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("LEM")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))


distanz<-(max(Aggdf$meanPosPixels)-min(Aggdf$meanPosPixels))/20

if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$meanPosPixels)+distanz/2+distanz,by=distanz*2,length.out = length(sigv)),
             #label=c(sprintf("p=%.4f",tscores.df$tscores[sigv])),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*"),colour="red",size=10)
}
plotnew

####Fibronectin ##### 

setwd("/Volumes/Macintosh HD/DataJoe/VESoutputfibronectin/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/VES_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcVES<-DACmc[complete.cases(DACmc[,32]),]

setwd("/Volumes/Macintosh HD/DataJoe/TMAoutputfibronectin/forR/")
Annot<-read.table("Annotations",header=T,sep = "\t", stringsAsFactors = F)
Detec<-read.table("Detections",header=T,sep = "\t", stringsAsFactors = F)
Cores<-read.table("TMAdata",header=T,sep = "\t", stringsAsFactors = F)

mapprov<-read.table("/Volumes/Macintosh HD/DataJoe/TMA26_sample_core_map.csv",header=T,sep = "\t", stringsAsFactors = F)
classes<-read.table("/Volumes/Macintosh HD/DataJoe/rna_cluster_sample_map.txt",header=T,sep = "\t", stringsAsFactors = F)

names(Annot)[3]<-"Centroid.X.Annot"
names(Annot)[4]<-"Centroid.Y.Annot"
names(Annot)
names(Detec)

Detec<-Detec[,-c(1)]
names(Detec)[2]<-"Name"

DA<-merge(x = Detec, y = Annot, by = "Name", all.x = TRUE)
names(Cores)[1]<-"Name"
DAC<-merge(x = DA, y = Cores, by = c("Name","Unique.ID"), all.x = TRUE)
DACm<-merge(x = DAC, y = mapprov, by = "Unique.ID", all.x = TRUE)
DACmc<-merge(x = DACm, y = classes, by = "Sample", all.x = TRUE)

DACmcTMA<-DACmc[complete.cases(DACmc[,32]),]

CombDACmc<-rbind(DACmcVES,DACmcTMA)
CombDACmc<-CombDACmc[!CombDACmc$Sample=="Prov 59",] # removes sample 59 which is falsely included 




Aggdf<-aggregate(CombDACmc$Positive...of.stained.pixels,by=list(CombDACmc$sample_id,CombDACmc$class), FUN=mean)
names(Aggdf)<-c("sample_id","class","meanPosPixels")

Fibronectin<-Aggdf

a<-Aggdf[Aggdf$class=="A","meanPosPixels"]
b1<-Aggdf[Aggdf$class=="B1","meanPosPixels"]
b2<-Aggdf[Aggdf$class=="B2","meanPosPixels"]
c<-Aggdf[Aggdf$class=="C","meanPosPixels"]

tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)

name <- Aggdf$class
name[grep("A", Aggdf$class)] <- "DEV"
name[grep("B1", Aggdf$class)] <- "LMS-LIKE"
name[grep("B2", Aggdf$class)] <- "ECM"
name[grep("C", Aggdf$class)] <- "PROLIF"
Aggdf<-cbind(Aggdf,name)
Aggdf$name <- factor(Aggdf$name,levels=c("DEV","LMS-LIKE","ECM","PROLIF"),ordered=T)



#############
plotnew<-ggplot(Aggdf,aes(name,meanPosPixels))+
  geom_boxplot()+
  theme_bw()+
  geom_point(aes(colour=name),size=5)+
  ylab("Positive % of stained pixels")+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("Fibronectin")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))


distanz<-(max(Aggdf$meanPosPixels)-min(Aggdf$meanPosPixels))/20

if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$meanPosPixels)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$meanPosPixels)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$meanPosPixels)+distanz/2+distanz,by=distanz*2,length.out = length(sigv)),
             #label=c(sprintf("p=%.4f",tscores.df$tscores[sigv])),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*"),size=10,
             colour="red")
}
plotnew


##### Nucleus Detection Plots #####

tab<-read.csv2("/Volumes/Macintosh HD/DataJoe/TMA_statistics.csv",sep=",",stringsAsFactors = F)
tabinfo<-tab[,c("sample","class","nucArea","cellsPerArea")]
tabinfo$nucArea<-as.numeric(tabinfo$nucArea)
tabinfo$cellsPerArea<-as.numeric(tabinfo$cellsPerArea)

Aggdf<-aggregate(.~sample+class,tabinfo, FUN=mean)

Cellscounts<-Aggdf
name <- Aggdf$class
name[grep("A", Aggdf$class)] <- "DEV"
name[grep("B1", Aggdf$class)] <- "LMS-LIKE"
name[grep("B2", Aggdf$class)] <- "ECM"
name[grep("C", Aggdf$class)] <- "PROLIF"
Aggdf<-cbind(Aggdf,name)
Aggdf$name <- factor(Aggdf$name,levels=c("DEV","LMS-LIKE","ECM","PROLIF"),ordered=T)

####
a<-Aggdf[Aggdf$class=="A","nucArea"]
b1<-Aggdf[Aggdf$class=="B1","nucArea"]
b2<-Aggdf[Aggdf$class=="B2","nucArea"]
c<-Aggdf[Aggdf$class=="C","nucArea"]


tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)




plotnew<-ggplot(Aggdf,aes(name,nucArea))+
  geom_boxplot()+
  theme_bw()+
  geom_point(aes(colour=name),size=5)+
  ylab("Nuclear Area")+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("Nuclear Area")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))


distanz<-(max(Aggdf$nucArea)-min(Aggdf$nucArea))/20
if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$nucArea)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$nucArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$nucArea)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$nucArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$nucArea)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$nucArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$nucArea)+distanz+distanz/2,by=distanz*2,length.out = length(sigv)),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*"),colour="red",size=10)
}
plotnew



###

a<-Aggdf[Aggdf$class=="A","cellsPerArea"]
b1<-Aggdf[Aggdf$class=="B1","cellsPerArea"]
b2<-Aggdf[Aggdf$class=="B2","cellsPerArea"]
c<-Aggdf[Aggdf$class=="C","cellsPerArea"]


tscores<-t.test(a,b1)$p.value
tscores<-c(tscores,t.test(a,b2)$p.value)
tscores<-c(tscores,t.test(a,c)$p.value)
tscores<-c(tscores,t.test(b1,b2)$p.value)
tscores<-c(tscores,t.test(b1,c)$p.value)
tscores<-c(tscores,t.test(b2,c)$p.value)
names(tscores)<-c("A-B1","A-B2","A-C","B1-B2","B1-C","B2-C")
tscoresig<-tscores<=0.05
coord1<-c(1,1,1,2,2,3)
coord2<-c(2,3,4,3,4,4)
tscores.df<-data.frame(cbind(tscores,tscoresig,coord1,coord2))

sigv<-which(tscores.df$tscoresig==T)

plotnew<-ggplot(Aggdf,aes(name,cellsPerArea))+
  geom_boxplot()+
  theme_bw()+
  geom_point(aes(colour=name),size=5)+
  ylab("Cells Per Area")+
  xlab("RNA Group")+
  theme(legend.title=element_blank(),legend.position="none",plot.title = element_text(hjust = 0.5))+ ggtitle("Cells Per Area")+
  scale_colour_manual(values=c("dodgerblue2","gold1","springgreen3","tomato2"))


distanz<-(max(Aggdf$cellsPerArea)-min(Aggdf$cellsPerArea))/20
if (length(sigv)>=1){
  plotnew<-plotnew+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$cellsPerArea)+distanz,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$cellsPerArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord1[sigv], xend = tscores.df$coord1[sigv], 
             y = seq(from=max(Aggdf$cellsPerArea)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$cellsPerArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("segment", x = tscores.df$coord2[sigv], xend = tscores.df$coord2[sigv], 
             y = seq(from=max(Aggdf$cellsPerArea)+distanz/2,by=distanz*2,length.out = length(sigv)), 
             yend = seq(from=max(Aggdf$cellsPerArea)+distanz,by=distanz*2,length.out = length(sigv)),colour="red")+
    annotate("text",x= (tscores.df$coord1[sigv]+tscores.df$coord2[sigv])/2 ,
             y=seq(from=max(Aggdf$cellsPerArea)+distanz+distanz/2,by=distanz*2,length.out = length(sigv)),
             #label=c(sprintf("p=%.2g",tscores.df$tscores[sigv])),
             label = ifelse(tscores.df$tscores[sigv]<=0.01,"**","*"),
             colour="red",size=10)
}

plotnew


###### correlations #######

names(Collagen1)[3]<-"Collagen1"
names(Collagen6)[3]<-"Collagen6"
names(LEM)[3]<-"LEM"
names(Fibronectin)[3]<-"Fibronectin"



 collecteddf<-merge(x=Collagen1, y=Collagen6, by=c("sample_id","class"),all.x=TRUE)
 collecteddf<-merge(x=collecteddf, y=LEM, by=c("sample_id","class"),all.x=TRUE)
 collecteddf<-merge(x=collecteddf, y=Fibronectin, by=c("sample_id","class"),all.x=TRUE)

 names(Cellscounts)[1]<-"sample_id"
cor(Cellscounts$nucArea,Cellscounts$cellsPerArea)

collecteddf<-merge(x=collecteddf, y=Cellscounts, by=c("sample_id","class"),all = T)
only<-collecteddf[complete.cases(collecteddf),]
cor(only$nucArea,only$Collagen1)
cor(only$nucArea,only$Collagen6)
cor(only$nucArea,only$Fibronectin)
cor(only$nucArea,only$LEM)
cor(only$cellsPerArea,only$Collagen1)
cor(only$cellsPerArea,only$Collagen6)
cor(only$cellsPerArea,only$Fibronectin)
cor(only$cellsPerArea,only$LEM)

library(reshape2)
onlynuc<-melt(only,id.vars=c("sample_id","class","nucArea"))

ggplot(data=onlynuc,aes(x=nucArea,y=value,col=variable))+
  geom_point()+
  geom_smooth(method=lm,se=F)

onlycells<-melt(only,id.vars=c("sample_id","class","cellsPerArea"))

ggplot(data=onlycells,aes(x=cellsPerArea,y=value,col=variable))+
  geom_point()+
  geom_smooth(method=lm,se=F)

ggcol1<-melt(collecteddf,id.vars=c("sample_id","class","Collagen1"))

ggplot(data=ggcol1,aes(x=Collagen1,y=value,col=variable))+
  geom_point()+
  geom_smooth(method=lm,se=F)

ggcol6<-melt(collecteddf,id.vars=c("sample_id","class","Collagen6"))

ggplot(data=ggcol6,aes(x=Collagen6,y=value,col=variable))+
  geom_point()

ggfn<-melt(collecteddf,id.vars=c("sample_id","class","Fibronectin"))

ggplot(data=ggfn,aes(x=Fibronectin,y=value,col=variable))+
  geom_point()

gglem<-melt(collecteddf,id.vars=c("sample_id","class","LEM"))

ggplot(data=gglem,aes(x=LEM,y=value,col=variable))+
  geom_point()

