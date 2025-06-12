##this code generates a heatmap of enrichment of GO terms of interest (TOI) and their associations with viral groups
##last updated LEF 6/12/25

##load in packages
library(pheatmap)
library(RColorBrewer)
library(wesanderson)
library(dplyr)
library(tidyr)
library(tidyverse)

##load in terms ##
  TOI<- read.csv('TOI.csv')

##collect enrichment data, selecting delta.ranks and term names for each file##
  Bing=read.delim("MWU_BP_new_Bingvirus_GOMWU.csv",sep=' ')
  Bing[Bing$p.adj > 0.1, 1] <- 0
  Bing=Bing[,c(1,6)]
  colnames(Bing)=c("Bingvirus","name")
  
  BL3=read.delim("MWU_BP_new_BrochothrixphageBL3_GOMWU.csv",sep=' ')
  BL3[BL3$p.adj > 0.1, 1] <- 0
  BL3=BL3[,c(1,6)]
  colnames(BL3)=c("Brochothrix phage BL3","name")
  
  Burro=read.delim("MWU_BP_new_Burrovirus_GOMWU.csv",sep=' ')
  Burro[Burro$p.adj > 0.1, 1] <- 0
  Burro=Burro[,c(1,6)]
  colnames(Burro)=c("Burrovirus","name")
  
  Carme=read.delim("MWU_BP_new_Carmenvirus_GOMWU.csv",sep=' ')
  Carme[Carme$p.adj > 0.1, 1] <- 0
  Carme=Carme[,c(1,6)]
  colnames(Carme)=c("Carmenvirus","name")
  
  CT19406C=read.delim("MWU_BP_new_ClostridiumphagephiCT19406C_GOMWU.csv",sep=' ')
  CT19406C[CT19406C$p.adj > 0.1, 1] <- 0
  CT19406C=CT19406C[,c(1,6)]
  colnames(CT19406C)=c("Clostridium phage phiCT19406C","name")
  
  Decur=read.delim("MWU_BP_new_Decurrovirus_GOMWU.csv",sep=' ')
  Decur[Decur$p.adj > 0.1, 1] <- 0
  Decur=Decur[,c(1,6)]
  colnames(Decur)=c("Decurrovirus","name")
  
  Godon=read.delim("MWU_BP_new_Godonkavirus_GOMWU.csv",sep=' ')
  Godon[Godon$p.adj > 0.1, 1] <- 0
  Godon=Godon[,c(1,6)]
  colnames(Godon)=c("Godonkavirus","name")
  
  Halcy=read.delim("MWU_BP_new_Halcyonevirus_GOMWU.csv",sep=' ')
  Halcy[Halcy$p.adj > 0.1, 1] <- 0
  Halcy=Halcy[,c(1,6)]
  colnames(Halcy)=c("Halcyonevirus","name")
  
  K02=read.delim("MWU_BP_new_KlebsiellaphagephiKO2_GOMWU.csv",sep=' ')
  K02[K02$p.adj > 0.1, 1] <- 0
  K02=K02[,c(1,6)]
  colnames(K02)=c("Klebsiella phage phiKO2","name")
  
  Lacu=read.delim("MWU_BP_new_Lacusarxvirus_GOMWU.csv",sep=' ')
  Lacu[Lacu$p.adj > 0.1, 1] <- 0
  Lacu=Lacu[,c(1,6)]
  colnames(Lacu)=c("Lacusarxvirus","name")
  
  Light=read.delim("MWU_BP_new_Lightbulbvirus_GOMWU.csv",sep=' ')
  Light[Light$p.adj > 0.1, 1] <- 0
  Light=Light[,c(1,6)]
  colnames(Light)=c("Lightbulbvirus","name")
  
  Nona=read.delim("MWU_BP_new_Nonagvirus_GOMWU.csv",sep=' ')
  Nona[Nona$p.adj > 0.1, 1] <- 0
  Nona=Nona[,c(1,6)]
  colnames(Nona)=c("Nonagvirus","name")
  
  H103=read.delim("MWU_BP_new_PseudoalteromonasphageH103_GOMWU.csv",sep=' ')
  H103[H103$p.adj > 0.1, 1] <- 0
  H103=H103[,c(1,6)]
  colnames(H103)=c("Pseudoalteromonas phage H103","name")
  
  JBD44=read.delim("MWU_BP_new_PseudomonasphageJBD44_GOMWU.csv",sep=' ')
  JBD44[JBD44$p.adj > 0.1, 1] <- 0
  JBD44=JBD44[,c(1,6)]
  colnames(JBD44)=c("Pseudomonas phage JBD44","name")
  
  Schit=read.delim("MWU_BP_new_Schitoviridae_GOMWU.csv",sep=' ')
  Schit[Schit$p.adj > 0.1, 1] <- 0
  Schit=Schit[,c(1,6)]
  colnames(Schit)=c("Schitoviridae","name")
  
  Spiz=read.delim("MWU_BP_new_Spizizenvirus_GOMWU.csv",sep=' ')
  Spiz[Spiz$p.adj > 0.1, 1] <- 0
  Spiz=Spiz[,c(1,6)]
  colnames(Spiz)=c("Spizizenvirus","name")
  
  ARI0746=read.delim("MWU_BP_new_StreptococcusphagephiARI0746_GOMWU.csv",sep=' ')
  ARI0746[ARI0746$p.adj > 0.1, 1] <- 0
  ARI0746=ARI0746[,c(1,6)]
  colnames(ARI0746)=c("Streptococcus phage phiARI0746","name")
  
  NJ2=read.delim("MWU_BP_new_StreptococcusphagephiNJ2_GOMWU.csv",sep=' ')
  NJ2[NJ2$p.adj > 0.1, 1] <- 0
  NJ2=NJ2[,c(1,6)]
  colnames(NJ2)=c("Streptococcus phage phiNJ2","name")
  
  Ute=read.delim("MWU_BP_new_Uetakevirus_GOMWU.csv",sep=' ')
  Ute[Ute$p.adj > 0.1, 1] <- 0
  Ute=Ute[,c(1,6)]
  colnames(Ute)=c("Uetakevirus","name")
  
  Verto=read.delim("MWU_BP_new_Vertoviridae_GOMWU.csv",sep=' ')
  Verto[Verto$p.adj > 0.1, 1] <- 0
  Verto=Verto[,c(1,6)]
  colnames(Verto)=c("Vertoviridae","name")
  
  Vivid=read.delim("MWU_BP_new_Vividuovirus_GOMWU.csv",sep=' ')
  Vivid[Vivid$p.adj > 0.1, 1] <- 0
  Vivid=Vivid[,c(1,6)]
  colnames(Vivid)=c("Vividuovirus","name")
  
  Zobel=read.delim("MWU_BP_new_Zobellviridae_GOMWU.csv",sep=' ')
  Zobel[Zobel$p.adj > 0.1, 1] <- 0
  Zobel=Zobel[,c(1,6)]
  colnames(Zobel)=c("Zobellviridae","name")

  
##merge collected data with our TOI frame to just pull data for our terms of interest##
  test=merge(TOI, Bing, by="name")
  test=merge(test,BL3, by="name")
  test=merge(test,Burro, by="name")  
  test=merge(test,Carme, by="name")
  test=merge(test,Decur, by="name")
  test=merge(test,Godon, by="name")
  test=merge(test,Halcy, by="name")
  test=merge(test,K02, by="name")
  test=merge(test,Lacu, by="name")
  test=merge(test,Light, by="name")
  test=merge(test,Nona, by="name")
  test=merge(test,H103, by="name")
  test=merge(test,JBD44, by="name")
  test=merge(test,Schit, by="name")
  test=merge(test,Spiz, by="name")
  test=merge(test,ARI0746, by="name")
  test=merge(test,NJ2, by="name")
  test=merge(test,Ute, by="name")
  test=merge(test,Vivid, by="name")
  test=merge(test,Verto, by="name")
  test=merge(test,Zobel, by="name")
  
  ##fix data frame for heatmapping##
  test=test %>% remove_rownames %>% column_to_rownames(var="name")
  #test=as.data.frame(t(test))

  
##and make the heatmap##
  
  pheatmap(test, show_rownames=T, color = colorRampPalette(c("#5a7bac","#dee6f6","#bd004a"))(100))  
  