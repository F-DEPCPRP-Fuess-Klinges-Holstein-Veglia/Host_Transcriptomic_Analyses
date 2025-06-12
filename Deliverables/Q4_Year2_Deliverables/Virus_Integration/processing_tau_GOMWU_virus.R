##this code processes our tau value results matrix to generate input for gene ontology enrichment analysis with GOMWU
##last updated LEF 6/12/25

##load in your necessary packages
library(data.table)
library(tidyr)
library(janitor)
library(tibble)
##load in the read data to get a list of all orthologs (ncessary for enrichment##
  genenames = read.csv("normalized_reads_corr.csv", header=TRUE)
  genenames=t(genenames)
  genenames=genenames %>%
    row_to_names(row_number = 1)
  genenames=as.data.frame(genenames)
  genenames <- tibble::rownames_to_column(genenames, "sample_R")
  genenames=genenames[,c(1:2)]
  colnames(genenames)[1] <- "gene"


##and parse each results file individully to generate a GOMWU input file##
##I've fully annotated the first iteration only, but it's just repeated code tweaked as needing for each data frame##

  ##Bing
    ##read in tau values##
    Bing = read.csv('Corr_results_tau_Bingvirus.csv', check.names = FALSE)
    Bing=setDT(Bing, keep.rownames = TRUE)[]
    dim(Bing)
    ##keep just the column with ortho IDs and taus of associaitons with viral group##
    Bing=Bing[,c(1,385)]
    ##rename columns
    colnames(Bing) <- c("gene", "tau")
    ##write out the tau values for supplemental file##
    write.csv(Bing, file="Bingvirus_taus.csv", row.names = FALSE)
    ##merge with list of all orthologs (necessary for GOMWU)
    Merged_Bing = merge(Bing, genenames, by = "gene", all.y = TRUE)
    ##replace NAs (i.e. orthos that weren't in your top 5%) with 0, following protocol from Fuess et al. mBio##
    Merged_Bing[is.na(Merged_Bing)] <- 0
    dim(Merged_Bing)
    ##select just the ortho ID and tau value columns
    Merged_Bing=Merged_Bing[,c(1:2)]
    ##write it all out##
    write.csv(Merged_Bing, file="Bingvirus_GOMWU.csv", row.names = FALSE)
  
  
  ##BL3
  BL3 = read.csv('Corr_results_tau_BrochothrixphageBL3.csv', check.names = FALSE)
  BL3=setDT(BL3, keep.rownames = TRUE)[]
  dim(BL3)
  BL3=BL3[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(BL3) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(BL3, file="BrochothrixphageBL3__taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_BL3 = merge(BL3, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_BL3[is.na(Merged_BL3)] <- 0
  dim(Merged_BL3)
  Merged_BL3=Merged_BL3[,c(1:2)]
  ##writeitout##
  write.csv(Merged_BL3, file="BrochothrixphageBL3_GOMWU.csv", row.names = FALSE)
  
  ##Burro
  Burro = read.csv('Corr_results_tau_Burrovirus.csv', check.names = FALSE)
  Burro=setDT(Burro, keep.rownames = TRUE)[]
  dim(Burro)
  Burro=Burro[,c(1,338)]
  ##make the GOMWU matrix##
  colnames(Burro) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Burro, file="Burrovirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Burro = merge(Burro, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Burro[is.na(Merged_Burro)] <- 0
  dim(Merged_Burro)
  Merged_Burro=Merged_Burro[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Burro, file="Burrovirus_GOMWU.csv", row.names = FALSE)
  
  ##Carm
  Carm = read.csv('Corr_results_tau_Carmenvirus.csv', check.names = FALSE)
  Carm=setDT(Carm, keep.rownames = TRUE)[]
  dim(Carm)
  Carm=Carm[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Carm) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Carm, file="Carmenvirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Carm = merge(Carm, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Carm[is.na(Merged_Carm)] <- 0
  dim(Merged_Carm)
  Merged_Carm=Merged_Carm[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Carm, file="Carmenvirus_GOMWU.csv", row.names = FALSE)
  
  ##CT19406C
  CT19406C = read.csv('Corr_results_tau_ClostridiumphagephiCT19406C.csv', check.names = FALSE)
  CT19406C=setDT(CT19406C, keep.rownames = TRUE)[]
  dim(CT19406C)
  CT19406C=CT19406C[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(CT19406C) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(CT19406C, file="ClostridiumphagephiCT19406C_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_CT19406C = merge(CT19406C, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_CT19406C[is.na(Merged_CT19406C)] <- 0
  dim(Merged_CT19406C)
  Merged_CT19406C=Merged_CT19406C[,c(1:2)]
  ##writeitout##
  write.csv(Merged_CT19406C, file="ClostridiumphagephiCT19406C_GOMWU.csv", row.names = FALSE)
  
  
 
  ##Decu
  Decu = read.csv('Corr_results_tau_Decurrovirus.csv', check.names = FALSE)
  Decu=setDT(Decu, keep.rownames = TRUE)[]
  dim(Decu)
  Decu=Decu[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Decu) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Decu, file="Decurrovirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Decu = merge(Decu, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Decu[is.na(Merged_Decu)] <- 0
  dim(Merged_Decu)
  Merged_Decu=Merged_Decu[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Decu, file="Decurrovirus_GOMWU.csv", row.names = FALSE)
  
  
  ##Godo
  Godo = read.csv('Corr_results_tau_Godonkavirus.csv', check.names = FALSE)
  Godo=setDT(Godo, keep.rownames = TRUE)[]
  dim(Godo)
  Godo=Godo[,c(1,330)]
  ##make the GOMWU matrix##
  colnames(Godo) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Godo, file="Godonkavirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Godo = merge(Godo, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Godo[is.na(Merged_Godo)] <- 0
  dim(Merged_Godo)
  Merged_Godo=Merged_Godo[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Godo, file="Godonkavirus_GOMWU.csv", row.names = FALSE)
  
  
  ##Halcy
  Halcy = read.csv('Corr_results_tau_Halcyonevirus.csv', check.names = FALSE)
  Halcy=setDT(Halcy, keep.rownames = TRUE)[]
  dim(Halcy)
  Halcy=Halcy[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Halcy) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Halcy, file="Halcyonevirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Halcy = merge(Halcy, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Halcy[is.na(Merged_Halcy)] <- 0
  dim(Merged_Halcy)
  Merged_Halcy=Merged_Halcy[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Halcy, file="Halcyonevirus_GOMWU.csv", row.names = FALSE)
  
  
  
  ##KO2
  KO2 = read.csv('Corr_results_tau_KlebsiellaphagephiKO2.csv', check.names = FALSE)
  KO2=setDT(KO2, keep.rownames = TRUE)[]
  dim(KO2)
  KO2=KO2[,c(1,242)]
  ##make the GOMWU matrix##
  colnames(KO2) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(KO2, file="KlebsiellaphagephiKO2_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_KO2 = merge(KO2, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_KO2[is.na(Merged_KO2)] <- 0
  dim(Merged_KO2)
  Merged_KO2=Merged_KO2[,c(1:2)]
  ##writeitout##
  write.csv(Merged_KO2, file="KlebsiellaphagephiKO2_GOMWU.csv", row.names = FALSE)
  
  
  ##Lacu
  Lacu = read.csv('Corr_results_tau_Lacusarxvirus.csv', check.names = FALSE)
  Lacu=setDT(Lacu, keep.rownames = TRUE)[]
  dim(Lacu)
  Lacu=Lacu[,c(1,186)]
  ##make the GOMWU matrix##
  colnames(Lacu) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Lacu, file="Lacusarxvirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Lacu = merge(Lacu, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Lacu[is.na(Merged_Lacu)] <- 0
  dim(Merged_Lacu)
  Merged_Lacu=Merged_Lacu[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Lacu, file="Lacusarxvirus_GOMWU.csv", row.names = FALSE)
  
  
  ##Light
  Light = read.csv('Corr_results_tau_Lightbulbvirus.csv', check.names = FALSE)
  Light=setDT(Light, keep.rownames = TRUE)[]
  dim(Light)
  Light=Light[,c(1,357)]
  ##make the GOMWU matrix##
  colnames(Light) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Light, file="Lightbulbvirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Light = merge(Light, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Light[is.na(Merged_Light)] <- 0
  dim(Merged_Light)
  Merged_Light=Merged_Light[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Light, file="Lightbulbvirus_GOMWU.csv", row.names = FALSE)
  
  
  ##Nona
  Nona = read.csv('Corr_results_tau_Nonagvirus.csv', check.names = FALSE)
  Nona=setDT(Nona, keep.rownames = TRUE)[]
  dim(Nona)
  Nona=Nona[,c(1,346)]
  ##make the GOMWU matrix##
  colnames(Nona) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Nona, file="Nonagvirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Nona = merge(Nona, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Nona[is.na(Merged_Nona)] <- 0
  dim(Merged_Nona)
  Merged_Nona=Merged_Nona[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Nona, file="Nonagvirus_GOMWU.csv", row.names = FALSE)
  
  
  ##H103
  H103 = read.csv('Corr_results_tau_PseudoalteromonasphageH103.csv', check.names = FALSE)
  H103=setDT(H103, keep.rownames = TRUE)[]
  dim(H103)
  H103=H103[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(H103) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(H103, file="PseudoalteromonasphageH103_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_H103 = merge(H103, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_H103[is.na(Merged_H103)] <- 0
  dim(Merged_H103)
  Merged_H103=Merged_H103[,c(1:2)]
  ##writeitout##
  write.csv(Merged_H103, file="PseudoalteromonasphageH103_GOMWU.csv", row.names = FALSE)
  
  
  ##JBD44
  JBD44 = read.csv('Corr_results_tau_PseudomonasphageJBD44.csv', check.names = FALSE)
  JBD44=setDT(JBD44, keep.rownames = TRUE)[]
  dim(JBD44)
  JBD44=JBD44[,c(1,381)]
  ##make the GOMWU matrix##
  colnames(JBD44) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(JBD44, file="PseudomonasphageJBD44_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_JBD44 = merge(JBD44, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_JBD44[is.na(Merged_JBD44)] <- 0
  dim(Merged_JBD44)
  Merged_JBD44=Merged_JBD44[,c(1:2)]
  ##writeitout##
  write.csv(Merged_JBD44, file="PseudomonasphageJBD44_GOMWU.csv", row.names = FALSE)
  
  
  ##Schit
  Schit = read.csv('Corr_results_tau_Schitoviridae.csv', check.names = FALSE)
  Schit=setDT(Schit, keep.rownames = TRUE)[]
  dim(Schit)
  Schit=Schit[,c(1,297)]
  ##make the GOMWU matrix##
  colnames(Schit) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Schit, file="Schitoviridae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Schit = merge(Schit, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Schit[is.na(Merged_Schit)] <- 0
  dim(Merged_Schit)
  Merged_Schit=Merged_Schit[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Schit, file="Schitoviridae_GOMWU.csv", row.names = FALSE)
  
  ##Spiz
  Spiz = read.csv('Corr_results_tau_Spizizenvirus.csv', check.names = FALSE)
  Spiz=setDT(Spiz, keep.rownames = TRUE)[]
  dim(Spiz)
  Spiz=Spiz[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Spiz) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Spiz, file="Spizizenvirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Spiz = merge(Spiz, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Spiz[is.na(Merged_Spiz)] <- 0
  dim(Merged_Spiz)
  Merged_Spiz=Merged_Spiz[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Spiz, file="Spizizenvirus_GOMWU.csv", row.names = FALSE)

  
  ##ARI0746
  ARI0746 = read.csv('Corr_results_tau_StreptococcusphagephiARI0746.csv', check.names = FALSE)
  ARI0746=setDT(ARI0746, keep.rownames = TRUE)[]
  dim(ARI0746)
  ARI0746=ARI0746[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(ARI0746) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(ARI0746, file="StreptococcusphagephiARI0746_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_ARI0746 = merge(ARI0746, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_ARI0746[is.na(Merged_ARI0746)] <- 0
  dim(Merged_ARI0746)
  Merged_ARI0746=Merged_ARI0746[,c(1:2)]
  ##writeitout##
  write.csv(Merged_ARI0746, file="StreptococcusphagephiARI0746_GOMWU.csv", row.names = FALSE)
  
  
  ##NJ2
  NJ2 = read.csv('Corr_results_tau_StreptococcusphagephiNJ2.csv', check.names = FALSE)
  NJ2=setDT(NJ2, keep.rownames = TRUE)[]
  dim(NJ2)
  NJ2=NJ2[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(NJ2) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(NJ2, file="StreptococcusphagephiNJ2_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_NJ2 = merge(NJ2, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_NJ2[is.na(Merged_NJ2)] <- 0
  dim(Merged_NJ2)
  Merged_NJ2=Merged_NJ2[,c(1:2)]
  ##writeitout##
  write.csv(Merged_NJ2, file="StreptococcusphagephiNJ2_GOMWU.csv", row.names = FALSE)

  
  ##Ueta
  Ueta = read.csv('Corr_results_tau_Uetakevirus.csv', check.names = FALSE)
  Ueta=setDT(Ueta, keep.rownames = TRUE)[]
  dim(Ueta)
  Ueta=Ueta[,c(1,320)]
  ##make the GOMWU matrix##
  colnames(Ueta) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Ueta, file="Uetakevirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Ueta = merge(Ueta, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Ueta[is.na(Merged_Ueta)] <- 0
  dim(Merged_Ueta)
  Merged_Ueta=Merged_Ueta[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Ueta, file="Uetakevirus_GOMWU.csv", row.names = FALSE)  
  
  
  ##Verto
  Verto = read.csv('Corr_results_tau_Vertoviridae.csv', check.names = FALSE)
  Verto=setDT(Verto, keep.rownames = TRUE)[]
  dim(Verto)
  Verto=Verto[,c(1,319)]
  ##make the GOMWU matrix##
  colnames(Verto) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Verto, file="Vertoviridae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Verto = merge(Verto, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Verto[is.na(Merged_Verto)] <- 0
  dim(Merged_Verto)
  Merged_Verto=Merged_Verto[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Verto, file="Vertoviridae_GOMWU.csv", row.names = FALSE)
  
  
  ##Vivid
  Vivid = read.csv('Corr_results_tau_Vividuovirus.csv', check.names = FALSE)
  Vivid=setDT(Vivid, keep.rownames = TRUE)[]
  dim(Vivid)
  Vivid=Vivid[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Vivid) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Vivid, file="Vividuovirus_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Vivid = merge(Vivid, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Vivid[is.na(Merged_Vivid)] <- 0
  dim(Merged_Vivid)
  Merged_Vivid=Merged_Vivid[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Vivid, file="Vividuovirus_GOMWU.csv", row.names = FALSE)

  
  ##Zobe
  Zobe = read.csv('Corr_results_tau_Zobellviridae.csv', check.names = FALSE)
  Zobe=setDT(Zobe, keep.rownames = TRUE)[]
  dim(Zobe)
  Zobe=Zobe[,c(1,385)]
  ##make the GOMWU matrix##
  colnames(Zobe) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Zobe, file="Zobellviridae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Zobe = merge(Zobe, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Zobe[is.na(Merged_Zobe)] <- 0
  dim(Merged_Zobe)
  Merged_Zobe=Merged_Zobe[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Zobe, file="Zobellviridae_GOMWU.csv", row.names = FALSE)  
  