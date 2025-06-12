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

  ##Blasto
  ##read in tau values##
  Blasto = read.csv('Corr_results_tau_Blastocatellaceae.csv', check.names = FALSE)
  Blasto=setDT(Blasto, keep.rownames = TRUE)[]
  dim(Blasto)
  ##keep just the column with ortho IDs and taus of associaitons with viral group##
  Blasto=Blasto[,c(1,261)]
  ##rename columns
  colnames(Blasto) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Blasto, file="Blastocatellaceae_taus.csv", row.names = FALSE)
  ###merge with list of all orthologs (necessary for GOMWU)
  Merged_Blasto = merge(Blasto, genenames, by = "gene", all.y = TRUE)
  ##replace NAs (i.e. orthos that weren't in your top 5%) with 0, following protocol from Fuess et al. mBio##
  Merged_Blasto[is.na(Merged_Blasto)] <- 0
  dim(Merged_Blasto)
  ##select just the ortho ID and tau value columns
  Merged_Blasto=Merged_Blasto[,c(1:2)]
  ##write it all out##
  write.csv(Merged_Blasto, file="Blastocatellaceae_GOMWU.csv", row.names = FALSE)
  
  
 
  ##Cami
  Cami = read.csv('Corr_results_tau_Caminicellaceae.csv', check.names = FALSE)
  Cami=setDT(Cami, keep.rownames = TRUE)[]
  dim(Cami)
  Cami=Cami[,c(1,262)]
  ##make the GOMWU matrix##
  colnames(Cami) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Cami, file="Caminicellaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Cami = merge(Cami, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Cami[is.na(Merged_Cami)] <- 0
  dim(Merged_Cami)
  Merged_Cami=Merged_Cami[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Cami, file="Caminicellaceae_GOMWU.csv", row.names = FALSE)
  
  ##Chloro
  Chloro = read.csv('Corr_results_tau_Chloroflexaceae.csv', check.names = FALSE)
  Chloro=setDT(Chloro, keep.rownames = TRUE)[]
  dim(Chloro)
  Chloro=Chloro[,c(1,227)]
  ##make the GOMWU matrix##
  colnames(Chloro) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Chloro, file="Chloroflexaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Chloro = merge(Chloro, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Chloro[is.na(Merged_Chloro)] <- 0
  dim(Merged_Chloro)
  Merged_Chloro=Merged_Chloro[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Chloro, file="Chloroflexaceae_GOMWU.csv", row.names = FALSE)
  
  ##Chris
  Chris = read.csv('Corr_results_tau_Christensenellaceae.csv', check.names = FALSE)
  Chris=setDT(Chris, keep.rownames = TRUE)[]
  dim(Chris)
  Chris=Chris[,c(1,295)]
  ##make the GOMWU matrix##
  colnames(Chris) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Chris, file="Christensenellaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Chris = merge(Chris, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Chris[is.na(Merged_Chris)] <- 0
  dim(Merged_Chris)
  Merged_Chris=Merged_Chris[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Chris, file="Christensenellaceae_GOMWU.csv", row.names = FALSE)
  
  ##Chrom
  Chrom = read.csv('Corr_results_tau_Chromatiaceae.csv', check.names = FALSE)
  Chrom=setDT(Chrom, keep.rownames = TRUE)[]
  dim(Chrom)
  Chrom=Chrom[,c(1,140)]
  ##make the GOMWU matrix##
  colnames(Chrom) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Chrom, file="Chromatiaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Chrom = merge(Chrom, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Chrom[is.na(Merged_Chrom)] <- 0
  dim(Merged_Chrom)
  Merged_Chrom=Merged_Chrom[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Chrom, file="Chromatiaceae_GOMWU.csv", row.names = FALSE)
  
  ##Desulf
  Desulf = read.csv('Corr_results_tau_Desulfobacteraceae.csv', check.names = FALSE)
  Desulf=setDT(Desulf, keep.rownames = TRUE)[]
  dim(Desulf)
  Desulf=Desulf[,c(1,181)]
  ##make the GOMWU matrix##
  colnames(Desulf) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Desulf, file="Desulfobacteraceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Desulf = merge(Desulf, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Desulf[is.na(Merged_Desulf)] <- 0
  dim(Merged_Desulf)
  Merged_Desulf=Merged_Desulf[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Desulf, file="Desulfobacteraceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Desulf
  Desulf = read.csv('Corr_results_tau_Desulfotomaculales.csv', check.names = FALSE)
  Desulf=setDT(Desulf, keep.rownames = TRUE)[]
  dim(Desulf)
  Desulf=Desulf[,c(1,282)]
  ##make the GOMWU matrix##
  colnames(Desulf) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Desulf, file="Desulfotomaculales_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Desulf = merge(Desulf, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Desulf[is.na(Merged_Desulf)] <- 0
  dim(Merged_Desulf)
  Merged_Desulf=Merged_Desulf[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Desulf, file="Desulfotomaculales_GOMWU.csv", row.names = FALSE)
  
  
  ##Deth
  Deth = read.csv('Corr_results_tau_Dethiobacteraceae.csv', check.names = FALSE)
  Deth=setDT(Deth, keep.rownames = TRUE)[]
  dim(Deth)
  Deth=Deth[,c(1,308)]
  ##make the GOMWU matrix##
  colnames(Deth) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Deth, file="Dethiobacteraceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Deth = merge(Deth, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Deth[is.na(Merged_Deth)] <- 0
  dim(Merged_Deth)
  Merged_Deth=Merged_Deth[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Deth, file="Dethiobacteraceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Endo
  Endo = read.csv('Corr_results_tau_Endozoicomonadaceae.csv', check.names = FALSE)
  Endo=setDT(Endo, keep.rownames = TRUE)[]
  dim(Endo)
  Endo=Endo[,c(1,230)]
  ##make the GOMWU matrix##
  colnames(Endo) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Endo, file="Endozoicomonadaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Endo = merge(Endo, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Endo[is.na(Merged_Endo)] <- 0
  dim(Merged_Endo)
  Merged_Endo=Merged_Endo[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Endo, file="Endozoicomonadaceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Halo
  Halo = read.csv('Corr_results_tau_Halobacteroidaceae.csv', check.names = FALSE)
  Halo=setDT(Halo, keep.rownames = TRUE)[]
  dim(Halo)
  Halo=Halo[,c(1,241)]
  ##make the GOMWU matrix##
  colnames(Halo) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Halo, file="Halobacteroidaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Halo = merge(Halo, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Halo[is.na(Merged_Halo)] <- 0
  dim(Merged_Halo)
  Merged_Halo=Merged_Halo[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Halo, file="Halobacteroidaceae_GOMWU.csv", row.names = FALSE)
  
  
  
  ##Inqui
  Inqui = read.csv('Corr_results_tau_Inquilinaceae.csv', check.names = FALSE)
  Inqui=setDT(Inqui, keep.rownames = TRUE)[]
  dim(Inqui)
  Inqui=Inqui[,c(1,163)]
  ##make the GOMWU matrix##
  colnames(Inqui) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Inqui, file="Inquilinaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Inqui = merge(Inqui, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Inqui[is.na(Merged_Inqui)] <- 0
  dim(Merged_Inqui)
  Merged_Inqui=Merged_Inqui[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Inqui, file="Inquilinaceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Kilon
  Kilon = read.csv('Corr_results_tau_Kiloniellaceae.csv', check.names = FALSE)
  Kilon=setDT(Kilon, keep.rownames = TRUE)[]
  dim(Kilon)
  Kilon=Kilon[,c(1,203)]
  ##make the GOMWU matrix##
  colnames(Kilon) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Kilon, file="Kiloniellaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Kilon = merge(Kilon, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Kilon[is.na(Merged_Kilon)] <- 0
  dim(Merged_Kilon)
  Merged_Kilon=Merged_Kilon[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Kilon, file="Kiloniellaceae_GOMWU.csv", row.names = FALSE)
  
  ##Len
  Len = read.csv('Corr_results_tau_Lentimicrobiaceae.csv', check.names = FALSE)
  Len=setDT(Len, keep.rownames = TRUE)[]
  dim(Len)
  Len=Len[,c(1,159)]
  ##make the GOMWU matrix##
  colnames(Len) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Len, file="Lenmicrobiaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Len = merge(Len, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Len[is.na(Merged_Len)] <- 0
  dim(Merged_Len)
  Merged_Len=Merged_Len[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Len, file="Lentimicrobiaceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Lenti
  Lenti = read.csv('Corr_results_tau_Lentisphaeraceae.csv', check.names = FALSE)
  Lenti=setDT(Lenti, keep.rownames = TRUE)[]
  dim(Lenti)
  Lenti=Lenti[,c(1,242)]
  ##make the GOMWU matrix##
  colnames(Lenti) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Lenti, file="Lentisphaeraceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Lenti = merge(Lenti, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Lenti[is.na(Merged_Lenti)] <- 0
  dim(Merged_Lenti)
  Merged_Lenti=Merged_Lenti[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Lenti, file="Lentisphaeraceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Nisa
  Nisa = read.csv('Corr_results_tau_Nisaeaceae.csv', check.names = FALSE)
  Nisa=setDT(Nisa, keep.rownames = TRUE)[]
  dim(Nisa)
  Nisa=Nisa[,c(1,259)]
  ##make the GOMWU matrix##
  colnames(Nisa) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Nisa, file="Nisaeaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Nisa = merge(Nisa, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Nisa[is.na(Merged_Nisa)] <- 0
  dim(Merged_Nisa)
  Merged_Nisa=Merged_Nisa[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Nisa, file="Nisaeaceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Oligo
  Oligo = read.csv('Corr_results_tau_Oligoflexaceae.csv', check.names = FALSE)
  Oligo=setDT(Oligo, keep.rownames = TRUE)[]
  dim(Oligo)
  Oligo=Oligo[,c(1,278)]
  ##make the GOMWU matrix##
  colnames(Oligo) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Oligo, file="Oligoflexaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Oligo = merge(Oligo, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Oligo[is.na(Merged_Oligo)] <- 0
  dim(Merged_Oligo)
  Merged_Oligo=Merged_Oligo[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Oligo, file="Oligoflexaceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Prolix
  Prolix = read.csv('Corr_results_tau_Prolixibacteraceae.csv', check.names = FALSE)
  Prolix=setDT(Prolix, keep.rownames = TRUE)[]
  dim(Prolix)
  Prolix=Prolix[,c(1,231)]
  ##make the GOMWU matrix##
  colnames(Prolix) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Prolix, file="Prolixibacteraceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Prolix = merge(Prolix, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Prolix[is.na(Merged_Prolix)] <- 0
  dim(Merged_Prolix)
  Merged_Prolix=Merged_Prolix[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Prolix, file="Prolixibacteraceae_GOMWU.csv", row.names = FALSE)
  
  
  ##Rhodo
  Rhodo = read.csv('Corr_results_tau_Rhodocyclaceae.csv', check.names = FALSE)
  Rhodo=setDT(Rhodo, keep.rownames = TRUE)[]
  dim(Rhodo)
  Rhodo=Rhodo[,c(1,70)]
  ##make the GOMWU matrix##
  colnames(Rhodo) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Rhodo, file="Rhodocyclaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Rhodo = merge(Rhodo, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Rhodo[is.na(Merged_Rhodo)] <- 0
  dim(Merged_Rhodo)
  Merged_Rhodo=Merged_Rhodo[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Rhodo, file="Rhodocyclaceae_GOMWU.csv", row.names = FALSE)
  
  ##Spir
  Spir = read.csv('Corr_results_tau_Spirochaetaceae.csv', check.names = FALSE)
  Spir=setDT(Spir, keep.rownames = TRUE)[]
  dim(Spir)
  Spir=Spir[,c(1,232)]
  ##make the GOMWU matrix##
  colnames(Spir) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Spir, file="Spirochaetaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Spir = merge(Spir, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Spir[is.na(Merged_Spir)] <- 0
  dim(Merged_Spir)
  Merged_Spir=Merged_Spir[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Spir, file="Spirochaetaceae_GOMWU.csv", row.names = FALSE)

  
  ##Tera
  Tera = read.csv('Corr_results_tau_Terasakiellaceae.csv', check.names = FALSE)
  Tera=setDT(Tera, keep.rownames = TRUE)[]
  dim(Tera)
  Tera=Tera[,c(1,240)]
  ##make the GOMWU matrix##
  colnames(Tera) <- c("gene", "tau")
  ##write out the tau values for supplemental file##
  write.csv(Tera, file="Terasakiellaceae_taus.csv", row.names = FALSE)
  ##finish up GOMWU input##
  Merged_Tera = merge(Tera, genenames, by = "gene", all.y = TRUE)
  ##replace NA with 0##
  Merged_Tera[is.na(Merged_Tera)] <- 0
  dim(Merged_Tera)
  Merged_Tera=Merged_Tera[,c(1:2)]
  ##writeitout##
  write.csv(Merged_Tera, file="Terasakiellaceae_GOMWU.csv", row.names = FALSE)
  
  