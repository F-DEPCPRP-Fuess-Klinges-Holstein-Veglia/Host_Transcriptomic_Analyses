##this code parses the output of our correlation code which generates pvalues for all pariwise comparisons
##specifically we determine our top 5% most correlated microbe families and genes and create new matrices for each so we can calucate tau values for downstream analyses
##last updated LEF 6/12/25

##load packages
library(janitor)
library(tibble)
library(dplyr)
library(data.table)

##read in p value matrix
  p=read.csv("Corr_results_pval.csv")
  colnames(p)
  ##reduce down to just grab those columns corresponding to microbe families
  p_reduced=p[,c(1:390)] 
  colnames(p_reduced)
  rownames(p_reduced)
  ##reduce down to just grap those rows corresponding to orthologs
  p_reduced=p_reduced[c(391:7852),]
  ##set non sig values to NA for processing
  p_reduced[p_reduced>0.05]=NA

##and now we look for those families in the top 5%##
  P_fam=as.data.frame(t(p_reduced))
  ##start by replacing non-sig (n/as) with 1000000)
  P_fam[is.na(P_fam)] <- 1000000
  ##create a column for row sums##
  P_fam$sums = rowSums(P_fam)
  ##order by sums##
  P_fam_sorted <- P_fam[order(P_fam$sums),] 
  P_fam_sorted$sums
  ##pull out the top 5% of ASV (roughly 20; pull 21 cause 1 is unknown)##
  P_fam_top_sorted = P_fam_sorted[1:21,]
  row.names(P_fam_top_sorted)
  is.na(P_fam_top_sorted) <- P_fam_top_sorted == 1000000
  ##generate a list of top families for reference##
  P_fam_familes = P_fam_top_sorted[1]
  write.csv(P_fam_familes, file="top5families.csv")

##ok now we want just the top 5% of genes that are correlated to any given family##
  All_genes=p_reduced
  ##and follow code from above##
  ##start by replacing non-sig (n/as) with 1000000)
  All_genes[is.na(All_genes)] <- 1000000
  lapply(All_genes,class)
  ##create a column for row sums##
  All_genes$sums = rowSums(All_genes)
  ##order by sums##
  All_sorted_genes <- All_genes[order(All_genes$sums),] 
  ##pull out the top 5% of genes (roughly 373)##
  All_sorted_genes$sums
  ##we actually do 376 because there is a large tie for genes which are correlated to the same number of genes; this is the best break##
  ##code which validates this point is shown at the end of this script##
  All_top_sorted_genes = All_sorted_genes[1:376,]
  top_genes = All_top_sorted_genes[1]
  write.csv(top_genes, file="top5genes.csv")
  is.na(All_top_sorted_genes) <- All_top_sorted_genes == 1000000

##now parse out lists of significant genes in top 5% of both family and genes##
##and create matrices for finding tau values
  ##processing files##
  reads=read.csv("normalized_reads_corr.csv")
  reads=t(reads)
  reads=reads %>%
    row_to_names(row_number = 1)
  reads=as.data.frame(reads)
  reads <- tibble::rownames_to_column(reads, "rn")
  
  ##read in Prop data and format##
  Prop = read.csv("families_prop.csv", check.names = FALSE)
  Prop=t(Prop)
  Prop=Prop %>%
    row_to_names(row_number = 1)
  Prop=as.data.frame(Prop)
  Prop <- tibble::rownames_to_column(Prop, "sample_R")
  
  ##read in sample names##
  meta=read.csv("metadata_matchasv.csv")
  meta=meta[,c(1:2)]
  
  ##now use the dataframes we've generated to create lists of sig genes
  ##also create matrices for getting tau values of associatons between each top family and associated top genes.
  
  ##Kiloniellaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Kiloniellaceae = All_top_sorted_genes[,c("Kiloniellaceae_Prop","sums")]
    Kiloniellaceae = Kiloniellaceae[complete.cases(Kiloniellaceae), ]
    Kiloniellaceae = setDT(Kiloniellaceae, keep.rownames = TRUE)[]
    Kiloniellaceae = Kiloniellaceae[,c(1)]
    write.csv(Kiloniellaceae, "Kiloniellaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    KilonMatrix = merge(Kiloniellaceae, reads,
                      by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    KilonMatrix = KilonMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    KilonMatrix = as.data.frame(t(KilonMatrix))
    KilonMatrix <- tibble::rownames_to_column(KilonMatrix, "sample_R")
    Kilon = Prop[c("sample_R","Kiloniellaceae_Prop")]
    Kilon <- merge(Kilon,meta, by.x="sample_R", by.y="sample.ID")
    Kilon <- Kilon[,c(3,2)]
    colnames(Kilon)[1] <- "sample_R"
    KilonMatrix_final = merge(KilonMatrix, Kilon,
                       by = "sample_R", all.x = TRUE,
                       sort = TRUE, no.dups = FALSE)
    dim(KilonMatrix_final)
    KilonMatrix_final = KilonMatrix_final[,c(2:203)]
    write.csv(KilonMatrix_final, file = "Kiloniellaceae_Tau_Matrix.csv", row.names = FALSE)
    
    

  ##Endozoicomonadaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Endozoicomonadaceae = All_top_sorted_genes[,c("Endozoicomonadaceae_Prop","sums")]
    Endozoicomonadaceae = Endozoicomonadaceae[complete.cases(Endozoicomonadaceae), ]
    Endozoicomonadaceae = setDT(Endozoicomonadaceae, keep.rownames = TRUE)[]
    Endozoicomonadaceae = Endozoicomonadaceae[,c(1)]
    write.csv(Endozoicomonadaceae, "Endozoicomonadaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    EndoMatrix = merge(Endozoicomonadaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    EndoMatrix = EndoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    EndoMatrix = as.data.frame(t(EndoMatrix))
    EndoMatrix <- tibble::rownames_to_column(EndoMatrix, "sample_R")
    Endo = Prop[c("sample_R","Endozoicomonadaceae_Prop")]
    Endo <- merge(Endo,meta, by.x="sample_R", by.y="sample.ID")
    Endo <- Endo[,c(3,2)]
    colnames(Endo)[1] <- "sample_R"
    EndoMatrix_final = merge(EndoMatrix, Endo,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(EndoMatrix_final)
    EndoMatrix_final = EndoMatrix_final[,c(2:230)]
    write.csv(EndoMatrix_final, file = "Endozoicomonadaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Nisaeaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Nisaeaceae = All_top_sorted_genes[,c("Nisaeaceae_Prop","sums")]
    Nisaeaceae = Nisaeaceae[complete.cases(Nisaeaceae), ]
    Nisaeaceae = setDT(Nisaeaceae, keep.rownames = TRUE)[]
    Nisaeaceae = Nisaeaceae[,c(1)]
    write.csv(Nisaeaceae, "Nisaeaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    NisaeMatrix = merge(Nisaeaceae, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    NisaeMatrix = NisaeMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    NisaeMatrix = as.data.frame(t(NisaeMatrix))
    NisaeMatrix <- tibble::rownames_to_column(NisaeMatrix, "sample_R")
    Nisae = Prop[c("sample_R","Nisaeaceae_Prop")]
    Nisae <- merge(Nisae,meta, by.x="sample_R", by.y="sample.ID")
    Nisae <- Nisae[,c(3,2)]
    colnames(Nisae)[1] <- "sample_R"
    NisaeMatrix_final = merge(NisaeMatrix, Nisae,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(NisaeMatrix_final)
    NisaeMatrix_final = NisaeMatrix_final[,c(2:259)]
    write.csv(NisaeMatrix_final, file = "Nisaeaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Prolixibacteraceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Prolixibacteraceae = All_top_sorted_genes[,c("Prolixibacteraceae_Prop","sums")]
    Prolixibacteraceae = Prolixibacteraceae[complete.cases(Prolixibacteraceae), ]
    Prolixibacteraceae = setDT(Prolixibacteraceae, keep.rownames = TRUE)[]
    Prolixibacteraceae = Prolixibacteraceae[,c(1)]
    write.csv(Prolixibacteraceae, "Prolixibacteraceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ProlixMatrix = merge(Prolixibacteraceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    ProlixMatrix = ProlixMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ProlixMatrix = as.data.frame(t(ProlixMatrix))
    ProlixMatrix <- tibble::rownames_to_column(ProlixMatrix, "sample_R")
    Prolix = Prop[c("sample_R","Prolixibacteraceae_Prop")]
    Prolix <- merge(Prolix,meta, by.x="sample_R", by.y="sample.ID")
    Prolix <- Prolix[,c(3,2)]
    colnames(Prolix)[1] <- "sample_R"
    ProlixMatrix_final = merge(ProlixMatrix, Prolix,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(ProlixMatrix_final)
    ProlixMatrix_final = ProlixMatrix_final[,c(2:231)]
    write.csv(ProlixMatrix_final, file = "Prolixibacteraceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Desulfobacteraceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Desulfobacteraceae = All_top_sorted_genes[,c("Desulfobacteraceae_Prop","sums")]
    Desulfobacteraceae = Desulfobacteraceae[complete.cases(Desulfobacteraceae), ]
    Desulfobacteraceae = setDT(Desulfobacteraceae, keep.rownames = TRUE)[]
    Desulfobacteraceae = Desulfobacteraceae[,c(1)]
    write.csv(Desulfobacteraceae, "Desulfobacteraceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    DesulfoMatrix = merge(Desulfobacteraceae, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    DesulfoMatrix = DesulfoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    DesulfoMatrix = as.data.frame(t(DesulfoMatrix))
    DesulfoMatrix <- tibble::rownames_to_column(DesulfoMatrix, "sample_R")
    Desulfo = Prop[c("sample_R","Desulfobacteraceae_Prop")]
    Desulfo <- merge(Desulfo,meta, by.x="sample_R", by.y="sample.ID")
    Desulfo <- Desulfo[,c(3,2)]
    colnames(Desulfo)[1] <- "sample_R"
    DesulfoMatrix_final = merge(DesulfoMatrix, Desulfo,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(DesulfoMatrix_final)
    DesulfoMatrix_final = DesulfoMatrix_final[,c(2:181)]
    write.csv(DesulfoMatrix_final, file = "Desulfobacteraceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Rhodocyclaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Rhodocyclaceae = All_top_sorted_genes[,c("Rhodocyclaceae_Prop","sums")]
    Rhodocyclaceae = Rhodocyclaceae[complete.cases(Rhodocyclaceae), ]
    Rhodocyclaceae = setDT(Rhodocyclaceae, keep.rownames = TRUE)[]
    Rhodocyclaceae = Rhodocyclaceae[,c(1)]
    write.csv(Rhodocyclaceae, "Rhodocyclaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    RhodoMatrix = merge(Rhodocyclaceae, reads,
                          by = "rn", all.x = TRUE,
                          sort = TRUE,  no.dups = FALSE,)
    RhodoMatrix = RhodoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    RhodoMatrix = as.data.frame(t(RhodoMatrix))
    RhodoMatrix <- tibble::rownames_to_column(RhodoMatrix, "sample_R")
    Rhodo = Prop[c("sample_R","Rhodocyclaceae_Prop")]
    Rhodo <- merge(Rhodo,meta, by.x="sample_R", by.y="sample.ID")
    Rhodo <- Rhodo[,c(3,2)]
    colnames(Rhodo)[1] <- "sample_R"
    RhodoMatrix_final = merge(RhodoMatrix, Rhodo,
                                by = "sample_R", all.x = TRUE,
                                sort = TRUE, no.dups = FALSE)
    dim(RhodoMatrix_final)
    RhodoMatrix_final = RhodoMatrix_final[,c(2:70)]
    write.csv(RhodoMatrix_final, file = "Rhodocyclaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  
  ##Blastocatellaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Blastocatellaceae = All_top_sorted_genes[,c("Blastocatellaceae_Prop","sums")]
    Blastocatellaceae = Blastocatellaceae[complete.cases(Blastocatellaceae), ]
    Blastocatellaceae = setDT(Blastocatellaceae, keep.rownames = TRUE)[]
    Blastocatellaceae = Blastocatellaceae[,c(1)]
    write.csv(Blastocatellaceae, "Blastocatellaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    BlastoMatrix = merge(Blastocatellaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    BlastoMatrix = BlastoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    BlastoMatrix = as.data.frame(t(BlastoMatrix))
    BlastoMatrix <- tibble::rownames_to_column(BlastoMatrix, "sample_R")
    Blasto = Prop[c("sample_R","Blastocatellaceae_Prop")]
    Blasto <- merge(Blasto,meta, by.x="sample_R", by.y="sample.ID")
    Blasto <- Blasto[,c(3,2)]
    colnames(Blasto)[1] <- "sample_R"
    BlastoMatrix_final = merge(BlastoMatrix, Blasto,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(BlastoMatrix_final)
    BlastoMatrix_final = BlastoMatrix_final[,c(2:261)]
    write.csv(BlastoMatrix_final, file = "Blastocatellaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Halobacteroidaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Halobacteroidaceae = All_top_sorted_genes[,c("Halobacteroidaceae_Prop","sums")]
    Halobacteroidaceae = Halobacteroidaceae[complete.cases(Halobacteroidaceae), ]
    Halobacteroidaceae = setDT(Halobacteroidaceae, keep.rownames = TRUE)[]
    Halobacteroidaceae = Halobacteroidaceae[,c(1)]
    write.csv(Halobacteroidaceae, "Halobacteroidaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    HaloMatrix = merge(Halobacteroidaceae, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    HaloMatrix = HaloMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    HaloMatrix = as.data.frame(t(HaloMatrix))
    HaloMatrix <- tibble::rownames_to_column(HaloMatrix, "sample_R")
    Halo = Prop[c("sample_R","Halobacteroidaceae_Prop")]
    Halo <- merge(Halo,meta, by.x="sample_R", by.y="sample.ID")
    Halo <- Halo[,c(3,2)]
    colnames(Halo)[1] <- "sample_R"
    HaloMatrix_final = merge(HaloMatrix, Halo,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(HaloMatrix_final)
    HaloMatrix_final = HaloMatrix_final[,c(2:241)]
    write.csv(HaloMatrix_final, file = "Halobacteroidaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Desulfotomaculales
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Desulfotomaculales = All_top_sorted_genes[,c("Desulfotomaculales.Incertae.Sedis_Prop","sums")]
    Desulfotomaculales = Desulfotomaculales[complete.cases(Desulfotomaculales), ]
    Desulfotomaculales = setDT(Desulfotomaculales, keep.rownames = TRUE)[]
    Desulfotomaculales = Desulfotomaculales[,c(1)]
    write.csv(Desulfotomaculales, "Desulfotomaculales_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    DesulfMatrix = merge(Desulfotomaculales, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    DesulfMatrix = DesulfMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    DesulfMatrix = as.data.frame(t(DesulfMatrix))
    DesulfMatrix <- tibble::rownames_to_column(DesulfMatrix, "sample_R")
    Desulf = Prop[c("sample_R","Desulfotomaculales Incertae Sedis_Prop")]
    Desulf <- merge(Desulf,meta, by.x="sample_R", by.y="sample.ID")
    Desulf <- Desulf[,c(3,2)]
    colnames(Desulf)[1] <- "sample_R"
    DesulfMatrix_final = merge(DesulfMatrix, Desulf,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(DesulfMatrix_final)
    DesulfMatrix_final = DesulfMatrix_final[,c(2:282)]
    write.csv(DesulfMatrix_final, file = "Desulfotomaculales_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Lentimicrobiaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Lentimicrobiaceae = All_top_sorted_genes[,c("Lentimicrobiaceae_Prop","sums")]
    Lentimicrobiaceae = Lentimicrobiaceae[complete.cases(Lentimicrobiaceae), ]
    Lentimicrobiaceae = setDT(Lentimicrobiaceae, keep.rownames = TRUE)[]
    Lentimicrobiaceae = Lentimicrobiaceae[,c(1)]
    write.csv(Lentimicrobiaceae, "Lentimicrobiaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    LentiMatrix = merge(Lentimicrobiaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    LentiMatrix = LentiMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    LentiMatrix = as.data.frame(t(LentiMatrix))
    LentiMatrix <- tibble::rownames_to_column(LentiMatrix, "sample_R")
    Lenti = Prop[c("sample_R","Lentimicrobiaceae_Prop")]
    Lenti <- merge(Lenti,meta, by.x="sample_R", by.y="sample.ID")
    Lenti <- Lenti[,c(3,2)]
    colnames(Lenti)[1] <- "sample_R"
    LentiMatrix_final = merge(LentiMatrix, Lenti,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(LentiMatrix_final)
    LentiMatrix_final = LentiMatrix_final[,c(2:159)]
    write.csv(LentiMatrix_final, file = "Lentimicrobiaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Dethiobacteraceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Dethiobacteraceae = All_top_sorted_genes[,c("Dethiobacteraceae_Prop","sums")]
    Dethiobacteraceae = Dethiobacteraceae[complete.cases(Dethiobacteraceae), ]
    Dethiobacteraceae = setDT(Dethiobacteraceae, keep.rownames = TRUE)[]
    Dethiobacteraceae = Dethiobacteraceae[,c(1)]
    write.csv(Dethiobacteraceae, "Dethiobacteraceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    DethMatrix = merge(Dethiobacteraceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    DethMatrix = DethMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    DethMatrix = as.data.frame(t(DethMatrix))
    DethMatrix <- tibble::rownames_to_column(DethMatrix, "sample_R")
    Deth = Prop[c("sample_R","Dethiobacteraceae_Prop")]
    Deth <- merge(Deth,meta, by.x="sample_R", by.y="sample.ID")
    Deth <- Deth[,c(3,2)]
    colnames(Deth)[1] <- "sample_R"
    DethMatrix_final = merge(DethMatrix, Deth,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(DethMatrix_final)
    DethMatrix_final = DethMatrix_final[,c(2:308)]
    write.csv(DethMatrix_final, file = "Dethiobacteraceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  
  ##Terasakiellaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Terasakiellaceae = All_top_sorted_genes[,c("Terasakiellaceae_Prop","sums")]
    Terasakiellaceae = Terasakiellaceae[complete.cases(Terasakiellaceae), ]
    Terasakiellaceae = setDT(Terasakiellaceae, keep.rownames = TRUE)[]
    Terasakiellaceae = Terasakiellaceae[,c(1)]
    write.csv(Terasakiellaceae, "Terasakiellaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    TerasMatrix = merge(Terasakiellaceae, reads,
                          by = "rn", all.x = TRUE,
                          sort = TRUE,  no.dups = FALSE,)
    TerasMatrix = TerasMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    TerasMatrix = as.data.frame(t(TerasMatrix))
    TerasMatrix <- tibble::rownames_to_column(TerasMatrix, "sample_R")
    Teras = Prop[c("sample_R","Terasakiellaceae_Prop")]
    Teras <- merge(Teras,meta, by.x="sample_R", by.y="sample.ID")
    Teras <- Teras[,c(3,2)]
    colnames(Teras)[1] <- "sample_R"
    TerasMatrix_final = merge(TerasMatrix, Teras,
                                by = "sample_R", all.x = TRUE,
                                sort = TRUE, no.dups = FALSE)
    dim(TerasMatrix_final)
    TerasMatrix_final = TerasMatrix_final[,c(2:240)]
    write.csv(TerasMatrix_final, file = "Terasakiellaceae_Tau_Matrix.csv", row.names = FALSE)
  
  ##Caminicellaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Caminicellaceae = All_top_sorted_genes[,c("Caminicellaceae_Prop","sums")]
    Caminicellaceae = Caminicellaceae[complete.cases(Caminicellaceae), ]
    Caminicellaceae = setDT(Caminicellaceae, keep.rownames = TRUE)[]
    Caminicellaceae = Caminicellaceae[,c(1)]
    write.csv(Caminicellaceae, "Caminicellaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    CaminMatrix = merge(Caminicellaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    CaminMatrix = CaminMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    CaminMatrix = as.data.frame(t(CaminMatrix))
    CaminMatrix <- tibble::rownames_to_column(CaminMatrix, "sample_R")
    Camin = Prop[c("sample_R","Caminicellaceae_Prop")]
    Camin <- merge(Camin,meta, by.x="sample_R", by.y="sample.ID")
    Camin <- Camin[,c(3,2)]
    colnames(Camin)[1] <- "sample_R"
    CaminMatrix_final = merge(CaminMatrix, Camin,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(CaminMatrix_final)
    CaminMatrix_final = CaminMatrix_final[,c(2:262)]
    write.csv(CaminMatrix_final, file = "Caminicellaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  
  
  ##Chromatiaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Chromatiaceae = All_top_sorted_genes[,c("Chromatiaceae_Prop","sums")]
    Chromatiaceae = Chromatiaceae[complete.cases(Chromatiaceae), ]
    Chromatiaceae = setDT(Chromatiaceae, keep.rownames = TRUE)[]
    Chromatiaceae = Chromatiaceae[,c(1)]
    write.csv(Chromatiaceae, "Chromatiaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ChromMatrix = merge(Chromatiaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    ChromMatrix = ChromMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ChromMatrix = as.data.frame(t(ChromMatrix))
    ChromMatrix <- tibble::rownames_to_column(ChromMatrix, "sample_R")
    Chrom = Prop[c("sample_R","Chromatiaceae_Prop")]
    Chrom <- merge(Chrom,meta, by.x="sample_R", by.y="sample.ID")
    Chrom <- Chrom[,c(3,2)]
    colnames(Chrom)[1] <- "sample_R"
    ChromMatrix_final = merge(ChromMatrix, Chrom,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(ChromMatrix_final)
    ChromMatrix_final = ChromMatrix_final[,c(2:140)]
    write.csv(ChromMatrix_final, file = "Chromatiaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Christensenellaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Christensenellaceae = All_top_sorted_genes[,c("Christensenellaceae_Prop","sums")]
    Christensenellaceae = Christensenellaceae[complete.cases(Christensenellaceae), ]
    Christensenellaceae = setDT(Christensenellaceae, keep.rownames = TRUE)[]
    Christensenellaceae = Christensenellaceae[,c(1)]
    write.csv(Christensenellaceae, "Christensenellaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ChristMatrix = merge(Christensenellaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    ChristMatrix = ChristMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ChristMatrix = as.data.frame(t(ChristMatrix))
    ChristMatrix <- tibble::rownames_to_column(ChristMatrix, "sample_R")
    Christ = Prop[c("sample_R","Christensenellaceae_Prop")]
    Christ <- merge(Christ,meta, by.x="sample_R", by.y="sample.ID")
    Christ <- Christ[,c(3,2)]
    colnames(Christ)[1] <- "sample_R"
    ChristMatrix_final = merge(ChristMatrix, Christ,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(ChristMatrix_final)
    ChristMatrix_final = ChristMatrix_final[,c(2:295)]
    write.csv(ChristMatrix_final, file = "Christensenellaceae_Tau_Matrix.csv", row.names = FALSE)

  
  ##Inquilinaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Inquilinaceae = All_top_sorted_genes[,c("Inquilinaceae_Prop","sums")]
    Inquilinaceae = Inquilinaceae[complete.cases(Inquilinaceae), ]
    Inquilinaceae = setDT(Inquilinaceae, keep.rownames = TRUE)[]
    Inquilinaceae = Inquilinaceae[,c(1)]
    write.csv(Inquilinaceae, "Inquilinaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    InquMatrix = merge(Inquilinaceae, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    InquMatrix = InquMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    InquMatrix = as.data.frame(t(InquMatrix))
    InquMatrix <- tibble::rownames_to_column(InquMatrix, "sample_R")
    Inqu = Prop[c("sample_R","Inquilinaceae_Prop")]
    Inqu <- merge(Inqu,meta, by.x="sample_R", by.y="sample.ID")
    Inqu <- Inqu[,c(3,2)]
    colnames(Inqu)[1] <- "sample_R"
    InquMatrix_final = merge(InquMatrix, Inqu,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(InquMatrix_final)
    InquMatrix_final = InquMatrix_final[,c(2:163)]
    write.csv(InquMatrix_final, file = "Inquilinaceae_Tau_Matrix.csv", row.names = FALSE)

  
  ##Spirochaetaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Spirochaetaceae = All_top_sorted_genes[,c("Spirochaetaceae_Prop","sums")]
    Spirochaetaceae = Spirochaetaceae[complete.cases(Spirochaetaceae), ]
    Spirochaetaceae = setDT(Spirochaetaceae, keep.rownames = TRUE)[]
    Spirochaetaceae = Spirochaetaceae[,c(1)]
    write.csv(Spirochaetaceae, "Spirochaetaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    SpirouMatrix = merge(Spirochaetaceae, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    SpirouMatrix = SpirouMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    SpirouMatrix = as.data.frame(t(SpirouMatrix))
    SpirouMatrix <- tibble::rownames_to_column(SpirouMatrix, "sample_R")
    Spirou = Prop[c("sample_R","Spirochaetaceae_Prop")]
    Spirou <- merge(Spirou,meta, by.x="sample_R", by.y="sample.ID")
    Spirou <- Spirou[,c(3,2)]
    colnames(Spirou)[1] <- "sample_R"
    SpirouMatrix_final = merge(SpirouMatrix, Spirou,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(SpirouMatrix_final)
    SpirouMatrix_final = SpirouMatrix_final[,c(2:232)]
    write.csv(SpirouMatrix_final, file = "Spirochaetaceae_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Chloroflexaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Chloroflexaceae = All_top_sorted_genes[,c("Chloroflexaceae_Prop","sums")]
    Chloroflexaceae = Chloroflexaceae[complete.cases(Chloroflexaceae), ]
    Chloroflexaceae = setDT(Chloroflexaceae, keep.rownames = TRUE)[]
    Chloroflexaceae = Chloroflexaceae[,c(1)]
    write.csv(Chloroflexaceae, "Chloroflexaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ChloroMatrix = merge(Chloroflexaceae, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    ChloroMatrix = ChloroMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ChloroMatrix = as.data.frame(t(ChloroMatrix))
    ChloroMatrix <- tibble::rownames_to_column(ChloroMatrix, "sample_R")
    Chloro = Prop[c("sample_R","Chloroflexaceae_Prop")]
    Chloro <- merge(Chloro,meta, by.x="sample_R", by.y="sample.ID")
    Chloro <- Chloro[,c(3,2)]
    colnames(Chloro)[1] <- "sample_R"
    ChloroMatrix_final = merge(ChloroMatrix, Chloro,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(ChloroMatrix_final)
    ChloroMatrix_final = ChloroMatrix_final[,c(2:227)]
    write.csv(ChloroMatrix_final, file = "Chloroflexaceae_Tau_Matrix.csv", row.names = FALSE)

  
  ##Lentisphaeraceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Lentisphaeraceae = All_top_sorted_genes[,c("Lentisphaeraceae_Prop","sums")]
    Lentisphaeraceae = Lentisphaeraceae[complete.cases(Lentisphaeraceae), ]
    Lentisphaeraceae = setDT(Lentisphaeraceae, keep.rownames = TRUE)[]
    Lentisphaeraceae = Lentisphaeraceae[,c(1)]
    write.csv(Lentisphaeraceae, "Lentisphaeraceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    LenMatrix = merge(Lentisphaeraceae, reads,
                      by = "rn", all.x = TRUE,
                      sort = TRUE,  no.dups = FALSE,)
    LenMatrix = LenMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    LenMatrix = as.data.frame(t(LenMatrix))
    LenMatrix <- tibble::rownames_to_column(LenMatrix, "sample_R")
    Len = Prop[c("sample_R","Lentisphaeraceae_Prop")]
    Len <- merge(Len,meta, by.x="sample_R", by.y="sample.ID")
    Len <- Len[,c(3,2)]
    colnames(Len)[1] <- "sample_R"
    LenMatrix_final = merge(LenMatrix, Len,
                            by = "sample_R", all.x = TRUE,
                            sort = TRUE, no.dups = FALSE)
    dim(LenMatrix_final)
    LenMatrix_final = LenMatrix_final[,c(2:242)]
    write.csv(LenMatrix_final, file = "Lentisphaeraceae_Tau_Matrix.csv", row.names = FALSE)
    
    
  ##Oligoflexaceae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Oligoflexaceae = All_top_sorted_genes[,c("Oligoflexaceae_Prop","sums")]
    Oligoflexaceae = Oligoflexaceae[complete.cases(Oligoflexaceae), ]
    Oligoflexaceae = setDT(Oligoflexaceae, keep.rownames = TRUE)[]
    Oligoflexaceae = Oligoflexaceae[,c(1)]
    write.csv(Oligoflexaceae, "Oligoflexaceae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    OligoMatrix = merge(Oligoflexaceae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    OligoMatrix = OligoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    OligoMatrix = as.data.frame(t(OligoMatrix))
    OligoMatrix <- tibble::rownames_to_column(OligoMatrix, "sample_R")
    Oligo = Prop[c("sample_R","Oligoflexaceae_Prop")]
    Oligo <- merge(Oligo,meta, by.x="sample_R", by.y="sample.ID")
    Oligo <- Oligo[,c(3,2)]
    colnames(Oligo)[1] <- "sample_R"
    OligoMatrix_final = merge(OligoMatrix, Oligo,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(OligoMatrix_final)
    OligoMatrix_final = OligoMatrix_final[,c(2:278)]
    write.csv(OligoMatrix_final, file = "Oligoflexaceae_Tau_Matrix.csv", row.names = FALSE)
  
  
  
  
  
  
  