##this code parses the output of our correlation code which generates pvalues for all pariwise comparisons
##specifically we determine our top 5% most correlated orders and genes and create new matrices for each so we can calucate tau values for downstream analyses
##last updated LEF 6/12/25

##load in necessary packages
library(janitor)
library(tibble)
library(data.table)

##read in p value matrix
  p=read.csv("Corr_results_pval_virus.csv")
  colnames(p)
  ##reduce down to just grab those columns corresponding to a viral order/group##
  p_reduced=p[,c(1:442)] 
  colnames(p_reduced)
  rownames(p_reduced)
  ##reduce down to just grap those rows corresponding to orthologs
  p_reduced=p_reduced[c(444:7905),]
  ##set non sig values to NA for processing
  p_reduced[p_reduced>0.05]=NA

##and now we look for those orders in the top 5%##
  P_fam=as.data.frame(t(p_reduced))
  ##start by replacing non-sig (n/as) with 1000000)
  P_fam[is.na(P_fam)] <- 1000000
  ##create a column for row sums##
  P_fam$sums = rowSums(P_fam)
  ##order by sums##
  P_fam_sorted <- P_fam[order(P_fam$sums),] 
  P_fam_sorted$sums
  ##pull out the top 5% of ASV (roughly 22)##
  P_fam_top_sorted = P_fam_sorted[1:22,]
  row.names(P_fam_top_sorted)
  is.na(P_fam_top_sorted) <- P_fam_top_sorted == 1000000
  ##generate a list of top orders for reference##
  P_fam_familes = P_fam_top_sorted[1]
  write.csv(P_fam_familes, file="top5orders.csv")

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
  ##pull out the top 5% of genes (roughly 373)-383 is best for tie##
  All_sorted_genes$sums
  All_top_sorted_genes = All_sorted_genes[1:383,]
  is.na(All_top_sorted_genes) <- All_top_sorted_genes == 1000000

##now parse out lists of significant genes in top 5% of both family and genes##
##and create matrices for figuring out tau values##
  ##processing files##
  ##read in normalized reads and format##
  reads=read.csv("normalized_reads_corr.csv")
  reads=t(reads)
  reads=reads %>%
    row_to_names(row_number = 1)
  reads=as.data.frame(reads)
  reads <- tibble::rownames_to_column(reads, "rn")
  
  ##read in count matrix data and format##
  Prop = read.csv("virus_norm_counts_grouped.csv", check.names = FALSE)
  ##read in sample names##
  meta=read.csv("metadata_matchvirus.csv")
  meta=meta[,c(1:2)]
  ##fix it to be correctly formatted##
  Prop=merge(Prop,meta, by="sample_ID")
  Prop=Prop[,c(446,2:444)]
  

  ##now use the dataframes we've generated to create lists of sig genes
  ##also create matrices for getting tau values of associatons between each top order and associated top genes.
  
  ##Streptococcus.phage.phiNJ2
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Streptococcus.phage.phiNJ2 = All_top_sorted_genes[,c("Streptococcus.phage.phiNJ2","sums")]
    Streptococcus.phage.phiNJ2 = Streptococcus.phage.phiNJ2[complete.cases(Streptococcus.phage.phiNJ2), ]
    Streptococcus.phage.phiNJ2 = setDT(Streptococcus.phage.phiNJ2, keep.rownames = TRUE)[]
    Streptococcus.phage.phiNJ2 = Streptococcus.phage.phiNJ2[,c(1)]
    write.csv(Streptococcus.phage.phiNJ2, "StreptococcusphagephiNJ2_Sig_Genes.csv", row.names = FALSE)
  
    ##then generate your matrix for determining tau values##
    StreNJMatrix = merge(Streptococcus.phage.phiNJ2, reads,
                      by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    StreNJMatrix = StreNJMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    StreNJMatrix = as.data.frame(t(StreNJMatrix))
    StreNJMatrix <- tibble::rownames_to_column(StreNJMatrix, "sample_R")
    StreNJ = Prop[c("sample_R","Streptococcus phage phiNJ2")]
    StreNJMatrix_final = merge(StreNJMatrix, StreNJ,
                       by = "sample_R", all.x = TRUE,
                       sort = TRUE, no.dups = FALSE)
    dim(StreNJMatrix_final)
    StreNJMatrix_final = StreNJMatrix_final[,c(2:385)]
    write.csv(StreNJMatrix_final, file = "StreptococcusphagephiNJ2_Tau_Matrix.csv", row.names = FALSE)
  
  

  ##Godonkavirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Godonkavirus = All_top_sorted_genes[,c("Godonkavirus","sums")]
    Godonkavirus = Godonkavirus[complete.cases(Godonkavirus), ]
    Godonkavirus = setDT(Godonkavirus, keep.rownames = TRUE)[]
    Godonkavirus = Godonkavirus[,c(1)]
    write.csv(Godonkavirus, "Godonkavirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    GodoMatrix = merge(Godonkavirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    GodoMatrix = GodoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    GodoMatrix = as.data.frame(t(GodoMatrix))
    GodoMatrix <- tibble::rownames_to_column(GodoMatrix, "sample_R")
    Godo = Prop[c("sample_R","Godonkavirus")]
    GodoMatrix_final = merge(GodoMatrix, Godo,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(GodoMatrix_final)
    GodoMatrix_final = GodoMatrix_final[,c(2:330)]
    write.csv(GodoMatrix_final, file = "Godonkavirus_Tau_Matrix.csv", row.names = FALSE)
  
  
  ##Clostridium.phage.phiCT19406C
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Clostridium.phage.phiCT19406C = All_top_sorted_genes[,c("Clostridium.phage.phiCT19406C","sums")]
    Clostridium.phage.phiCT19406C = Clostridium.phage.phiCT19406C[complete.cases(Clostridium.phage.phiCT19406C), ]
    Clostridium.phage.phiCT19406C = setDT(Clostridium.phage.phiCT19406C, keep.rownames = TRUE)[]
    Clostridium.phage.phiCT19406C = Clostridium.phage.phiCT19406C[,c(1)]
    write.csv(Clostridium.phage.phiCT19406C, "ClostridiumphagephiCT19406C_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ClostMatrix = merge(Clostridium.phage.phiCT19406C, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    ClostMatrix = ClostMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ClostMatrix = as.data.frame(t(ClostMatrix))
    ClostMatrix <- tibble::rownames_to_column(ClostMatrix, "sample_R")
    Clost = Prop[c("sample_R","Clostridium phage phiCT19406C")]
    ClostMatrix_final = merge(ClostMatrix, Clost,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(ClostMatrix_final)
    ClostMatrix_final = ClostMatrix_final[,c(2:385)]
    write.csv(ClostMatrix_final, file = "ClostridiumphagephiCT19406C_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Decurrovirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Decurrovirus = All_top_sorted_genes[,c("Decurrovirus","sums")]
    Decurrovirus = Decurrovirus[complete.cases(Decurrovirus), ]
    Decurrovirus = setDT(Decurrovirus, keep.rownames = TRUE)[]
    Decurrovirus = Decurrovirus[,c(1)]
    write.csv(Decurrovirus, "Decurrovirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    DecurMatrix = merge(Decurrovirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    DecurMatrix = DecurMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    DecurMatrix = as.data.frame(t(DecurMatrix))
    DecurMatrix <- tibble::rownames_to_column(DecurMatrix, "sample_R")
    Decur = Prop[c("sample_R","Decurrovirus")]
    DecurMatrix_final = merge(DecurMatrix, Decur,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(DecurMatrix_final)
    DecurMatrix_final = DecurMatrix_final[,c(2:385)]
    write.csv(DecurMatrix_final, file = "Decurrovirus_Tau_Matrix.csv", row.names = FALSE)
    
    
  ##Bingvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Bingvirus = All_top_sorted_genes[,c("Bingvirus","sums")]
    Bingvirus = Bingvirus[complete.cases(Bingvirus), ]
    Bingvirus = setDT(Bingvirus, keep.rownames = TRUE)[]
    Bingvirus = Bingvirus[,c(1)]
    write.csv(Bingvirus, "Bingvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    BingMatrix = merge(Bingvirus, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    BingMatrix = BingMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    BingMatrix = as.data.frame(t(BingMatrix))
    BingMatrix <- tibble::rownames_to_column(BingMatrix, "sample_R")
    Bing = Prop[c("sample_R","Bingvirus")]
    BingMatrix_final = merge(BingMatrix, Bing,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(BingMatrix_final)
    BingMatrix_final = BingMatrix_final[,c(2:385)]
    write.csv(BingMatrix_final, file = "Bingvirus_Tau_Matrix.csv", row.names = FALSE)
    
    
  ##Pseudoalteromonas.phage.H103
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Pseudoalteromonas.phage.H103 = All_top_sorted_genes[,c("Pseudoalteromonas.phage.H103","sums")]
    Pseudoalteromonas.phage.H103 = Pseudoalteromonas.phage.H103[complete.cases(Pseudoalteromonas.phage.H103), ]
    Pseudoalteromonas.phage.H103 = setDT(Pseudoalteromonas.phage.H103, keep.rownames = TRUE)[]
    Pseudoalteromonas.phage.H103 = Pseudoalteromonas.phage.H103[,c(1)]
    write.csv(Pseudoalteromonas.phage.H103, "PseudoalteromonasphageH103_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    PseudoMatrix = merge(Pseudoalteromonas.phage.H103, reads,
                          by = "rn", all.x = TRUE,
                          sort = TRUE,  no.dups = FALSE,)
    PseudoMatrix = PseudoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    PseudoMatrix = as.data.frame(t(PseudoMatrix))
    PseudoMatrix <- tibble::rownames_to_column(PseudoMatrix, "sample_R")
    Pseudo = Prop[c("sample_R","Pseudoalteromonas phage H103")]
    PseudoMatrix_final = merge(PseudoMatrix, Pseudo,
                                by = "sample_R", all.x = TRUE,
                                sort = TRUE, no.dups = FALSE)
    dim(PseudoMatrix_final)
    PseudoMatrix_final = PseudoMatrix_final[,c(2:385)]
    write.csv(PseudoMatrix_final, file = "PseudoalteromonasphageH103_Tau_Matrix.csv", row.names = FALSE)
    
  
  
  
  ##Brochothrix.phage.BL3
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Brochothrix.phage.BL3 = All_top_sorted_genes[,c("Brochothrix.phage.BL3","sums")]
    Brochothrix.phage.BL3 = Brochothrix.phage.BL3[complete.cases(Brochothrix.phage.BL3), ]
    Brochothrix.phage.BL3 = setDT(Brochothrix.phage.BL3, keep.rownames = TRUE)[]
    Brochothrix.phage.BL3 = Brochothrix.phage.BL3[,c(1)]
    write.csv(Brochothrix.phage.BL3, "BrochothrixphageBL3_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    BL3Matrix = merge(Brochothrix.phage.BL3, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    BL3Matrix = BL3Matrix %>% remove_rownames %>% column_to_rownames(var="rn")
    BL3Matrix = as.data.frame(t(BL3Matrix))
    BL3Matrix <- tibble::rownames_to_column(BL3Matrix, "sample_R")
    BL3 = Prop[c("sample_R","Brochothrix phage BL3")]
    BL3Matrix_final = merge(BL3Matrix, BL3,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(BL3Matrix_final)
    BL3Matrix_final = BL3Matrix_final[,c(2:385)]
    write.csv(BL3Matrix_final, file = "BrochothrixphageBL3_Tau_Matrix.csv", row.names = FALSE)
    
    
  ##Burrovirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Burrovirus = All_top_sorted_genes[,c("Burrovirus","sums")]
    Burrovirus = Burrovirus[complete.cases(Burrovirus), ]
    Burrovirus = setDT(Burrovirus, keep.rownames = TRUE)[]
    Burrovirus = Burrovirus[,c(1)]
    write.csv(Burrovirus, "Burrovirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    BurrMatrix = merge(Burrovirus, reads,
                         by = "rn", all.x = TRUE,
                         sort = TRUE,  no.dups = FALSE,)
    BurrMatrix = BurrMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    BurrMatrix = as.data.frame(t(BurrMatrix))
    BurrMatrix <- tibble::rownames_to_column(BurrMatrix, "sample_R")
    Burr = Prop[c("sample_R","Burrovirus")]
    BurrMatrix_final = merge(BurrMatrix, Burr,
                               by = "sample_R", all.x = TRUE,
                               sort = TRUE, no.dups = FALSE)
    dim(BurrMatrix_final)
    BurrMatrix_final = BurrMatrix_final[,c(2:338)]
    write.csv(BurrMatrix_final, file = "Burrovirus_Tau_Matrix.csv", row.names = FALSE)
    
  
  ##Carmenvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Carmenvirus = All_top_sorted_genes[,c("Carmenvirus","sums")]
    Carmenvirus = Carmenvirus[complete.cases(Carmenvirus), ]
    Carmenvirus = setDT(Carmenvirus, keep.rownames = TRUE)[]
    Carmenvirus = Carmenvirus[,c(1)]
    write.csv(Carmenvirus, "Carmenvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    CarmMatrix = merge(Carmenvirus, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    CarmMatrix = CarmMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    CarmMatrix = as.data.frame(t(CarmMatrix))
    CarmMatrix <- tibble::rownames_to_column(CarmMatrix, "sample_R")
    Carm = Prop[c("sample_R","Carmenvirus")]
    CarmMatrix_final = merge(CarmMatrix, Carm,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(CarmMatrix_final)
    CarmMatrix_final = CarmMatrix_final[,c(2:385)]
    write.csv(CarmMatrix_final, file = "Carmenvirus_Tau_Matrix.csv", row.names = FALSE)
    
    
  
  ##Halcyonevirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Halcyonevirus = All_top_sorted_genes[,c("Halcyonevirus","sums")]
    Halcyonevirus = Halcyonevirus[complete.cases(Halcyonevirus), ]
    Halcyonevirus = setDT(Halcyonevirus, keep.rownames = TRUE)[]
    Halcyonevirus = Halcyonevirus[,c(1)]
    write.csv(Halcyonevirus, "Halcyonevirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    HalcyMatrix = merge(Halcyonevirus, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    HalcyMatrix = HalcyMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    HalcyMatrix = as.data.frame(t(HalcyMatrix))
    HalcyMatrix <- tibble::rownames_to_column(HalcyMatrix, "sample_R")
    Halcy = Prop[c("sample_R","Halcyonevirus")]
    HalcyMatrix_final = merge(HalcyMatrix, Halcy,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(HalcyMatrix_final)
    HalcyMatrix_final = HalcyMatrix_final[,c(2:385)]
    write.csv(HalcyMatrix_final, file = "Halcyonevirus_Tau_Matrix.csv", row.names = FALSE)
    
    
  
  ##Klebsiella.phage.phiKO2
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Klebsiella.phage.phiKO2 = All_top_sorted_genes[,c("Klebsiella.phage.phiKO2","sums")]
    Klebsiella.phage.phiKO2 = Klebsiella.phage.phiKO2[complete.cases(Klebsiella.phage.phiKO2), ]
    Klebsiella.phage.phiKO2 = setDT(Klebsiella.phage.phiKO2, keep.rownames = TRUE)[]
    Klebsiella.phage.phiKO2 = Klebsiella.phage.phiKO2[,c(1)]
    write.csv(Klebsiella.phage.phiKO2, "KlebsiellaphagephiKO2_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    KO2Matrix = merge(Klebsiella.phage.phiKO2, reads,
                          by = "rn", all.x = TRUE,
                          sort = TRUE,  no.dups = FALSE,)
    KO2Matrix = KO2Matrix %>% remove_rownames %>% column_to_rownames(var="rn")
    KO2Matrix = as.data.frame(t(KO2Matrix))
    KO2Matrix <- tibble::rownames_to_column(KO2Matrix, "sample_R")
    KO2 = Prop[c("sample_R","Klebsiella phage phiKO2")]
    KO2Matrix_final = merge(KO2Matrix, KO2,
                                by = "sample_R", all.x = TRUE,
                                sort = TRUE, no.dups = FALSE)
    dim(KO2Matrix_final)
    KO2Matrix_final = KO2Matrix_final[,c(2:242)]
    write.csv(KO2Matrix_final, file = "KlebsiellaphagephiKO2_Tau_Matrix.csv", row.names = FALSE)
  
    
    
  ##Uetakevirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Uetakevirus = All_top_sorted_genes[,c("Uetakevirus","sums")]
    Uetakevirus = Uetakevirus[complete.cases(Uetakevirus), ]
    Uetakevirus = setDT(Uetakevirus, keep.rownames = TRUE)[]
    Uetakevirus = Uetakevirus[,c(1)]
    write.csv(Uetakevirus, "Uetakevirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    UetakMatrix = merge(Uetakevirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    UetakMatrix = UetakMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    UetakMatrix = as.data.frame(t(UetakMatrix))
    UetakMatrix <- tibble::rownames_to_column(UetakMatrix, "sample_R")
    Uetak = Prop[c("sample_R","Uetakevirus")]
    UetakMatrix_final = merge(UetakMatrix, Uetak,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(UetakMatrix_final)
    UetakMatrix_final = UetakMatrix_final[,c(2:320)]
    write.csv(UetakMatrix_final, file = "Uetakevirus_Tau_Matrix.csv", row.names = FALSE)  
  
  
  ##Vertoviridae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Vertoviridae = All_top_sorted_genes[,c("Vertoviridae","sums")]
    Vertoviridae = Vertoviridae[complete.cases(Vertoviridae), ]
    Vertoviridae = setDT(Vertoviridae, keep.rownames = TRUE)[]
    Vertoviridae = Vertoviridae[,c(1)]
    write.csv(Vertoviridae, "Vertoviridae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    VertMatrix = merge(Vertoviridae, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    VertMatrix = VertMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    VertMatrix = as.data.frame(t(VertMatrix))
    VertMatrix <- tibble::rownames_to_column(VertMatrix, "sample_R")
    Vert = Prop[c("sample_R","Vertoviridae")]
    VertMatrix_final = merge(VertMatrix, Vert,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(VertMatrix_final)
    VertMatrix_final = VertMatrix_final[,c(2:319)]
    write.csv(VertMatrix_final, file = "Vertoviridae_Tau_Matrix.csv", row.names = FALSE)

  ##Zobellviridae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Zobellviridae = All_top_sorted_genes[,c("Zobellviridae","sums")]
    Zobellviridae = Zobellviridae[complete.cases(Zobellviridae), ]
    Zobellviridae = setDT(Zobellviridae, keep.rownames = TRUE)[]
    Zobellviridae = Zobellviridae[,c(1)]
    write.csv(Zobellviridae, "Zobellviridae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    ZobeMatrix = merge(Zobellviridae, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    ZobeMatrix = ZobeMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    ZobeMatrix = as.data.frame(t(ZobeMatrix))
    ZobeMatrix <- tibble::rownames_to_column(ZobeMatrix, "sample_R")
    Zobe = Prop[c("sample_R","Zobellviridae")]
    ZobeMatrix_final = merge(ZobeMatrix, Zobe,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(ZobeMatrix_final)
    ZobeMatrix_final = ZobeMatrix_final[,c(2:385)]
    write.csv(ZobeMatrix_final, file = "Zobellviridae_Tau_Matrix.csv", row.names = FALSE)  
    
  
  ##Lacusarxvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Lacusarxvirus = All_top_sorted_genes[,c("Lacusarxvirus","sums")]
    Lacusarxvirus = Lacusarxvirus[complete.cases(Lacusarxvirus), ]
    Lacusarxvirus = setDT(Lacusarxvirus, keep.rownames = TRUE)[]
    Lacusarxvirus = Lacusarxvirus[,c(1)]
    write.csv(Lacusarxvirus, "Lacusarxvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    LacuMatrix = merge(Lacusarxvirus, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    LacuMatrix = LacuMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    LacuMatrix = as.data.frame(t(LacuMatrix))
    LacuMatrix <- tibble::rownames_to_column(LacuMatrix, "sample_R")
    Lacu = Prop[c("sample_R","Lacusarxvirus")]
    LacuMatrix_final = merge(LacuMatrix, Lacu,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(LacuMatrix_final)
    LacuMatrix_final = LacuMatrix_final[,c(2:186)]
    write.csv(LacuMatrix_final, file = "Lacusarxvirus_Tau_Matrix.csv", row.names = FALSE) 
    
  
  ##Lightbulbvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Lightbulbvirus = All_top_sorted_genes[,c("Lightbulbvirus","sums")]
    Lightbulbvirus = Lightbulbvirus[complete.cases(Lightbulbvirus), ]
    Lightbulbvirus = setDT(Lightbulbvirus, keep.rownames = TRUE)[]
    Lightbulbvirus = Lightbulbvirus[,c(1)]
    write.csv(Lightbulbvirus, "Lightbulbvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    LightMatrix = merge(Lightbulbvirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    LightMatrix = LightMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    LightMatrix = as.data.frame(t(LightMatrix))
    LightMatrix <- tibble::rownames_to_column(LightMatrix, "sample_R")
    Light = Prop[c("sample_R","Lightbulbvirus")]
    LightMatrix_final = merge(LightMatrix, Light,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(LightMatrix_final)
    LightMatrix_final = LightMatrix_final[,c(2:357)]
    write.csv(LightMatrix_final, file = "Lightbulbvirus_Tau_Matrix.csv", row.names = FALSE) 

  ##Nonagvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Nonagvirus = All_top_sorted_genes[,c("Nonagvirus","sums")]
    Nonagvirus = Nonagvirus[complete.cases(Nonagvirus), ]
    Nonagvirus = setDT(Nonagvirus, keep.rownames = TRUE)[]
    Nonagvirus = Nonagvirus[,c(1)]
    write.csv(Nonagvirus, "Nonagvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    NonaMatrix = merge(Nonagvirus, reads,
                       by = "rn", all.x = TRUE,
                       sort = TRUE,  no.dups = FALSE,)
    NonaMatrix = NonaMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    NonaMatrix = as.data.frame(t(NonaMatrix))
    NonaMatrix <- tibble::rownames_to_column(NonaMatrix, "sample_R")
    Nona = Prop[c("sample_R","Nonagvirus")]
    NonaMatrix_final = merge(NonaMatrix, Nona,
                             by = "sample_R", all.x = TRUE,
                             sort = TRUE, no.dups = FALSE)
    dim(NonaMatrix_final)
    NonaMatrix_final = NonaMatrix_final[,c(2:346)]
    write.csv(NonaMatrix_final, file = "Nonagvirus_Tau_Matrix.csv", row.names = FALSE)  
    
  
  ##Pseudomonas.phage.JBD44
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Pseudomonas.phage.JBD44 = All_top_sorted_genes[,c("Pseudomonas.phage.JBD44","sums")]
    Pseudomonas.phage.JBD44 = Pseudomonas.phage.JBD44[complete.cases(Pseudomonas.phage.JBD44), ]
    Pseudomonas.phage.JBD44 = setDT(Pseudomonas.phage.JBD44, keep.rownames = TRUE)[]
    Pseudomonas.phage.JBD44 = Pseudomonas.phage.JBD44[,c(1)]
    write.csv(Pseudomonas.phage.JBD44, "PseudomonasphageJBD44_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    JBD44Matrix = merge(Pseudomonas.phage.JBD44, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    JBD44Matrix = JBD44Matrix %>% remove_rownames %>% column_to_rownames(var="rn")
    JBD44Matrix = as.data.frame(t(JBD44Matrix))
    JBD44Matrix <- tibble::rownames_to_column(JBD44Matrix, "sample_R")
    JBD44 = Prop[c("sample_R","Pseudomonas phage JBD44")]
    JBD44Matrix_final = merge(JBD44Matrix, JBD44,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(JBD44Matrix_final)
    JBD44Matrix_final = JBD44Matrix_final[,c(2:381)]
    write.csv(JBD44Matrix_final, file = "PseudomonasphageJBD44_Tau_Matrix.csv", row.names = FALSE) 
    
  
  ##Schitoviridae
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Schitoviridae = All_top_sorted_genes[,c("Schitoviridae","sums")]
    Schitoviridae = Schitoviridae[complete.cases(Schitoviridae), ]
    Schitoviridae = setDT(Schitoviridae, keep.rownames = TRUE)[]
    Schitoviridae = Schitoviridae[,c(1)]
    write.csv(Schitoviridae, "Schitoviridae_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    SchitMatrix = merge(Schitoviridae, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    SchitMatrix = SchitMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    SchitMatrix = as.data.frame(t(SchitMatrix))
    SchitMatrix <- tibble::rownames_to_column(SchitMatrix, "sample_R")
    Schit = Prop[c("sample_R","Schitoviridae")]
    SchitMatrix_final = merge(SchitMatrix, Schit,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(SchitMatrix_final)
    SchitMatrix_final = SchitMatrix_final[,c(2:297)]
    write.csv(SchitMatrix_final, file = "Schitoviridae_Tau_Matrix.csv", row.names = FALSE) 
  
    
  ##Spizizenvirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Spizizenvirus = All_top_sorted_genes[,c("Spizizenvirus","sums")]
    Spizizenvirus = Spizizenvirus[complete.cases(Spizizenvirus), ]
    Spizizenvirus = setDT(Spizizenvirus, keep.rownames = TRUE)[]
    Spizizenvirus = Spizizenvirus[,c(1)]
    write.csv(Spizizenvirus, "Spizizenvirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    SpizMatrix = merge(Spizizenvirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    SpizMatrix = SpizMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    SpizMatrix = as.data.frame(t(SpizMatrix))
    SpizMatrix <- tibble::rownames_to_column(SpizMatrix, "sample_R")
    Spiz = Prop[c("sample_R","Spizizenvirus")]
    SpizMatrix_final = merge(SpizMatrix, Spiz,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(SpizMatrix_final)
    SpizMatrix_final = SpizMatrix_final[,c(2:385)]
    write.csv(SpizMatrix_final, file = "Spizizenvirus_Tau_Matrix.csv", row.names = FALSE)  
    
  
  ##Streptococcus.phage.phiARI0746
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Streptococcus.phage.phiARI0746 = All_top_sorted_genes[,c("Streptococcus.phage.phiARI0746","sums")]
    Streptococcus.phage.phiARI0746 = Streptococcus.phage.phiARI0746[complete.cases(Streptococcus.phage.phiARI0746), ]
    Streptococcus.phage.phiARI0746 = setDT(Streptococcus.phage.phiARI0746, keep.rownames = TRUE)[]
    Streptococcus.phage.phiARI0746 = Streptococcus.phage.phiARI0746[,c(1)]
    write.csv(Streptococcus.phage.phiARI0746, "StreptococcusphagephiARI0746_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    StrepMatrix = merge(Streptococcus.phage.phiARI0746, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    StrepMatrix = StrepMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    StrepMatrix = as.data.frame(t(StrepMatrix))
    StrepMatrix <- tibble::rownames_to_column(StrepMatrix, "sample_R")
    Strep = Prop[c("sample_R","Streptococcus phage phiARI0746")]
    StrepMatrix_final = merge(StrepMatrix, Strep,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(StrepMatrix_final)
    StrepMatrix_final = StrepMatrix_final[,c(2:385)]
    write.csv(StrepMatrix_final, file = "StreptococcusphagephiARI0746_Tau_Matrix.csv", row.names = FALSE) 
    
  
  ##Vividuovirus
    ##start with outputting a general list of which of the top 5% genes are correlated to this order
    Vividuovirus = All_top_sorted_genes[,c("Vividuovirus","sums")]
    Vividuovirus = Vividuovirus[complete.cases(Vividuovirus), ]
    Vividuovirus = setDT(Vividuovirus, keep.rownames = TRUE)[]
    Vividuovirus = Vividuovirus[,c(1)]
    write.csv(Vividuovirus, "Vividuovirus_Sig_Genes.csv", row.names = FALSE)
    
    ##then generate your matrix for determining tau values##
    VividMatrix = merge(Vividuovirus, reads,
                        by = "rn", all.x = TRUE,
                        sort = TRUE,  no.dups = FALSE,)
    VividMatrix = VividMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
    VividMatrix = as.data.frame(t(VividMatrix))
    VividMatrix <- tibble::rownames_to_column(VividMatrix, "sample_R")
    Vivid = Prop[c("sample_R","Vividuovirus")]
    VividMatrix_final = merge(VividMatrix, Vivid,
                              by = "sample_R", all.x = TRUE,
                              sort = TRUE, no.dups = FALSE)
    dim(VividMatrix_final)
    VividMatrix_final = VividMatrix_final[,c(2:385)]
    write.csv(VividMatrix_final, file = "Vividuovirus_Tau_Matrix.csv", row.names = FALSE) 
    