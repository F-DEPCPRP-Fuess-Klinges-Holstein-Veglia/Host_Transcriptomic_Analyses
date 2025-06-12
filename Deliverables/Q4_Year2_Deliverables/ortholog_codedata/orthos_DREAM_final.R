##this code is for ortholog analysis of differential gene expression associated with our FDEP project, using DREAM
##last updated LEF 6/12/25

##load libraries
library(tximport)
library(DESeq2)
library(tibble)
library(tidyr)
library(limma)
library(edgeR)
library(variancePartition)

##to start with need to make new tx2gene for our orthologs for each species.
  ##start by loading our master sc ortholog datasheet
  sc_ortho=read.csv("SC_orthologs.csv",header=TRUE)
  mcav=sc_ortho[c(5,1)]
  mcav$mcav_ref_proteome <- gsub(".{3}$","",mcav$mcav_ref_proteome)
  cnat=sc_ortho[,c(4,1)]  
  cnat$cnat_ref_proteome <- gsub(".{3}$","",cnat$cnat_ref_proteome)
  ofav=sc_ortho[,c(6,1)] 
  ofav$ofav_ref_proteome <- gsub(".{3}$","",ofav$ofav_ref_proteome)
  ofra=sc_ortho[,c(7,1)]  
  ofra$ofra_ref_proteome <- gsub(".{3}$","",ofra$ofra_ref_proteome)
  write.csv(mcav, "mcav_orthogroup_tx2gene.csv",row.names=FALSE)  
  write.csv(cnat, "cnat_orthogroup_tx2gene.csv",row.names=FALSE)  
  write.csv(ofav, "ofav_orthogroup_tx2gene.csv",row.names=FALSE)    
  write.csv(ofra, "ofra_orthogroup_tx2gene.csv",row.names=FALSE)   
  
##now we'll run some models- starting with predictors from pre-disease data##
  ##lets make our count matrixes for each species##
  meta=read.csv("allspec_meta_pre.csv")
  cnat <- subset(meta, species == "CNAT")
  mcav <- subset(meta, species == "MCAV")
  ofav <- subset(meta, species == "OFAV")
  ofra <- subset(meta, species == "OFRA")
  ##mcav ##
    ##transcript to genes list##
    tx2genes=read.csv("mcav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/mcav_quants", mcav$salmon_ID, "quant.sf")
    names(files) = paste0("sample",10:16)
    all(file.exists(files))
    
    ##import##
    txi_mcav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    head(txi_mcav$counts)
    mcav_counts <- txi_mcav$counts
    colnames(mcav_counts) <- paste0("mcav",1:7)
    mcav_counts <- as.data.frame(mcav_counts)
    #write.csv(mcav_counts, file = "mcav_counts_orthos_pre.csv",quote = FALSE)
    
  ##cnat ##
    ##transcript to genes list##
    tx2genes=read.csv("cnat_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/cnat_quants", cnat$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 1:9)
    all(file.exists(files))
    
    ##import##
    txi_cnat = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    #head(txi$counts)
    cnat_counts <- txi_cnat$counts
    colnames(cnat_counts) <- paste0("cnat",1:9)
    cnat_counts <- as.data.frame(cnat_counts)
    #write.csv(cnat_counts, file = "cnat_counts_orthos_pre.csv",quote = FALSE)
    
  ##ofra ##
    ##transcript to genes list##
    tx2genes=read.csv("ofra_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofra_quants", ofra$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 27:35)
    all(file.exists(files))
    
    ##import##
    txi_ofra = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    head(txi_ofra$counts)
    ofra_counts <- txi_ofra$counts
    colnames(ofra_counts) <- paste0("ofra",1:9)
    ofra_counts <- as.data.frame(ofra_counts)
    #write.csv(ofra_counts, file = "ofra_counts_orthos_pre.csv",quote = FALSE)
    
  ##ofav ##
    ##transcript to genes list##
    tx2genes=read.csv("ofav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofav_quants", ofav$salmon_ID, "quant.sf")
    names(files) = paste("sample",18:26)
    all(file.exists(files))
    
    ##import##
    txi_ofav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    head(txi_ofav$counts)
    ofav_counts <- txi_ofav$counts
    colnames(ofav_counts) <- paste0("ofav",1:9)
    ofav_counts <- as.data.frame(ofav_counts)
    #write.csv(ofav_counts, file = "ofav_counts_orthos_pre.csv",quote = FALSE)
    
  ##normalization- if needed for downstream figs and analysis
  ##mcav
    #dds_mcav <- DESeqDataSetFromTximport(txi_mcav, mcav, ~category)
    ##normalizing reads##
    #dds_mcav <- estimateSizeFactors(dds_mcav)         
    #dds_mcav <- estimateDispersions(dds_mcav)
    #vst_mcav <- getVarianceStabilizedData(dds_mcav)
    #write.csv(vst_mcav, file = "mcav_normalized_pre.csv")
  ##cnat
    #dds_cnat <- DESeqDataSetFromTximport(txi_cnat, cnat, ~category)
    ##normalizing reads##
    #dds_cnat <- estimateSizeFactors(dds_cnat)         
    #dds_cnat <- estimateDispersions(dds_cnat)
    #vst_cnat <- getVarianceStabilizedData(dds_cnat)
    #write.csv(vst_cnat, file = "cnat_normalized_pre.csv")
  ##ofav
    #dds_ofav <- DESeqDataSetFromTximport(txi_ofav, ofav, ~category)
    ##normalizing reads##
    #dds_ofav <- estimateSizeFactors(dds_ofav)         
    #dds_ofav <- estimateDispersions(dds_ofav)
    #vst_ofav <- getVarianceStabilizedData(dds_ofav)
    #write.csv(vst_ofav, file = "ofav_normalized_pre.csv")
  ##ofra
    #dds_ofra <- DESeqDataSetFromTximport(txi_ofra, ofra, ~category)
    ##normalizing reads##
    #dds_ofra <- estimateSizeFactors(dds_ofra)         
    #dds_ofra <- estimateDispersions(dds_ofra)
    #vst_ofra <- getVarianceStabilizedData(dds_ofra)
    #write.csv(vst_ofra, file = "ofra_normalized_pre.csv")
  
  
##now we merge together all dataframes##
  mcav_counts <- tibble::rownames_to_column(mcav_counts, "ortho")
  cnat_counts <- tibble::rownames_to_column(cnat_counts, "ortho")
  ofra_counts <- tibble::rownames_to_column(ofra_counts, "ortho")
  ofav_counts <- tibble::rownames_to_column(ofav_counts, "ortho")
  int=merge(cnat_counts,mcav_counts, by="ortho")
  int2=merge(int,ofav_counts, by="ortho")
  final_counts=merge(int2,ofra_counts, by="ortho")  
  final_meta=meta[,c(5,4,6:7)]
  final_meta=final_meta %>% remove_rownames %>% column_to_rownames(var="sample_ortho")
  
  
  ##dream analysis-interaction##
    dge=DGEList(final_counts)
    ##normalize##
    dge=calcNormFactors(dge)
    ##set model for interaction of species and disease reistance category##
    test<- ~0+species*category+(1 | site)
    vobjDream <- voomWithDreamWeights(dge,test, final_meta)
    ##set up your contrasts##
    L=makeContrastsDream(test, final_meta, contrasts = c( MCAV_OFAVint="categorysusceptible:speciesMCAV-categorysusceptible:speciesOFAV",
                                                         MCAV_OFRAint="categorysusceptible:speciesMCAV-categorysusceptible:speciesOFRA",
                                                         OFAV_OFRAint="categorysusceptible:speciesOFAV-categorysusceptible:speciesOFRA"))
    ##check contrasts##
    plotContrasts(L)
    ##run model
    fitmm<-dream(vobjDream, test, final_meta, L)
    fitmm<-eBayes(fitmm)
    ##for reference
    colnames(fitmm)
    ##process and parse data##
    
    ##mcav vs. ofav interaction term##
    MCAVOFAVint=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="MCAV_OFAVint")
    write.csv(MCAVOFAVint, "ortho_pre_int_results_MCAVOFAV.csv")
    MCAVOFAVint=MCAVOFAVint[,c(1:2,6)]
    colnames(MCAVOFAVint)<- c("ortho","LFC_MCAVOFAV","padj_MCAVOFAV")
    
    ##mcav vs. ofra interaction term##
    MCAVOFRAint=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="MCAV_OFRAint")
    write.csv(MCAVOFRAint, "ortho_pre_int_results_MCAVOFRA.csv")
    MCAVOFRAint=MCAVOFRAint[,c(1:2,6)]
    colnames(MCAVOFRAint)<- c("ortho","LFC_MCAVOFRA","padj_MCAVOFRA")
    
    ##ofav vs. ofra interaction term##
    OFAVOFRAint=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="OFAV_OFRAint")
    write.csv(OFAVOFRAint, "ortho_pre_int_results_OFAVOFRA.csv")
    OFAVOFRAint=OFAVOFRAint[,c(1:2,6)]
    colnames(OFAVOFRAint)<- c("ortho","LFC_OFAVOFRA","padj_OFAVOFRA")
    
    ##cnat vs. mcav interaction term##
    CNATMCAVint=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="categorysusceptible:speciesMCAV")
    write.csv(CNATMCAVint, "ortho_pre_int_results_CNATMCAV.csv")
    CNATMCAVint=CNATMCAVint[,c(1:2,6)]
    colnames(CNATMCAVint)<- c("ortho","LFC_CNATMCAV","padj_CNATMCAV")
    
    ##cnat vs. ofav interaction term##
    CNATOFAVint=topTable(fitmm, number = 7462, adjust = "BH", p.value = 1, coef="categorysusceptible:speciesOFAV")
    write.csv(CNATOFAVint, "ortho_pre_int_results_CNATOFAV.csv")
    CNATOFAVint=CNATOFAVint[,c(1:2,6)]
    colnames(CNATOFAVint)<- c("ortho","LFC_CNATOFAV","padj_CNATOFAV")
    
    ##cnat vs. ofra interaction term##
    CNATOFRAint=topTable(fitmm, number = 7462, adjust = "BH", p.value = 1, coef="categorysusceptible:speciesOFRA")
    write.csv(CNATOFRAint, "ortho_pre_int_results_CNATOFRA.csv")
    CNATOFRAint=CNATOFRAint[,c(1:2,6)]
    colnames(CNATOFRAint)<- c("ortho","LFC_CNATOFRA","padj_CNATOFRA")
    
    ##create one summary file containing LFC and padj for each interaction term
    test=merge(MCAVOFAVint,MCAVOFRAint, by="ortho")
    test=merge(test,OFAVOFRAint, by="ortho")
    test=merge(test,CNATMCAVint, by="ortho")
    test=merge(test,CNATOFAVint, by="ortho")
    test=merge(test,CNATOFRAint, by="ortho")
    write.csv(test,"ortho_pre_int_results_sum.csv")

  
  ##now do the main effects model##
  dge=DGEList(final_counts)
  ##normalize##
  dge=calcNormFactors(dge)
  ##set up model##
  test<- ~0+category+species+(1 | site)
  vobjDream <- voomWithDreamWeights(dge,test, final_meta)
  ##set up your contrasts##
  L=makeContrastsDream(test, final_meta, contrasts = c(sus="categoryresistant -categorysusceptible"))
  ##check contrasts
  plotContrasts(L)
  ##run##
  fitmm<-dream(vobjDream, test, final_meta,L)
  fitmm<-eBayes(fitmm)
  colnames(fitmm)
  ##pull main effects of resistance category##
  sus=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="sus")
  write.csv(sus, "ortho_pre_main_results_sus.csv", row.names = FALSE)
  
  
##next let's model during disease to find our response and epidemic markers##
  ##lets make our count matrixes for each species (same general code as before)##
  meta=read.csv("allspec_meta_2022.csv")
  cnat <- subset(meta, species == "CNAT")
  mcav <- subset(meta, species == "MCAV")
  ofav <- subset(meta, species == "OFAV")
  ofra <- subset(meta, species == "OFRA")
  ##mcav ##
    ##transcript to genes list##
    tx2genes=read.csv("mcav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/mcav_quants", mcav$salmon_ID, "quant.sf")
    names(files) = paste0("sample",22:38)
    all(file.exists(files))
    
    ##import##
    txi_mcav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    mcav_counts <- txi_mcav$counts
    colnames(mcav_counts) <- paste0("mcav",1:17)
    mcav_counts <- as.data.frame(mcav_counts)
    write.csv(mcav_counts, file = "mcav_counts_orthos_post.csv",quote = FALSE)
  
  ##cnat ##
    ##transcript to genes list##
    tx2genes=read.csv("cnat_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/cnat_quants", cnat$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 1:21)
    all(file.exists(files))
    
    ##import##
    txi_cnat = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    #head(txi$counts)
    cnat_counts <- txi_cnat$counts
    colnames(cnat_counts) <- paste0("cnat",1:21)
    cnat_counts <- as.data.frame(cnat_counts)
    write.csv(cnat_counts, file = "cnat_counts_orthos_post.csv",quote = FALSE)
  
  ##ofra ##
    ##transcript to genes list##
    tx2genes=read.csv("ofra_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofra_quants", ofra$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 59:76)
    all(file.exists(files))
    
    ##import##
    txi_ofra = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    ofra_counts <- txi_ofra$counts
    colnames(ofra_counts) <- paste0("ofra",1:18)
    ofra_counts <- as.data.frame(ofra_counts)
    write.csv(ofra_counts, file = "ofra_counts_orthos_post.csv",quote = FALSE)
  
  ##ofav ##
    ##transcript to genes list##
    tx2genes=read.csv("ofav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofav_quants", ofav$salmon_ID, "quant.sf")
    names(files) = paste("sample",39:58)
    all(file.exists(files))
    
    ##import##
    txi_ofav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    ofav_counts <- txi_ofav$counts
    colnames(ofav_counts) <- paste0("ofav",1:20)
    ofav_counts <- as.data.frame(ofav_counts)
    write.csv(ofav_counts, file = "ofav_counts_orthos_post.csv",quote = FALSE)
  
  ##normalization- if needed for downstream plots/analyses##
  ##mcav
    #dds_mcav <- DESeqDataSetFromTximport(txi_mcav, mcav, ~month_tis_spec)
    ##normalizing reads##
    #dds_mcav <- estimateSizeFactors(dds_mcav)         
    #dds_mcav <- estimateDispersions(dds_mcav)
    #vst_mcav <- getVarianceStabilizedData(dds_mcav)
    #write.csv(vst_mcav, file = "mcav_normalized_post.csv")
  ##cnat
    #dds_cnat <- DESeqDataSetFromTximport(txi_cnat, cnat, ~month_tis_spec)
    ##normalizing reads##
    #dds_cnat <- estimateSizeFactors(dds_cnat)         
    #dds_cnat <- estimateDispersions(dds_cnat)
    #vst_cnat <- getVarianceStabilizedData(dds_cnat)
    #write.csv(vst_cnat, file = "cnat_normalized_post.csv")
  ##ofav
    #dds_ofav <- DESeqDataSetFromTximport(txi_ofav, ofav, ~month_tis_spec)
    ##normalizing reads##
    #dds_ofav <- estimateSizeFactors(dds_ofav)         
    #dds_ofav <- estimateDispersions(dds_ofav)
    #vst_ofav <- getVarianceStabilizedData(dds_ofav)
    #write.csv(vst_ofav, file = "ofav_normalized_post.csv")
  ##ofra
    #dds_ofra <- DESeqDataSetFromTximport(txi_ofra, ofra, ~month_tis_spec)
    ##normalizing reads##
    #dds_ofra <- estimateSizeFactors(dds_ofra)         
    #dds_ofra <- estimateDispersions(dds_ofra)
    #vst_ofra <- getVarianceStabilizedData(dds_ofra)
    #write.csv(vst_ofra, file = "ofra_normalized_post.csv")
  
  ##now we merge together all dataframes##
  mcav_counts <- tibble::rownames_to_column(mcav_counts, "ortho")
  cnat_counts <- tibble::rownames_to_column(cnat_counts, "ortho")
  ofra_counts <- tibble::rownames_to_column(ofra_counts, "ortho")
  ofav_counts <- tibble::rownames_to_column(ofav_counts, "ortho")
  int=merge(cnat_counts,mcav_counts, by="ortho")
  int2=merge(int,ofav_counts, by="ortho")
  final_counts=merge(int2,ofra_counts, by="ortho")  
  final_meta=meta[,c(5,4,6,9,11,12)]
  final_meta=final_meta %>% remove_rownames %>% column_to_rownames(var="sample_ortho")
  ##set factors up##
  final_meta$species=as.factor(final_meta$species)
  final_meta$site=as.factor(final_meta$site)
  final_meta$colony_ID=as.factor(final_meta$colony_ID)
  final_meta$month_tis_type=as.factor(final_meta$month_tis_type)
  final_meta$month_tis_spec=as.factor(final_meta$month_tis_spec)
  
  
  ##dream analysis-first one to get first set of interaction terms##
  dge=DGEList(final_counts)
  ##normalize##
  dge=calcNormFactors(dge)
  ##set model; it's a little confusing but our contrast allows us to look at species specific effects##
  test<- ~0+month_tis_spec+(1 | site)+(1 | colony_ID)
  ##run it##
  vobjDream <- voomWithDreamWeights(dge,test, final_meta)
  colnames(vobjDream)
  ##set up your contrasts##
  L=makeContrastsDream(test, final_meta, contrasts = c(J_AHvH_MCAVvOFAV="(month_tis_specJune_AH_MCAV-month_tis_specJune_H_MCAV)-
                                                       (month_tis_specJune_AH_OFAV-month_tis_specJune_H_OFAV)",
                                                       J_AHvH_MCAVvOFRA="(month_tis_specJune_AH_MCAV-month_tis_specJune_H_MCAV)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_H_OFRA)",
                                                       J_AHvH_OFAVvOFRA="(month_tis_specJune_AH_OFAV-month_tis_specJune_H_OFAV)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_H_OFRA)",
                                                       J_AHvH_CNATvMCAV="(month_tis_specJune_AH_CNAT-month_tis_specJune_H_CNAT)-
                                                       (month_tis_specJune_AH_MCAV-month_tis_specJune_H_MCAV)",
                                                       J_AHvH_CNATvOFAV="(month_tis_specJune_AH_CNAT-month_tis_specJune_H_CNAT)-
                                                       (month_tis_specJune_AH_OFAV-month_tis_specJune_H_OFAV)",
                                                       J_AHvH_CNATvOFRA="(month_tis_specJune_AH_CNAT-month_tis_specJune_H_CNAT)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_H_OFRA)",
                                                       J_AHvDM_MCAVvOFAV="(month_tis_specJune_AH_MCAV-month_tis_specJune_DM_MCAV)-
                                                       (month_tis_specJune_AH_OFAV-month_tis_specJune_DM_OFAV)",
                                                       J_AHvDM_MCAVvOFRA="(month_tis_specJune_AH_MCAV-month_tis_specJune_DM_MCAV)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_DM_OFRA)",
                                                       J_AHvDM_OFAVvOFRA="(month_tis_specJune_AH_OFAV-month_tis_specJune_DM_OFAV)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_DM_OFRA)",
                                                       J_AHvDM_CNATvMCAV="(month_tis_specJune_AH_CNAT-month_tis_specJune_DM_CNAT)-
                                                       (month_tis_specJune_AH_MCAV-month_tis_specJune_DM_MCAV)",
                                                       J_AHvDM_CNATvOFAV="(month_tis_specJune_AH_CNAT-month_tis_specJune_DM_CNAT)-
                                                       (month_tis_specJune_AH_OFAV-month_tis_specJune_H_OFAV)",
                                                       J_AHvDM_CNATvOFRA="(month_tis_specJune_AH_CNAT-month_tis_specJune_DM_CNAT)-
                                                       (month_tis_specJune_AH_OFRA-month_tis_specJune_DM_OFRA)",
                                                       A_AHvH_MCAVvOFAV="(month_tis_specAug_AH_MCAV-month_tis_specAug_H_MCAV)-
                                                       (month_tis_specAug_AH_OFAV-month_tis_specAug_H_OFAV)",
                                                       A_AHvH_MCAVvOFRA="(month_tis_specAug_AH_MCAV-month_tis_specAug_H_MCAV)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_H_OFRA)",
                                                       A_AHvH_OFAVvOFRA="(month_tis_specAug_AH_OFAV-month_tis_specAug_H_OFAV)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_H_OFRA)",
                                                       A_AHvH_CNATvMCAV="(month_tis_specAug_AH_CNAT-month_tis_specAug_H_CNAT)-
                                                       (month_tis_specAug_AH_MCAV-month_tis_specAug_H_MCAV)",
                                                       A_AHvH_CNATvOFAV="(month_tis_specAug_AH_CNAT-month_tis_specAug_H_CNAT)-
                                                       (month_tis_specAug_AH_OFAV-month_tis_specAug_H_OFAV)",
                                                       A_AHvH_CNATvOFRA="(month_tis_specAug_AH_CNAT-month_tis_specAug_H_CNAT)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_H_OFRA)",
                                                       A_AHvDM_MCAVvOFAV="(month_tis_specAug_AH_MCAV-month_tis_specAug_DM_MCAV)-
                                                       (month_tis_specAug_AH_OFAV-month_tis_specAug_DM_OFAV)",
                                                       A_AHvDM_MCAVvOFRA="(month_tis_specAug_AH_MCAV-month_tis_specAug_DM_MCAV)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_DM_OFRA)",
                                                       A_AHvDM_OFAVvOFRA="(month_tis_specAug_AH_OFAV-month_tis_specAug_DM_OFAV)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_DM_OFRA)",
                                                       A_AHvDM_CNATvMCAV="(month_tis_specAug_AH_CNAT-month_tis_specAug_DM_CNAT)-
                                                       (month_tis_specAug_AH_MCAV-month_tis_specAug_DM_MCAV)",
                                                       A_AHvDM_CNATvOFAV="(month_tis_specAug_AH_CNAT-month_tis_specAug_DM_CNAT)-
                                                       (month_tis_specAug_AH_OFAV-month_tis_specAug_DM_OFAV)",
                                                       A_AHvDM_CNATvOFRA="(month_tis_specAug_AH_CNAT-month_tis_specAug_DM_CNAT)-
                                                       (month_tis_specAug_AH_OFRA-month_tis_specAug_DM_OFRA)"))
                                                    
  ##check contrasts##                                             
  plotContrasts(L)
  ##run it
  fitmm<-dream(vobjDream, test, final_meta,L)
  fitmm<-eBayes(fitmm)
  ##for reference
  colnames(fitmm)
  
  ##parse results for each contrast of interest##
  J_AHvH_MCAVvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_MCAVvOFAV")
  write.csv(J_AHvH_MCAVvOFAV, "ortho_epi_int_results_J_AHH_MCOF.csv")
  J_AHvH_MCAVvOFAV=J_AHvH_MCAVvOFAV[,c(1:2,6)]
  colnames(J_AHvH_MCAVvOFAV)<- c("ortho","LFC_J_AHH_MCOF","padj_J_AHH_MCOF")
  
  J_AHvH_MCAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_MCAVvOFRA")
  write.csv(J_AHvH_MCAVvOFRA, "ortho_epi_int_results_J_AHH_MCOR.csv")
  J_AHvH_MCAVvOFRA=J_AHvH_MCAVvOFRA[,c(1:2,6)]
  colnames(J_AHvH_MCAVvOFRA)<- c("ortho","LFC_J_AHH_MCOR","padj_J_AHH_MCOR")
  
  J_AHvH_OFAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_OFAVvOFRA")
  write.csv(J_AHvH_OFAVvOFRA, "ortho_epi_int_results_J_AHH_OFOR.csv")
  J_AHvH_OFAVvOFRA=J_AHvH_OFAVvOFRA[,c(1:2,6)]
  colnames(J_AHvH_OFAVvOFRA)<- c("ortho","LFC_J_AHH_OFOR","padj_J_AHH_OFOR")
  
  J_AHvH_CNATvMCAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_CNATvMCAV")
  write.csv(J_AHvH_CNATvMCAV, "ortho_epi_int_results_J_AHH_CNMC.csv")
  J_AHvH_CNATvMCAV=J_AHvH_CNATvMCAV[,c(1:2,6)]
  colnames(J_AHvH_CNATvMCAV)<- c("ortho","LFC_J_AHH_CNMC","padj_J_AHH_CNMC")
  
  J_AHvH_CNATvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_CNATvOFAV")
  write.csv(J_AHvH_CNATvOFAV, "ortho_epi_int_results_J_AHH_CNOF.csv")
  J_AHvH_CNATvOFAV=J_AHvH_CNATvOFAV[,c(1:2,6)]
  colnames(J_AHvH_CNATvOFAV)<- c("ortho","LFC_J_AHH_CNMC","padj_J_AHH_CNOF")
  
  J_AHvH_CNATvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvH_CNATvOFRA")
  write.csv(J_AHvH_CNATvOFRA, "ortho_epi_int_results_J_AHH_CNOR.csv")
  J_AHvH_CNATvOFRA=J_AHvH_CNATvOFRA[,c(1:2,6)]
  colnames(J_AHvH_CNATvOFRA)<- c("ortho","LFC_J_AHH_CNMC","padj_J_AHH_CNOR")
  
  A_AHvH_MCAVvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_MCAVvOFAV")
  write.csv(A_AHvH_MCAVvOFAV, "ortho_epi_int_results_A_AHH_MCOF.csv")
  A_AHvH_MCAVvOFAV=A_AHvH_MCAVvOFAV[,c(1:2,6)]
  colnames(A_AHvH_MCAVvOFAV)<- c("ortho","LFC_A_AHH_MCOF","padj_A_AHH_MCOF")
  
  A_AHvH_MCAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_MCAVvOFRA")
  write.csv(A_AHvH_MCAVvOFRA, "ortho_epi_int_results_A_AHH_MCOR.csv")
  A_AHvH_MCAVvOFRA=A_AHvH_MCAVvOFRA[,c(1:2,6)]
  colnames(A_AHvH_MCAVvOFRA)<- c("ortho","LFC_A_AHH_MCOR","padj_A_AHH_MCOR")
  
  A_AHvH_OFAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_OFAVvOFRA")
  write.csv(A_AHvH_OFAVvOFRA, "ortho_epi_int_results_A_AHH_OFOR.csv")
  A_AHvH_OFAVvOFRA=A_AHvH_OFAVvOFRA[,c(1:2,6)]
  colnames(A_AHvH_OFAVvOFRA)<- c("ortho","LFC_A_AHH_OFOR","padj_A_AHH_OFOR")
  
  A_AHvH_CNATvMCAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_CNATvMCAV")
  write.csv(A_AHvH_CNATvMCAV, "ortho_epi_int_results_A_AHH_CNMC.csv")
  A_AHvH_CNATvMCAV=A_AHvH_CNATvMCAV[,c(1:2,6)]
  colnames(A_AHvH_CNATvMCAV)<- c("ortho","LFC_A_AHH_CNMC","padj_A_AHH_CNMC")
  
  A_AHvH_CNATvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_CNATvOFAV")
  write.csv(A_AHvH_CNATvOFAV, "ortho_epi_int_results_A_AHH_CNOF.csv")
  A_AHvH_CNATvOFAV=A_AHvH_CNATvOFAV[,c(1:2,6)]
  colnames(A_AHvH_CNATvOFAV)<- c("ortho","LFC_A_AHH_CNMC","padj_A_AHH_CNOF")
  
  A_AHvH_CNATvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvH_CNATvOFRA")
  write.csv(A_AHvH_CNATvOFRA, "ortho_epi_int_results_A_AHH_CNOR.csv")
  A_AHvH_CNATvOFRA=A_AHvH_CNATvOFRA[,c(1:2,6)]
  colnames(A_AHvH_CNATvOFRA)<- c("ortho","LFC_A_AHH_CNMC","padj_A_AHH_CNOR")
  
  J_AHvDM_MCAVvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_MCAVvOFAV")
  write.csv(J_AHvDM_MCAVvOFAV, "ortho_epi_int_results_J_AHDM_MCOF.csv")
  J_AHvDM_MCAVvOFAV=J_AHvDM_MCAVvOFAV[,c(1:2,6)]
  colnames(J_AHvDM_MCAVvOFAV)<- c("ortho","LFC_J_AHDM_MCOF","padj_J_AHDM_MCOF")
  
  J_AHvDM_MCAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_MCAVvOFRA")
  write.csv(J_AHvDM_MCAVvOFRA, "ortho_epi_int_results_J_AHDM_MCOR.csv")
  J_AHvDM_MCAVvOFRA=J_AHvDM_MCAVvOFRA[,c(1:2,6)]
  colnames(J_AHvDM_MCAVvOFRA)<- c("ortho","LFC_J_AHDM_MCOR","padj_J_AHDM_MCOR")
  
  J_AHvDM_OFAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_OFAVvOFRA")
  write.csv(J_AHvDM_OFAVvOFRA, "ortho_epi_int_results_J_AHDM_OFOR.csv")
  J_AHvDM_OFAVvOFRA=J_AHvDM_OFAVvOFRA[,c(1:2,6)]
  colnames(J_AHvDM_OFAVvOFRA)<- c("ortho","LFC_J_AHDM_OFOR","padj_J_AHDM_OFOR")
  
  J_AHvDM_CNATvMCAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_CNATvMCAV")
  write.csv(J_AHvDM_CNATvMCAV, "ortho_epi_int_results_J_AHDM_CNMC.csv")
  J_AHvDM_CNATvMCAV=J_AHvDM_CNATvMCAV[,c(1:2,6)]
  colnames(J_AHvDM_CNATvMCAV)<- c("ortho","LFC_J_AHDM_CNMC","padj_J_AHDM_CNMC")
  
  J_AHvDM_CNATvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_CNATvOFAV")
  write.csv(J_AHvDM_CNATvOFAV, "ortho_epi_int_results_J_AHDM_CNOF.csv")
  J_AHvDM_CNATvOFAV=J_AHvDM_CNATvOFAV[,c(1:2,6)]
  colnames(J_AHvDM_CNATvOFAV)<- c("ortho","LFC_J_AHDM_CNMC","padj_J_AHDM_CNOF")
  
  J_AHvDM_CNATvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="J_AHvDM_CNATvOFRA")
  write.csv(J_AHvDM_CNATvOFRA, "ortho_epi_int_results_J_AHDM_CNOR.csv")
  J_AHvDM_CNATvOFRA=J_AHvDM_CNATvOFRA[,c(1:2,6)]
  colnames(J_AHvDM_CNATvOFRA)<- c("ortho","LFC_J_AHDM_CNMC","padj_J_AHDM_CNOR")
  
  A_AHvDM_MCAVvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_MCAVvOFAV")
  write.csv(A_AHvDM_MCAVvOFAV, "ortho_epi_int_results_A_AHDM_MCOF.csv")
  A_AHvDM_MCAVvOFAV=A_AHvDM_MCAVvOFAV[,c(1:2,6)]
  colnames(A_AHvDM_MCAVvOFAV)<- c("ortho","LFC_A_AHDM_MCOF","padj_A_AHDM_MCOF")
  
  A_AHvDM_MCAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_MCAVvOFRA")
  write.csv(A_AHvDM_MCAVvOFRA, "ortho_epi_int_results_A_AHDM_MCOR.csv")
  A_AHvDM_MCAVvOFRA=A_AHvDM_MCAVvOFRA[,c(1:2,6)]
  colnames(A_AHvDM_MCAVvOFRA)<- c("ortho","LFC_A_AHDM_MCOR","padj_A_AHDM_MCOR")
  
  A_AHvDM_OFAVvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_OFAVvOFRA")
  write.csv(A_AHvDM_OFAVvOFRA, "ortho_epi_int_results_A_AHDM_OFOR.csv")
  A_AHvDM_OFAVvOFRA=A_AHvDM_OFAVvOFRA[,c(1:2,6)]
  colnames(A_AHvDM_OFAVvOFRA)<- c("ortho","LFC_A_AHDM_OFOR","padj_A_AHDM_OFOR")
  
  A_AHvDM_CNATvMCAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_CNATvMCAV")
  write.csv(A_AHvDM_CNATvMCAV, "ortho_epi_int_results_A_AHDM_CNMC.csv")
  A_AHvDM_CNATvMCAV=A_AHvDM_CNATvMCAV[,c(1:2,6)]
  colnames(A_AHvDM_CNATvMCAV)<- c("ortho","LFC_A_AHDM_CNMC","padj_A_AHDM_CNMC")
  
  A_AHvDM_CNATvOFAV=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_CNATvOFAV")
  write.csv(A_AHvDM_CNATvOFAV, "ortho_epi_int_results_A_AHDM_CNOF.csv")
  A_AHvDM_CNATvOFAV=A_AHvDM_CNATvOFAV[,c(1:2,6)]
  colnames(A_AHvDM_CNATvOFAV)<- c("ortho","LFC_A_AHDM_CNMC","padj_A_AHDM_CNOF")
  
  A_AHvDM_CNATvOFRA=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="A_AHvDM_CNATvOFRA")
  write.csv(A_AHvDM_CNATvOFRA, "ortho_epi_int_results_A_AHDM_CNOR.csv")
  A_AHvDM_CNATvOFRA=A_AHvDM_CNATvOFRA[,c(1:2,6)]
  colnames(A_AHvDM_CNATvOFRA)<- c("ortho","LFC_A_AHDM_CNMC","padj_A_AHDM_CNOR")
  
  
  ##merging to get our summary dataframe
  test=merge(J_AHvH_MCAVvOFAV,J_AHvH_MCAVvOFRA, by="ortho")
  test=merge(test,J_AHvH_OFAVvOFRA, by="ortho")
  test=merge(test,J_AHvH_CNATvMCAV, by="ortho")
  test=merge(test,J_AHvH_CNATvOFAV, by="ortho")
  test=merge(test,J_AHvH_CNATvOFRA, by="ortho")
  test=merge(test,A_AHvH_MCAVvOFAV, by="ortho")
  test=merge(test,A_AHvH_MCAVvOFRA, by="ortho")
  test=merge(test,A_AHvH_OFAVvOFRA, by="ortho")
  test=merge(test,A_AHvH_CNATvMCAV, by="ortho")
  test=merge(test,A_AHvH_CNATvOFAV, by="ortho")
  test=merge(test,A_AHvH_CNATvOFRA, by="ortho")
  test=merge(test,J_AHvDM_MCAVvOFAV, by="ortho")
  test=merge(test,J_AHvDM_MCAVvOFRA, by="ortho")
  test=merge(test,J_AHvDM_OFAVvOFRA, by="ortho")
  test=merge(test,J_AHvDM_CNATvMCAV, by="ortho")
  test=merge(test,J_AHvDM_CNATvOFAV, by="ortho")
  test=merge(test,J_AHvDM_CNATvOFRA, by="ortho")
  test=merge(test,A_AHvDM_MCAVvOFAV, by="ortho")
  test=merge(test,A_AHvDM_MCAVvOFRA, by="ortho")
  test=merge(test,A_AHvDM_OFAVvOFRA, by="ortho")
  test=merge(test,A_AHvDM_CNATvMCAV, by="ortho")
  test=merge(test,A_AHvDM_CNATvOFAV, by="ortho")
  test=merge(test,A_AHvDM_CNATvOFRA, by="ortho")
  
  ##writing out our summary dataframe
  write.csv(test,"ortho_post_int_results_sum.csv", row.names=FALSE)
  
  
  ##dream analysis-let's get the second set of contrasts (main)##
  dge=DGEList(final_counts)
  ##normalize##
  dge=calcNormFactors(dge)
  ##set model, this time looking at main effects of tissue type contrasts within months, accounting for species##
  test<- ~0+month_tis_type+species+(1 | site)+(1 | colony_ID)
  vobjDream <- voomWithDreamWeights(dge,test, final_meta)
  ##set up contrasts##
  L=makeContrastsDream(test, final_meta, contrasts = c(AugAHvDM="month_tis_typeAug_AH-month_tis_typeAug_DM",
                                                       AugAHvH="month_tis_typeAug_AH-month_tis_typeAug_H",
                                                       JuneAHvDM="month_tis_typeJune_AH-month_tis_typeJune_DM",
                                                       JuneAHvH="month_tis_typeJune_AH-month_tis_typeJune_H"))
  ##check contrasts##
  plotContrasts(L)
  ##run it##
  fitmm<-dream(vobjDream, test, final_meta,L)
  fitmm<-eBayes(fitmm)
  ##for reference
  colnames(fitmm)
  
  ##process each contrast of interest##
  AugAHvDM=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="AugAHvDM")
  write.csv(AugAHvDM, "ortho_epi_int_results_AugAHvDM.csv")
  AugAHvDM=AugAHvDM[,c(1:2,6)]
  colnames(AugAHvDM)<- c("ortho","LFC_AugAHvDM","padj_AugAHvDM")
  
  AugAHvH=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="AugAHvH")
  write.csv(AugAHvH, "ortho_epi_int_results_AugAHvH.csv")
  AugAHvH=AugAHvH[,c(1:2,6)]
  colnames(AugAHvH)<- c("ortho","LFC_AugAHvH","padj_AugAHvH")
  
  JuneAHvDM=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="JuneAHvDM")
  write.csv(JuneAHvDM, "ortho_epi_int_results_JuneAHvDM.csv")
  JuneAHvDM=JuneAHvDM[,c(1:2,6)]
  colnames(JuneAHvDM)<- c("ortho","LFC_JuneAHvDM","padj_JuneAHvDM")
  
  JuneAHvH=topTable(fitmm, number = 8000, adjust = "BH", p.value = 1, coef="JuneAHvH")
  write.csv(JuneAHvH, "ortho_epi_int_results_JuneAHvH.csv")
  JuneAHvH=JuneAHvH[,c(1:2,6)]
  colnames(JuneAHvH)<- c("ortho","LFC_JuneAHvH","padj_JuneAHvH")
  
  ##merge and write out our summary data##
  test=merge(AugAHvDM,AugAHvH, by="ortho")
  test=merge(test,JuneAHvDM, by="ortho")
  test=merge(test,JuneAHvH, by="ortho")
  write.csv(test,"ortho_epi_main_results_sum.csv", row.names=FALSE)