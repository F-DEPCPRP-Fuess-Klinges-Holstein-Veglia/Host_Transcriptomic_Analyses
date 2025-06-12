##this code prepares our microbe data for correlating with our ortholog data##
##last updated LEF 6/12/25

##packages
library(tidyverse)
library(data.table)
library(janitor)
library(dplyr)
library(tibble)
library(DESeq2)
library(tximport)

##start by collapsing families##
##prior to this step we used BB_edit to reduce taxonomy to just family##
##the regex used was as follows: 
#Find: ,[^,^_]*+_[^,^_]*+_[^,^_]*+_[^,^_]*+_([^,^_]*+)_[^,^_]*+_[^,^_]*+
#Return: ,\1
##we needed to manipulate the first sample name to do this (added '_' to the front); was removed prior to saving

  ##data read in & processed##
  test=read.csv("drto_otus_test.csv")
  test=test %>% remove_rownames %>% column_to_rownames(var="X")
  ##transpose##
  test_1<-t(test)
  ##remove NA##
  test_2=test_1[complete.cases(test_1), ]
  ##set up for aggregating##
  test_2=as.data.frame(test_2)
  test_2$Tax=as.factor(test_2$Tax)
  ##make numeric##
  cols = c(2:134);    
  test_2[,cols] = apply(test_2[,cols], 2, function(x) as.numeric(as.character(x)));
  ##aggregate as factors##
  Summed = aggregate(. ~ Tax, data=test_2, FUN=sum)
  names(Summed)
  Summed=Summed %>% remove_rownames %>% column_to_rownames(var="Tax")
  ##check##
  Summed$Yersiniaceae
  
  ##now let's iterate over columns to create a proportion value##
  Summed=as.data.frame(t(Summed))
  Summed$Total=rowSums(Summed)
  ##need to get rid of a few columns cause their names are weird/they aren't highly abundant anyways##
  Summed=Summed[,c(4:394)]
  Summed <- tibble::rownames_to_column(Summed, "sample.ID")
  summed_2=setDT(Summed)[, paste0(names(Summed)[2:391],"_Prop") := lapply(.SD, `/`, Summed$Total), .SDcols = A21b:Yersiniaceae]
  summed_3=summed_2[,c(1,393:782)]
  ##magic!! now we filter for just those families with an average prop greater than .001##
  prop_t=t(summed_3)
  prop_t=prop_t %>%
    row_to_names(row_number = 1)
  prop_t=as.data.frame(prop_t)
  prop_t<- mutate_all(prop_t, function(x) as.numeric(as.character(x)))
  #write.csv(prop_t,"families_prop.csv")
  
  ##now we process for combining with read matrices##
  ##mostly need to pair sample names##
  meta=read.csv("metadata_matchasv.csv")
  meta=meta[1:2]  
  prop_t=t(prop_t)  
  prop_t=as.data.frame(prop_t)
  ##row names to column
  prop_t <- tibble::rownames_to_column(prop_t, "sample.ID")
  ##merge##
  asvdat=merge(meta,prop_t, by="sample.ID")
  asvdat=asvdat[,c(2:392)]

##now we process and normalize RNAseq data##
  ##lets make our count matrixes for each species##
  meta=read.csv("metadata_matchasv.csv")
  cnat <- subset(meta, species == "CNAT")
  mcav <- subset(meta, species == "MCAV")
  ofav <- subset(meta, species == "OFAV")
  ofra <- subset(meta, species == "OFRA")
  ##mcav ##
    ##transcript to genes list##
    tx2genes=read.csv("mcav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/shilling/FDEP/ortholog_analysis/mcav_quants", mcav$salmon_ID, "quant.sf")
    names(files) = paste0("sample",35:65)
    all(file.exists(files))
    
    ##import##
    txi_mcav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
  
  ##cnat ##
    ##transcript to genes list##
    tx2genes=read.csv("cnat_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/shilling/FDEP/ortholog_analysis/cnat_quants", cnat$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 1:34)
    all(file.exists(files))
    
    ##import##
    txi_cnat = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
  
    ##ofra ##
    ##transcript to genes list##
    tx2genes=read.csv("ofra_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/shilling/FDEP/ortholog_analysis/ofra_quants", ofra$salmon_ID, "quant.sf")
    names(files) = paste0("sample", 100:129)
    all(file.exists(files))
    
    ##import##
    txi_ofra = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    
  
  ##ofav ##
    ##transcript to genes list##
    tx2genes=read.csv("ofav_orthogroup_tx2gene.csv")
    head(tx2genes)
    ##files names##
    files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/shilling/FDEP/ortholog_analysis/ofav_quants", ofav$salmon_ID, "quant.sf")
    names(files) = paste("sample",66:99)
    all(file.exists(files))
    
    ##import##
    txi_ofav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    
  ##normalization
  ##mcav
    dds_mcav <- DESeqDataSetFromTximport(txi_mcav, mcav, ~sample_month)
    ##normalizing reads##
    dds_mcav <- estimateSizeFactors(dds_mcav)         
    dds_mcav <- estimateDispersions(dds_mcav)
    vst_mcav <- getVarianceStabilizedData(dds_mcav)
    vst_mcav=as.data.frame(vst_mcav)
    vst_mcav <- tibble::rownames_to_column(vst_mcav, "ortho")
  ##cnat
    dds_cnat <- DESeqDataSetFromTximport(txi_cnat, cnat, ~sample_month)
    ##normalizing reads##
    dds_cnat <- estimateSizeFactors(dds_cnat)         
    dds_cnat <- estimateDispersions(dds_cnat)
    vst_cnat <- getVarianceStabilizedData(dds_cnat)
    vst_cnat=as.data.frame(vst_cnat)
    vst_cnat <- tibble::rownames_to_column(vst_cnat, "ortho")
  ##ofav
    dds_ofav <- DESeqDataSetFromTximport(txi_ofav, ofav, ~sample_month)
    ##normalizing reads##
    dds_ofav <- estimateSizeFactors(dds_ofav)         
    dds_ofav <- estimateDispersions(dds_ofav)
    vst_ofav <- getVarianceStabilizedData(dds_ofav)
    vst_ofav=as.data.frame(vst_ofav)
    vst_ofav <- tibble::rownames_to_column(vst_ofav, "ortho")
  ##ofra
    dds_ofra <- DESeqDataSetFromTximport(txi_ofra, ofra, ~sample_month)
    ##normalizing reads##
    dds_ofra <- estimateSizeFactors(dds_ofra)         
    dds_ofra <- estimateDispersions(dds_ofra)
    vst_ofra <- getVarianceStabilizedData(dds_ofra)
    vst_ofra=as.data.frame(vst_ofra)
    vst_ofra <- tibble::rownames_to_column(vst_ofra, "ortho")
  
  ##combining into one frame
  reads=merge(vst_cnat,vst_mcav, by="ortho")
  reads2=merge(reads,vst_ofav, by="ortho")  
  final_reads=merge(reads2, vst_ofra, by="ortho")  
  
  ##formatting to match other data
  final_reads_merge=t(final_reads)
  final_reads_merge=final_reads_merge %>%
    row_to_names(row_number = 1)
  final_reads_merge=as.data.frame(final_reads_merge)
  final_reads_merge <- tibble::rownames_to_column(final_reads_merge, "sample_R")
  write.csv(final_reads_merge, "normalized_reads_corr.csv", row.names=FALSE)
  
  #merge to create our final correlation matrix
  final_matrix <- merge(asvdat,final_reads_merge, by="sample_R")
  final_matrix=final_matrix %>% remove_rownames %>% column_to_rownames(var="sample_R")  
  write.csv(final_matrix, "final_correlation_matrix_families.csv")
  lapply(final_matrix,class)
  final_matrix<- mutate_all(final_matrix, function(x) as.numeric(as.character(x)))
  