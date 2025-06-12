##this code prepares our virus data for correlating with our ortholog data##
##last updated LEF 6/12/25

##load necessary packages##
library(janitor)
library(dplyr)
library(tibble)
library(tximport)
library(DESeq2)

##read in virus classifcation data and format##
class=read.csv("Virus_classification_530_cleaned.csv")
  ##format to remove text##
  class$Domain=gsub(".*_","",class$Domain)
  class$domain_score=gsub(".*_","",class$domain_score)
  class$Realm=gsub(".*_","",class$Realm)
  class$realm_score=gsub(".*_","",class$realm_score)
  class$Kingdom=gsub(".*_","",class$Kingdom)
  class$kingdom_score=gsub(".*_","",class$kingdom_score)
  class$Phylum=gsub(".*_","",class$Phylum)
  class$phylum_score=gsub(".*_","",class$phylum_score)
  class$Class=gsub(".*_","",class$Class)
  class$class_score=gsub(".*_","",class$class_score)
  class$Order=gsub(".*_","",class$Order)
  class$order_score=gsub(".*_","",class$order_score)
  class$Family=gsub(".*_","",class$Family)
  class$family_score=gsub(".*_","",class$family_score)
  ##make numeric##
  class$domain_score=as.numeric(class$domain_score)
  class$realm_score=as.numeric(class$realm_score)
  class$kingdom_score=as.numeric(class$kingdom_score)
  class$phylum_score=as.numeric(class$phylum_score)
  class$class_score=as.numeric(class$class_score)
  class$order_score=as.numeric(class$order_score)
  class$family_score=as.numeric(class$family_score)
  ##select those that are confidently viruses
  clean=subset(class, domain_score >=1)
  clean<- clean[which(clean$Order != 'U'),]
  
  ##get rid of extra columns##
  clean=clean[,c(1:21)]
  clean=clean[,c(1,6:11)]

##now we'll bring in our normalized counts##
  counts=read.csv("Counts_TPMnorm.csv")
  ##good counts
    goodcounts=merge(counts,clean, by.x="X",by.y="sequence.name")
    goodcounts=goodcounts[,c(1:133,139)]  
  ##sum up at the level of Order##
    summed=goodcounts %>% 
      group_by(Order) %>%
      summarise(across(starts_with('DT'), sum))
    
    
##and finally we can create our matrix for correlation analyses##
    counts=t(summed)
    counts<-counts %>%
      row_to_names(row_number = 1)
    counts=as.data.frame(counts)
    counts <- tibble::rownames_to_column(counts, "sample_ID")
    ##input meta##
    meta=read.csv("metadata_matchvirus.csv")  
    meta=meta[,c(1:2)] 
    ##merge
    counts_final=merge(counts,meta, by="sample_ID")    
    counts_final=counts_final[,c(445,2:444)] 
    write.csv(counts,"virus_norm_counts_grouped.csv")
    
    
    ##now we process and normalize RNAseq data##
    ##lets make our count matrixes for each species##
    meta=read.csv("metadata_matchvirus.csv")
    cnat <- subset(meta, species == "CNAT")
    mcav <- subset(meta, species == "MCAV")
    ofav <- subset(meta, species == "OFAV")
    ofra <- subset(meta, species == "OFRA")
    
    ##mcav ##
      ##transcript to genes list- generated in original ortho file##
      tx2genes=read.csv("mcav_orthogroup_tx2gene.csv")
      head(tx2genes)
      ##files names##
      files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/mcav_quants", mcav$salmon_ID, "quant.sf")
      names(files) = paste0("sample",35:65)
      all(file.exists(files))
      
      ##import##
      txi_mcav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
   
    ##cnat ##
      ##transcript to genes list##
      tx2genes=read.csv("cnat_orthogroup_tx2gene.csv")
      head(tx2genes)
      ##files names##
      files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/cnat_quants", cnat$salmon_ID, "quant.sf")
      names(files) = paste0("sample", 1:34)
      all(file.exists(files))
      
      ##import##
      txi_cnat = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    
    ##ofra ##
      ##transcript to genes list##
      tx2genes=read.csv("ofra_orthogroup_tx2gene.csv")
      head(tx2genes)
      ##files names##
      files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofra_quants", ofra$salmon_ID, "quant.sf")
      names(files) = paste0("sample", 100:129)
      all(file.exists(files))
      
      ##import##
      txi_ofra = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    ##ofav ##
      ##transcript to genes list##
      tx2genes=read.csv("ofav_orthogroup_tx2gene.csv")
      head(tx2genes)
      ##files names##
      files <- file.path("/Users/e_f232/Library/CloudStorage/OneDrive-TexasStateUniversity/SymbiommunityLab/FDEP_omic/ortholog_analysis/ofav_quants", ofav$salmon_ID, "quant.sf")
      names(files) = paste("sample", 66:99)
      all(file.exists(files))
      
      ##import##
      txi_ofav = tximport(files, type = "salmon", tx2gene = tx2genes, countsFromAbundance = "lengthScaledTPM")
    
    ##normalization with DESeq
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
    
    ##combining it all together
    reads=merge(vst_cnat,vst_mcav, by="ortho")
    reads2=merge(reads,vst_ofav, by="ortho")  
    final_reads=merge(reads2, vst_ofra, by="ortho")  
    
    ##formatting it to make the matrix
    final_reads_merge=t(final_reads)
    final_reads_merge=final_reads_merge %>%
      row_to_names(row_number = 1)
    final_reads_merge=as.data.frame(final_reads_merge)
    final_reads_merge <- tibble::rownames_to_column(final_reads_merge, "sample_R")
    final_reads_merge$sample_R <- gsub("sample ", "sample", final_reads_merge$sample_R)
    write.csv(final_reads_merge, "normalized_reads_corr.csv", row.names=FALSE)
    
    #final merge to make the matrix
    final_matrix <- merge(counts_final,final_reads_merge, by="sample_R")
    final_matrix=final_matrix %>% remove_rownames %>% column_to_rownames(var="sample_R")  
    write.csv(final_matrix, "final_correlation_matrix_order.csv")
    