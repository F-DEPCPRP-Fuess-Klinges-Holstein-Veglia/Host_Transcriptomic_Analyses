##this code is for making graphs displaying trends of our identified response markers for our FDEP ortholog analysis
##last updated LEF 6/12/25

##load in necessary pacakges##
library(ggVennDiagram)
library(ggplot2)
library(janitor)
library(tibble)
library(reshape2)
library(ggpubr)
 
 
 ##epidemic peroid only
   ##load in our first set of orthos of interest- downregulated response markers
   GOI=read.csv("GOI_down_response.csv")
 
   #input 2022 meta data
   meta=read.csv("allspec_meta_2022.csv")
   meta <- meta[which(meta$tis != 'H'),]
   meta=meta[,c(3:4,7,9,13)]
 
 
 ##input our normalized reads for each species
   reads_mcav=read.csv("mcav_normalized_post.csv")
   reads_cnat=read.csv("cnat_normalized_post.csv")
   reads_ofav=read.csv("ofav_normalized_post.csv")
   reads_ofra=read.csv("ofra_normalized_post.csv")
   
   ##create one dataframe of all normalized reads
   reads=merge(reads_cnat,reads_mcav, by="X")
   int=merge(reads,reads_ofav, by="X")
   final_reads=merge(int,reads_ofra, by="X")  
 
   ##select out just our orthos of interest and melt as needed
   GOI_reads=merge(final_reads, GOI, by.x="X", by.y="ortho")
   GOI_data=t(GOI_reads)  
   GOI_data=GOI_data %>%
     row_to_names(row_number = 1)
   GOI_data=as.data.frame(GOI_data)
   GOI_data <- tibble::rownames_to_column(GOI_data, "sample_R")
   ##melt it##
   GOI_data=melt(GOI_data, id="sample_R")
 
   ##merge it with our meta data##
   graph=merge(meta,GOI_data,by="sample_R")
   graph$value=as.numeric(graph$value)
   graph=merge(graph,GOI, by.x="variable", by.y="ortho")
 
   ##order factors for graphing
   graph$tis <- factor(graph$tis , levels = c("AH", "DM"))
   graph$sample_month <- factor(graph$sample_month , levels = c("June", "August"))
   graph$gene <- factor(graph$gene , levels = c("AA2AR","AA2BR.1","AA2BR.2","ATS7","CY24B",
                                                "ERAP1","NFAT5","TNF15","KLH20","RN216",
                                                "S22A5","SETD2","TXD12","XERD"))
 
   ##graph it##
   ##graph it##
   ggplot(graph, aes(x=sample_month, y=value, fill=tis)) + 
     scale_fill_manual(values = c("AH" = "#5a7bac", "DM"="#ffcde6")) +
     geom_boxplot()+facet_wrap(~gene, scales="free") +
     theme_classic() +theme(legend.position="bottom")
 
   
   ##load in our second set of orthos of interest- up regulated response markers
   GOI=read.csv("GOI_up_response.csv")
   
   #input 2022 meta data
   meta=read.csv("allspec_meta_2022.csv")
   meta <- meta[which(meta$tis != 'H'),]
   meta=meta[,c(3:4,7,9,13)]
   
   
   ##input our normalized reads for each species
   reads_mcav=read.csv("mcav_normalized_post.csv")
   reads_cnat=read.csv("cnat_normalized_post.csv")
   reads_ofav=read.csv("ofav_normalized_post.csv")
   reads_ofra=read.csv("ofra_normalized_post.csv")
   
   ##create one dataframe of all normalized reads
   reads=merge(reads_cnat,reads_mcav, by="X")
   int=merge(reads,reads_ofav, by="X")
   final_reads=merge(int,reads_ofra, by="X")  
   
   ##select out just our orthos of interest and melt as needed
   GOI_reads=merge(final_reads, GOI, by.x="X", by.y="ortho")
   GOI_data=t(GOI_reads)  
   GOI_data=GOI_data %>%
     row_to_names(row_number = 1)
   GOI_data=as.data.frame(GOI_data)
   GOI_data <- tibble::rownames_to_column(GOI_data, "sample_R")
   ##melt it##
   GOI_data=melt(GOI_data, id="sample_R")
   
   ##merge it with our meta data##
   graph=merge(meta,GOI_data,by="sample_R")
   graph$value=as.numeric(graph$value)
   graph=merge(graph,GOI, by.x="variable", by.y="ortho")
   
   ##order factors for graphing
   graph$tis <- factor(graph$tis , levels = c("AH", "DM"))
   graph$sample_month <- factor(graph$sample_month , levels = c("June", "August"))
   graph$gene <- factor(graph$gene , levels = c("SAL.1","SAL.2","SML","CHCH2","CY24B",
                                                "NCF2","PDK3","PSB5","SESN1","STOX1","Y3154"))
   
   ##graph it##
   ##graph it##
   ggplot(graph, aes(x=sample_month, y=value, fill=tis)) + 
     scale_fill_manual(values = c("AH" = "#5a7bac", "DM"="#ffcde6")) +
     geom_boxplot()+facet_wrap(~gene, scales="free") +
     theme_classic() +theme(legend.position="bottom")
   
 

