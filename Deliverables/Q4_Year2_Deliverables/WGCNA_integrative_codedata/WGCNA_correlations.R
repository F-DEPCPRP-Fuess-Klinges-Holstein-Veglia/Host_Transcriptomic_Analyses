##This is a general WGCNA Script to take the ortholog network we created and correlate it to traits of interest##
##in this case we're correlating to microbial family abundance (for our top 5 families), and histological traits of interest
##last updated LEF 6/12/25

##we need to load some packages
library(janitor)
library(tibble)
library(WGCNA)
library(wesanderson)
library(extrafont)

##start by format our correlation data##
  ##microbiome##
  ##read in proportion data##
  prop=read.csv("families_prop.csv")
  ##read in families of interest##
  FOI=read.csv("top5families.csv")
  ##and finally metadata
  micro_met=read.csv("metadata_matchasv.csv")
  ##select only families of interest from propdata frame##
  prop_int=merge(FOI,prop, by="X")
  prop_int=prop_int[,c(1,3:135)]
  prop_int=t(prop_int)
  prop_int=prop_int %>%
    row_to_names(row_number = 1)
  prop_int=as.data.frame(prop_int)
  ##and merge with metadata to fix names##
  prop_int <- tibble::rownames_to_column(prop_int, "sample.ID")
  prop_fin=merge(prop_int,micro_met, by="sample.ID")
  prop_fin=prop_fin[,c(21,2:20)]
  
  ##histo
  ##read in histology data##
  hist=read.csv("histodata_clean_noblank.csv")
  ##select relavent columns##
  hist=hist[,c(2,7:9)]
  ##merge with microbiome data to create our final dataframe
  traitData=merge(prop_fin,hist, by="sample_R")


## then we load up our data from the previous WGCNA script
load("WGCNANetworkConstruction.RData")
#make sure to check it#
table(moduleColors)

## Next check traitData
dim(traitData)
names(traitData)
##Next we combine the sets (WGCNA data and trait data)
datExpr=datExpr0
Samples = rownames(datExpr);
traitRows = match(Samples, traitData$sample_R);
datTraits = traitData[traitRows, 1:23];
rownames(datTraits) = traitData[traitRows, 1];
names(datTraits)
datTraits = datTraits[c(2:23)]
#check
datTraits
collectGarbage();

##Make the Matrix! This is verbatum from WGCNA tutorials##

  # Define numbers of genes and samples
  nGenes = ncol(datExpr0);
  nSamples = nrow(datExpr0);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  MEsnames=MEs
  moduleTraitCor = bicor(MEs, datTraits, use = "pairwise.complete.obs");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  #par(mar = c(6, 8.5, 3, 3));
  # Make the tombstone plot
  colors = c("#5a7bac","#dee6f6","#bd004a")
  colorRampPalette(colors)
  names(MEsnames) <- substring(names(MEs), 3)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEsnames),
                 ySymbols = names(MEs),
                 xLabelsAngle = 60,
                 colorLabels = FALSE,
                 colors =  colorRampPalette(colors)(100),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 yColorWidth = 100,
                 cex.text = 0.5,
                 cex.lab = 0.6,
                 yColorOffset = 0.01,
                 zlim = c(-1,1))
  #plotLegend = FALSE)
  
  ##Now we generate gene significance for vaculoization
  #again taken straight from tutorials
  vac = as.data.frame(datTraits$mean_symb_vac);
  names(vac) = "SymbVac"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, vac, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(vac), sep="");
  names(GSPvalue) = paste("p.GS.", names(vac), sep="");


  ##Finally we export all of this info (gene significant, module membership, etc) to a csv##
  ##this will be used to create our matrices for GOMWU analysis of interesting modules
  # again taken straight from WGCNA tutorial
  locus = names(datExpr)
  geneInfo0 = data.frame(locus = locus, moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  modOrder = order(-abs(cor(MEs, vac, use = "p")));
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  geneInfo = geneInfo0
  write.csv(geneInfo, file = "WGCNA_Corrs.csv", row.names = FALSE)

