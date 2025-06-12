##This is a general WGCNA Script to construct a single network for all our samples across all species##
##it uses normalized ortholog read counts to generate this network
##last updated LEF 6/12/25

library(WGCNA)
library(tidyverse)
options(stringsAsFactors = FALSE);
##Upload the data from samples this data has 0 filtering before normalization, only genes removed are those with 0 variance (see matrix generation code)##
Data = read.csv("normalized_reads_corr.csv", check.names = FALSE); ##change this to match the group of interest##
##Check the data##
dim(Data);
names(Data);

##format the data##
datExpr0 = as.data.frame(Data);
datExpr0<-datExpr0 %>% remove_rownames %>% column_to_rownames(var="sample_R")
##remove 0 variance rows##
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

##print the genes that were removed##

if (!gsg$allOK)
{
  #Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
##there are none##

##generate final dimensions##
dim(datExpr0)

collectGarbage();

##let's start building a network##
##choose soft threshold##
spt <- pickSoftThreshold(datExpr0) 
spt

softPower = 6; ##based on spt##
adjacency = adjacency(datExpr0,
                      selectCols = NULL,
                      type = "signed",
                      power = softPower,
                      corFnc = "bicor")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType="signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
##unremark this if you want to see initial dendro
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#labels = FALSE, hang = 0.04);


# Set min module size pretty small since were only working with ~7500 orthos##
minModuleSize = 15;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Unremark to Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#dendroLabels = FALSE, hang = 0.03,
#addGuide = TRUE, guideHang = 0.05,
#main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# unremark to Plot the result
#sizeGrWindow(7, 6)
#plot(METree, main = "Clustering of module eigengenes",
#xlab = "", sub = "")
MEDissThres = 0.2
# Plot the cut line into the dendrogram
#abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors (a lot get merged)
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

table(mergedColors)


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent correlation analyses##
save(MEs, moduleLabels, datExpr0, moduleColors, geneTree, file = "WGCNANetworkConstruction.RData")
