##this is an R script to parse  WGCNA output and generate files for GOMWU and the other investigations##
##last updated by LEF 6/12/25
##first we load our package##
library(tidyr)
library(tibble)
library(janitor)

##initial formatting##
  ##input WGCNA output##
  data = read.csv("WGCNA_Corrs.csv")
  names(data)
  ##pull out just the rows we need (meta data and module membership values####
  data = data[c(1:2,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33)]
  ##import all your genes##
  genenames = read.csv("normalized_reads_corr.csv", check.names = FALSE)
  genenames=t(genenames)
  genenames<- genenames %>%
    row_to_names(row_number = 1)
  genenames=as.data.frame(genenames)
  genenames <- tibble::rownames_to_column(genenames, "locus")
  genenames=genenames[c(1:2)]

##now we will iterate through our modules and create two outputs for each: a list of all transcripts and GOMWU input##
##I will fully annotate the first chunk of code, the rest are repetitive##
  
  ##brown module##
  ##pull out just the rows we want##
  brown = data[data$moduleColor %in% c("brown"), ] 
  names(brown)
  brown = brown[c(1,13)] ##check names to select correct column##
  ##write this out to a csv##
  write.csv(brown, "brown_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(brown, genenames, by="locus" all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "brown_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##cyan module##
  ##pull out just the rows we want##
  cyan = data[data$moduleColor %in% c("cyan"), ] 
  names(cyan)
  cyan = cyan[c(1,3)]
  ##write this out to a csv##
  write.csv(cyan, "cyan_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(cyan, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "cyan_kme_GOMWU.csv", row.names = FALSE)

  
  ##darkgreen module##
  ##pull out just the rows we want##
  darkgreen = data[data$moduleColor %in% c("darkgreen"), ] 
  names(darkgreen)
  darkgreen = darkgreen[c(1,11)]
  ##write this out to a csv##
  write.csv(darkgreen, "darkgreen_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(darkgreen, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "darkgreen_kme_GOMWU.csv", row.names = FALSE)

  ##darkgrey module##
  ##pull out just the rows we want##
  darkgrey = data[data$moduleColor %in% c("darkgrey"), ] 
  names(darkgrey)
  darkgrey = darkgrey[c(1,16)]
  ##write this out to a csv##
  write.csv(darkgrey, "darkgrey_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(darkgrey, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "darkgrey_kme_GOMWU.csv", row.names = FALSE)

  ##ok start with the darkturquoise module##
  ##pull out just the rows we want##
  darkturquoise = data[data$moduleColor %in% c("darkturquoise"), ] 
  names(darkturquoise)
  darkturquoise = darkturquoise[c(1,17)]
  ##write this out to a csv##
  write.csv(darkturquoise, "darkturquoise_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(darkturquoise, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "darkturquoise_kme_GOMWU.csv", row.names = FALSE)  

  ##green module##
  ##pull out just the rows we want##
  green = data[data$moduleColor %in% c("green"), ] 
  names(green)
  green = green[c(1,9)]
  ##write this out to a csv##
  write.csv(green, "green_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(green, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "green_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##greenyellow module##
  ##pull out just the rows we want##
  greenyellow = data[data$moduleColor %in% c("greenyellow"), ] 
  names(greenyellow)
  greenyellow = greenyellow[c(1,7)]
  ##write this out to a csv##
  write.csv(greenyellow, "greenyellow_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(greenyellow, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "greenyellow_kme_GOMWU.csv", row.names = FALSE)

  
  ##lightcyan module##
  ##pull out just the rows we want##
  lightcyan = data[data$moduleColor %in% c("lightcyan"), ] 
  names(lightcyan)
  lightcyan = lightcyan[c(1,5)]
  ##write this out to a csv##
  write.csv(lightcyan, "lightcyan_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(lightcyan, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "lightcyan_kme_GOMWU.csv", row.names = FALSE)  
  
  
  ##lightyellow module##
  ##pull out just the rows we want##
  lightyellow = data[data$moduleColor %in% c("lightyellow"), ] 
  names(lightyellow)
  lightyellow = lightyellow[c(1,10)]
  ##write this out to a csv##
  write.csv(lightyellow, "lightyellow_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(lightyellow, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "lightyellow_kme_GOMWU.csv", row.names = FALSE)

  
  ##magenta module##
  ##pull out just the rows we want##
  magenta = data[data$moduleColor %in% c("magenta"), ] 
  names(magenta)
  magenta = magenta[c(1,8)]
  ##write this out to a csv##
  write.csv(magenta, "magenta_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(magenta, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "magenta_kme_GOMWU.csv", row.names = FALSE)
  
  ##midnightblue module##
  ##pull out just the rows we want##
  midnightblue = data[data$moduleColor %in% c("midnightblue"), ] 
  names(midnightblue)
  midnightblue = midnightblue[c(1,6)]
  ##write this out to a csv##
  write.csv(midnightblue, "midnightblue_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(midnightblue, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "midnightblue_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##orange module##
  ##pull out just the rows we want##
  orange = data[data$moduleColor %in% c("orange"), ] 
  names(orange)
  orange = orange[c(1,4)]
  ##write this out to a csv##
  write.csv(orange, "orange_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(orange, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "orange_kme_GOMWU.csv", row.names = FALSE)

  
  ##red module##
  ##pull out just the rows we want##
  red = data[data$moduleColor %in% c("red"), ] 
  names(red)
  red = red[c(1,14)]
  ##write this out to a csv##
  write.csv(red, "red_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(red, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "red_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##skyblue module##
  ##pull out just the rows we want##
  skyblue = data[data$moduleColor %in% c("skyblue"), ] 
  names(skyblue)
  skyblue = skyblue[c(1,15)]
  ##write this out to a csv##
  write.csv(skyblue, "skyblue_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(skyblue, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "skyblue_kme_GOMWU.csv", row.names = FALSE)

  
  ##white module##
  ##pull out just the rows we want##
  white = data[data$moduleColor %in% c("white"), ] 
  names(white)
  white = white[c(1,12)]
  ##write this out to a csv##
  write.csv(white, "white_transcripts.csv", row.names =FALSE)
  
  ##create the GOMWU input##
  Merged = merge(white, genenames, by="locus", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "white_kme_GOMWU.csv", row.names = FALSE)
  