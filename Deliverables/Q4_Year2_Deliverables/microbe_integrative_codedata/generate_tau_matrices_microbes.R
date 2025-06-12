##this code generates tau values for the associations between our top 5% microbial families, and top 5% orthologs
##last updated LEF 6/12/25

##define a function to pull tau values
cor.test.tau <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="kendall", use="pairwise")[["estimate"]]
  z <- outer(
    colnames(x),
    colnames(x),
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

  
  ##use said function to generate tau values for each of our top 5% families one at a time
  
  Terasak_Corr = read.csv("Terasakiellaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Terasak=cor.test.tau(Terasak_Corr)
  warnings()
  write.table(t_Terasak, file ="Corr_results_tau_Terasakiellaceae.csv",sep=",")
  
  
  KilonCorr = read.csv("Kiloniellaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Kilon=cor.test.tau(KilonCorr)
  warnings()
  write.table(t_Kilon, file ="Corr_results_tau_Kiloniellaceae.csv",sep=",")
  
  
  EndoCorr = read.csv("Endozoicomonadaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Endo=cor.test.tau(EndoCorr)
  warnings()
  write.table(t_Endo, file ="Corr_results_tau_Endozoicomonadaceae.csv",sep=",")
  
  
  NisaeCorr = read.csv("Nisaeaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Nisae=cor.test.tau(NisaeCorr)
  warnings()
  write.table(t_Nisae, file ="Corr_results_tau_Nisaeaceae.csv",sep=",")
  
  
  ProlixCorr = read.csv("Prolixibacteraceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Prolix=cor.test.tau(ProlixCorr)
  warnings()
  write.table(t_Prolix, file ="Corr_results_tau_Prolixibacteraceae.csv",sep=",")
  
  
  DesulfoCorr = read.csv("Desulfobacteraceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Desulfo=cor.test.tau(DesulfoCorr)
  warnings()
  write.table(t_Desulfo, file ="Corr_results_tau_Desulfobacteraceae.csv",sep=",")
  
  
  RhodoCorr = read.csv("Rhodocyclaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Rhodo=cor.test.tau(RhodoCorr)
  warnings()
  write.table(t_Rhodo, file ="Corr_results_tau_Rhodocyclaceae.csv",sep=",")
  
  
  BlastoCorr = read.csv("Blastocatellaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Blasto=cor.test.tau(BlastoCorr)
  warnings()
  write.table(t_Blasto, file ="Corr_results_tau_Blastocatellaceae.csv",sep=",")
  
  
  HaloCorr = read.csv("Halobacteroidaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Halo=cor.test.tau(HaloCorr)
  warnings()
  write.table(t_Halo, file ="Corr_results_tau_Halobacteroidaceae.csv",sep=",")
  
  
  DesulfCorr = read.csv("Desulfotomaculales_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Desulf=cor.test.tau(DesulfCorr)
  warnings()
  write.table(t_Desulf, file ="Corr_results_tau_Desulfotomaculales.csv",sep=",")
  
  
  LentiCorr = read.csv("Lentisphaeraceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Lenti=cor.test.tau(LentiCorr)
  warnings()
  write.table(t_Lenti, file ="Corr_results_tau_Lentisphaeraceae.csv",sep=",")
  
  
  DethCorr = read.csv("Dethiobacteraceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Deth=cor.test.tau(DethCorr)
  warnings()
  write.table(t_Deth, file ="Corr_results_tau_Dethiobacteraceae.csv",sep=",")
  
  
  CamiCorr = read.csv("Caminicellaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Cami=cor.test.tau(CamiCorr)
  warnings()
  write.table(t_Cami, file ="Corr_results_tau_Caminicellaceae.csv",sep=",")
  
  
  ChloroCorr = read.csv("Chloroflexaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Chloro=cor.test.tau(ChloroCorr)
  warnings()
  write.table(t_Chloro, file ="Corr_results_tau_Chloroflexaceae.csv",sep=",")

  
  ChristCorr = read.csv("Christensenellaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Christ=cor.test.tau(ChristCorr)
  warnings()
  write.table(t_Christ, file ="Corr_results_tau_Christensenellaceae.csv",sep=",")
  
  
  ChromCorr = read.csv("Chromatiaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Chrom=cor.test.tau(ChromCorr)
  warnings()
  write.table(t_Chrom, file ="Corr_results_tau_Chromatiaceae.csv",sep=",")

  
  InquiCorr = read.csv("Inquilinaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Inqui=cor.test.tau(InquiCorr)
  warnings()
  write.table(t_Inqui, file ="Corr_results_tau_Inquilinaceae.csv",sep=",")  

  LentiCorr = read.csv("Lentimicrobiaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Lenti=cor.test.tau(LentiCorr)
  warnings()
  write.table(t_Lenti, file ="Corr_results_tau_Lentimicrobiaceae.csv",sep=",")  
  
  SpiroCorr = read.csv("Spirochaetaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Spiro=cor.test.tau(SpiroCorr)
  warnings()
  write.table(t_Spiro, file ="Corr_results_tau_Spirochaetaceae.csv",sep=",")
  
  
  OligoCorr = read.csv("Oligoflexaceae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Oligo=cor.test.tau(OligoCorr)
  warnings()
  write.table(t_Oligo, file ="Corr_results_tau_Oligoflexaceae.csv",sep=",")
  