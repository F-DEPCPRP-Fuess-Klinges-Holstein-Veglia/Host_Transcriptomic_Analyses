##this code generates tau values for the associations between our top 5% viral orders, and top 5% orthologs
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

  Bing_Corr = read.csv("Bingvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Bing=cor.test.tau(Bing_Corr)
  warnings()
  write.table(t_Bing, file ="Corr_results_tau_Bingvirus.csv",sep=",")
  
  BL3_Corr = read.csv("BrochothrixphageBL3_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_BL3=cor.test.tau(BL3_Corr)
  warnings()
  write.table(t_BL3, file ="Corr_results_tau_BrochothrixphageBL3.csv",sep=",")
  
  Burro_Corr = read.csv("Burrovirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Burro=cor.test.tau(Burro_Corr)
  warnings()
  write.table(t_Burro, file ="Corr_results_tau_Burrovirus.csv",sep=",")
  
  Carm_Corr = read.csv("Carmenvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Carm=cor.test.tau(Carm_Corr)
  warnings()
  write.table(t_Carm, file ="Corr_results_tau_Carmenvirus.csv",sep=",")
  
  CT19406C_Corr = read.csv("ClostridiumphagephiCT19406C_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_CT19406C=cor.test.tau(CT19406C_Corr)
  warnings()
  write.table(t_CT19406C, file ="Corr_results_tau_ClostridiumphagephiCT19406C.csv",sep=",")
  
  Decur_Corr = read.csv("Decurrovirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Decur=cor.test.tau(Decur_Corr)
  warnings()
  write.table(t_Decur, file ="Corr_results_tau_Decurrovirus.csv",sep=",")
  
  Godo_Corr = read.csv("Godonkavirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Godo=cor.test.tau(Godo_Corr)
  warnings()
  write.table(t_Godo, file ="Corr_results_tau_Godonkavirus.csv",sep=",")
  
  Halcy_Corr = read.csv("Halcyonevirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Halcy=cor.test.tau(Halcy_Corr)
  warnings()
  write.table(t_Halcy, file ="Corr_results_tau_Halcyonevirus.csv",sep=",")
  
  KO2_Corr = read.csv("KlebsiellaphagephiKO2_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_KO2=cor.test.tau(KO2_Corr)
  warnings()
  write.table(t_KO2, file ="Corr_results_tau_KlebsiellaphagephiKO2.csv",sep=",")
  
  Lacu_Corr = read.csv("Lacusarxvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Lacu=cor.test.tau(Lacu_Corr)
  warnings()
  write.table(t_Lacu, file ="Corr_results_tau_Lacusarxvirus.csv",sep=",")
  
  Light_Corr = read.csv("Lightbulbvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Light=cor.test.tau(Light_Corr)
  warnings()
  write.table(t_Light, file ="Corr_results_tau_Lightbulbvirus.csv",sep=",")
  
  Nona_Corr = read.csv("Nonagvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Nona=cor.test.tau(Nona_Corr)
  warnings()
  write.table(t_Nona, file ="Corr_results_tau_Nonagvirus.csv",sep=",")
  
  
  PseudoCorr = read.csv("PseudoalteromonasphageH103_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Pseudo=cor.test.tau(PseudoCorr)
  warnings()
  write.table(t_Pseudo, file ="Corr_results_tau_PseudoalteromonasphageH103.csv",sep=",")
  
  
  JBD44Corr = read.csv("PseudomonasphageJBD44_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_JBD44=cor.test.tau(JBD44Corr)
  warnings()
  write.table(t_JBD44, file ="Corr_results_tau_PseudomonasphageJBD44.csv",sep=",")
  
  
  SchitCorr = read.csv("Schitoviridae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Schit=cor.test.tau(SchitCorr)
  warnings()
  write.table(t_Schit, file ="Corr_results_tau_Schitoviridae.csv",sep=",")
  
  
  SpizCorr = read.csv("Spizizenvirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Spiz=cor.test.tau(SpizCorr)
  warnings()
  write.table(t_Spiz, file ="Corr_results_tau_Spizizenvirus.csv",sep=",")
  
  
  StreptoCorr = read.csv("StreptococcusphagephiARI0746_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Strepto=cor.test.tau(StreptoCorr)
  warnings()
  write.table(t_Strepto, file ="Corr_results_tau_StreptococcusphagephiARI0746.csv",sep=",")
  
  
  NJ2Corr = read.csv("StreptococcusphagephiNJ2_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_NJ2=cor.test.tau(NJ2Corr)
  warnings()
  write.table(t_NJ2, file ="Corr_results_tau_StreptococcusphagephiNJ2.csv",sep=",")
  
  
  UetakCorr = read.csv("Uetakevirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Uetak=cor.test.tau(UetakCorr)
  warnings()
  write.table(t_Uetak, file ="Corr_results_tau_Uetakevirus.csv",sep=",")
  
  
  VertoCorr = read.csv("Vertoviridae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Verto=cor.test.tau(VertoCorr)
  warnings()
  write.table(t_Verto, file ="Corr_results_tau_Vertoviridae.csv",sep=",")
  
  Vivid_Corr = read.csv("Vividuovirus_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Vivid=cor.test.tau(Vivid_Corr)
  warnings()
  write.table(t_Vivid, file ="Corr_results_tau_Vividuovirus.csv",sep=",")
  
  
  ZobeCorr = read.csv("Zobellviridae_Tau_Matrix.csv", check.names = FALSE)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  t_Zobe=cor.test.tau(ZobeCorr)
  warnings()
  write.table(t_Zobe, file ="Corr_results_tau_Zobellviridae.csv",sep=",")
  