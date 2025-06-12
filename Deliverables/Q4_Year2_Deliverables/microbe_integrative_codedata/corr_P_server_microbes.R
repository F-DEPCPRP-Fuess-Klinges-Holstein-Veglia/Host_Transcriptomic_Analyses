##this code runs our microbiome correlations generating a pvalue for each comparison
##it is best run on a super computer to enhance computational power and speed
##last updated LEF 6/12/25

##define correlation function##
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="kendall", use="pairwise")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

  ##input data##
  MB_Corr = read.csv("final_correlation_matrix_families.csv", row.names=1)
  ##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
  p=cor.test.p(MB_Corr)
  warnings()
  write.table(p, file ="Corr_results_pval.csv",sep=",")