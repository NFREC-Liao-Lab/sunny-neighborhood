library(readxl)

data_raw <- data.frame(read_xlsx("genes assign to taxa.xlsx", skip=2))

index <- grep("taxa", names(data_raw))
data2 <- data_raw[,(index-1):(index+8)]
data2 <- data_raw[,(index+1):(index+8)]

nrow <- NROW(data2)
result <- matrix(nrow=nrow-1, ncol=8, dimnames=list(data_raw[-nrow, index], colnames(data2)))
fisher_list <- list()
groups <- list(1:2, 3:4, 5:6, 7:8)
for(g in groups){
  data_sub <- data2[,g]

  m <- data_sub[nrow, 1]
  n <- data_sub[nrow, 2]
  data_sub <- data_sub[-nrow, ]
  result[,g] <- t(apply(data_sub, 1, function(x){
      total <- sum(x)
      p1 <- phyper(x[1], m, n, total, lower.tail=F)
      p2 <- phyper(x[2], n, m, total, lower.tail=F)
      return(c(p1, p2))
  }))

  fisher_list[[g[1]]] <- apply(data_sub, 1, function(x){
    ctable <- matrix(c(x, m-x[1], n-x[2]), nrow=2)
    fisher.test(ctable)$p.value
  })

}

cname <- c('A. FunGene_P.tri (9988)_P.tri (P.tri x P.tri)_gene counts_vs_P. tri (Ptri x P. taeda)_gene_counts',
'B. FunGene_P.tri (4771)_P. taeda (P.tri x P.taeda)_gene counts_vs_P. taeda (P.taedax P.taeda)_gene counts',
'C. FunGene_P.taeda (1527)_P.tri (P.tri x P.tri)_gene counts_vs_P. tri (Ptri x P. taeda)_gene_counts',
'D. FunGene_P.taeda (2312)_P. taeda (P.tri x P.taeda)_gene counts_vs_P. taeda (P.taedax P.taeda)_gene counts')

result_fisher <- do.call("cbind", fisher_list)
dimnames(result_fisher) <- list(rownames(result), cname)

write.csv(result_fisher, "genes_to_taxa_fisher.csv")
write.csv(result, "genes_to_taxa_hyper.csv")
