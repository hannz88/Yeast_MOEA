# R 3.6.1
# to use PANTHER Overrepresentation test

# load dependencies----
if(!require(httr)){
  install.packages("httr")
  library(httr)
}
if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
}

# load the RData----
load("file path for RData from NSGA2_Postprocessing.R")


# use a loop to go through all the solutions to check with PANTHER API----
library(httr)
bin_list <- pareto_front_genome  # use this if jaccard coefficient is not used
bg_genes <- rownames(data)
ref_gene_list <- paste(bg_genes, collapse = ",")  # reference genes/ background genes
k <- length(bg_genes)
sig_gen <- c()
for (i in 1:length(bin_list)) {
  #print(i)
  gene_vect <- bin_list[[i]][1:k]  # get the list of genome
  test_candidate_list <- bg_genes[c(gene_vect==1)]  # candidate genes for testing
  new_test_candidate <- paste(test_candidate_list, collapse=",")
  response <- POST("http://pantherdb.org/services/oai/pantherdb/enrich/overrep", query=list(geneInputList=new_test_candidate, 
                                                                                            organism="559292", refInputList=ref_gene_list, 
                                                                                            annotDataSet="GO:0003674", refOrganism="559292",
                                                                                            enrichmentTestType="FISHER", 
                                                                                            correction="FDR"))
  resp <- httr::content(response, "parsed")  # the result is ordered by ascending pvalue
  for (j in 1:length(resp$results$result)) {  # return the significant bicluster number
    if (resp$results$result[[j]]$fdr < 0.05) {
      #print(resp$results$result[[j]]$fdr)
      sig_gen <- c(sig_gen, i)
      break
    }
  }
}
sig_gen

# then the results of the significant biclusters could be checked by either going to to
# PANTHER site or here using

# to check on the site, you'll need the background list and the candidate list
# double check the results with pantherdb
# write the candidate genes and background genes into separate files----
gene_vect <- bin_list[[i]][1:k]  # replace i with the ith biclusters that is significant
candidate_list <- bg_genes[c(gene_vect==1)]
lapply(bg_bg_genes, write, "file path/background_genes.txt", append=TRUE)
lapply(candidate_list, write, "file path/candidate_genes.txt", append=TRUE)


# to check the results of a significant bicluster programatically----
bg_genes <- rownames(data)
k <- length(bg_genes)
bin_list <- pool$population  # uncomment this if jaccard coefficient is not used
gene_vect <- bin_list[[73]][1:k]  # get the list of genome
test_candidate_list <- bg_genes[c(gene_vect==1)]  # candidate genes for testing
new_test_candidate <- paste(test_candidate_list, collapse=",")
ref_gene_list <- paste(bg_genes, collapse = ",")  # reference genes/ background genes
# n2 <- dput(as.character(test_candidate_list))
response <- POST("http://pantherdb.org/services/oai/pantherdb/enrich/overrep", query=list(geneInputList=new_test_candidate, 
                                                                                          organism="559292", refInputList=ref_gene_list, 
                                                                                          annotDataSet="GO:0003674", refOrganism="559292",
                                                                                          enrichmentTestType="FISHER", 
                                                                                          correction="None"))
resp <- httr::content(response, "parsed")  # the result is ordered by ascending pvalue
resp$results
# to draw the heatmap of a bicluster----
target_genome <- bin_list[[i]]  # change the i to the ith bicluster you want to draw
row_vect <- target_genome[1:nrow(data)]
col_vect <- target_genome[-(1:nrow(data))]
bic <- data[c(row_vect==1), c(col_vect==1)]
heatmap.2(as.matrix(bic), col = redgreen(75), scale="row", key=T, keysize=1.5, density.info="none",trace="none",cexCol=0.9)