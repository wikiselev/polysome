source("functions.R")

##### data preparation
process_data()

# alternative way - using DESeq2
diff_expr_polysome()

##### differential expression analysis of mRNAs using DESeq2
diff_expr_rna()

##### posthoc analysis of polysome fraction distributions
posthoc_analysis()

##### initialize gene sets
initialize_gene_sets()

##### for the rest of analysis - see report.Rmd and its output - report.pdf
for_report()

# # clustering example - not using it at the moment
# cluster_correlations()
# # get the clusters
# clusts.all <- get_clust()
# # GO analysis of mutational gene sets and their clusters
# go_clust(clusts.all)
