source("functions.R")

##### data preparation
process_data()

##### ANOVA analysis
anova_analysis()
# alternative way - using DESeq2
diff_expr_two_time_courses_cond("L", "LE", "all")
diff_expr_two_time_courses_cond("LE", "LEKU", "all")

diff_expr_two_time_courses_cond("L", "LE", "monosome")
diff_expr_two_time_courses_cond("LE", "LEKU", "monosome")
diff_expr_two_time_courses_cond("L", "LEKU", "monosome")

diff_expr_two_time_courses_cond("L", "LE", "light")
diff_expr_two_time_courses_cond("LE", "LEKU", "light")
diff_expr_two_time_courses_cond("L", "LEKU", "light")

diff_expr_two_time_courses_cond("L", "LE", "heavy")
diff_expr_two_time_courses_cond("LE", "LEKU", "heavy")
diff_expr_two_time_courses_cond("L", "LEKU", "heavy")

##### differential expression analysis of mRNAs
diff_expr()

##### initialize gene sets
initialize_gene_sets()

##### VENN diagrams of gene sets based on RNA-Seq
venn(list("LE-vs-L" = rownames(le), "LEKU-vs-L" = rownames(leku),
	"LE-vs-LEKU" = rownames(le.leku)), T, "rna-seq")

# venn(list(LE = rownames(le), LEKU = rownames(leku)), T, "LE-LEKU")

# venn_ellipses(list(LE.up = rownames(le.up),
# 		  LE.down = rownames(le.down),
# 		  LEKU.up = rownames(leku.up),
# 		  LEKU.down = rownames(leku.down)), F, "LE-LEKU-up-down")

##### VENN diagrams of gene sets based ANOVA analysis
venn(list("LE-vs-L" = res.l.le$gene_id,
	"LEKU-vs-L" = res.l.leku$gene_id,
	"LE-vs-LEKU" = res.le.leku$gene_id), T, "polysome")

venn(list("RNA-Seq" = rownames(le), "Polysome" = res.l.le$gene_id),
	T, "correlations-LE-L")
venn(list("RNA-Seq" = rownames(leku), "Polysome" = res.l.leku$gene_id),
	T, "correlations-LEKU-L")
venn(list("RNA-Seq" = rownames(le.leku), "Polysome" = res.le.leku$gene_id),
	T, "correlations-LE-LEKU")

go_correlations()

# expression of these genes hasn't changed in le condition
t <- res.l.le.deseq[!(res.l.le.deseq$ensembl_gene_id %in% rownames(le)), ]
t <- t[order(t$gene_length),]
plot_genes(head(t$ensembl_gene_id, 40), "test")

# expression of these genes hasn't changed in le condition
t <- res.le.leku.deseq[!(res.le.leku.deseq$ensembl_gene_id %in% rownames(le.leku)), ]
plot_genes(head(t[order(t$padj),]$ensembl_gene_id, 40), "test")





p <- ggplot(t, aes(gene_length)) +
	geom_histogram() +
	scale_x_log10()
p

s <- t[t$gene_length <= 1e4, 1]
m <- t[t$gene_length > 1e4 & t$gene_length < 1e5, 1]
l <- t[t$gene_length >= 1e5, 1]

plot_genes(sample(s, 10), "gene-len-short-L-LE")
plot_genes(sample(l, 10), "gene-len-long-L-LE")






res <- data.frame(ensembl_gene_id = rownames(count.matrix))
res <- merge(res, gene.len[ , c(1, 4)])











cluster_correlations()

# get the clusters
clusts.all <- get_clust()



# GO analysis of mutational gene sets and their clusters
go_clust(clusts.all)
# plot GO plots for clusters
plot_go_clust()
# plot GO trees for clusters
quickgo_clust()


##### L-LE comparison
# > length(setA.filt[setA.filt %in% res.l.le$gene_id])
# [1] 412
# > length(setB.filt[setB.filt %in% res.l.le$gene_id])
# [1] 256
# > length(setC.filt[setC.filt %in% res.l.le$gene_id])
# [1] 134
# > length(setD.filt[setD.filt %in% res.l.le$gene_id])
# [1] 78
# > length(setE.filt[setE.filt %in% res.l.le$gene_id])
# [1] 1360
# > length(setF.filt[setF.filt %in% res.l.le$gene_id])
# [1] 989
# > length(setG.filt[setG.filt %in% res.l.le$gene_id])
# [1] 18
# > length(setH.filt[setH.filt %in% res.l.le$gene_id])
# [1] 2

##### LE-LEKU comparison
# > length(setA.filt[setA.filt %in% res.le.leku$gene_id])
# [1] 249
# > length(setB.filt[setB.filt %in% res.le.leku$gene_id])
# [1] 72
# > length(setC.filt[setC.filt %in% res.le.leku$gene_id])
# [1] 12
# > length(setD.filt[setD.filt %in% res.le.leku$gene_id])
# [1] 217
# > length(setE.filt[setE.filt %in% res.le.leku$gene_id])
# [1] 114
# > length(setF.filt[setF.filt %in% res.le.leku$gene_id])
# [1] 228
# > length(setG.filt[setG.filt %in% res.le.leku$gene_id])
# [1] 33
# > length(setH.filt[setH.filt %in% res.le.leku$gene_id])
# [1] 2

##### L-LEKU comparison
# > length(setA.filt[setA.filt %in% res.l.leku$gene_id])
# [1] 37
# > length(setB.filt[setB.filt %in% res.l.leku$gene_id])
# [1] 144
# > length(setC.filt[setC.filt %in% res.l.leku$gene_id])
# [1] 74
# > length(setD.filt[setD.filt %in% res.l.leku$gene_id])
# [1] 257
# > length(setE.filt[setE.filt %in% res.l.leku$gene_id])
# [1] 438
# > length(setF.filt[setF.filt %in% res.l.leku$gene_id])
# [1] 1161
# > length(setG.filt[setG.filt %in% res.l.leku$gene_id])
# [1] 16
# > length(setH.filt[setH.filt %in% res.l.leku$gene_id])
# [1] 0

plot_groups()
go_groups()
cluster_groups()

get_clust_genes("clusts-A-L-LE", 3)
get_clust_genes("clusts-A-L-LEKU", 5)
get_clust_genes("clusts-A-LE-LEKU", 6)
get_clust_genes("clusts-B-L-LE", 5)
get_clust_genes("clusts-B-L-LEKU", 6)
get_clust_genes("clusts-B-LE-LEKU", 4)
get_clust_genes("clusts-C-L-LE", 4)
get_clust_genes("clusts-C-L-LEKU", 6)
get_clust_genes("clusts-C-LE-LEKU", 6)
get_clust_genes("clusts-D-L-LE", 5)
get_clust_genes("clusts-D-L-LEKU", 5)
get_clust_genes("clusts-D-LE-LEKU", 5)
get_clust_genes("clusts-E-L-LE", 4)
get_clust_genes("clusts-E-L-LEKU", 5)
get_clust_genes("clusts-E-LE-LEKU", 4)
get_clust_genes("clusts-F-L-LE", 4)
get_clust_genes("clusts-F-L-LEKU", 6)
get_clust_genes("clusts-F-LE-LEKU", 5)
get_clust_genes("clusts-G-L-LE", 3)
get_clust_genes("clusts-G-L-LEKU", 3)
get_clust_genes("clusts-G-LE-LEKU", 3)
get_clust_genes("clusts-H-L-LE", 4)
get_clust_genes("clusts-H-LE-LEKU", 4)
