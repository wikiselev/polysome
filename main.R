source("functions.R")

##### data preparation
process_data()

##### ANOVA analysis
anova_analysis()
# alternative way - using DESeq2
diff_expr_polysome()

##### differential expression analysis of mRNAs using DESeq2
diff_expr()

##### posthoc analysis of polysome fraction distributions
posthoc_analysis()

##### initialize gene sets
initialize_gene_sets()

##### VENN diagrams of gene sets based on RNA-Seq
venn(list("LE-vs-L" = rownames(le), "LEKU-vs-L" = rownames(leku),
	"LE-vs-LEKU" = rownames(le.leku)), T, "rna-seq")


t <- pol.all[!is.na(L.LE.all.padj) & L.LE.all.padj < 0.01 & !ensembl_gene_id %in% rownames(le)]
genes <- t[, ensembl_gene_id]
t1 <- posthoc.l.le[ensembl_gene_id %in% genes]
setkey(t1, "ensembl_gene_id")
min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
setkey(min.pvals, "ensembl_gene_id")
t1 <- t1[min.pvals]
t1 <- t1[sig.pf == min.pval]
GO(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), all.genes, "L-LE-monosome-change", 0.05)
GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "L-LE-light-change", 0.05)
GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "L-LE-heavy-change", 0.05)
t1 <- t1[order(min.pval)]
plot_genes(head(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), 20), "L-LE-monosome-change")
plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "L-LE-light-change")
plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "L-LE-heavy-change")
plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "L-LE-heavy16-change")

t <- pol.all[!is.na(L.LEKU.all.padj) & L.LEKU.all.padj < 0.01 & !ensembl_gene_id %in% rownames(leku)]
genes <- t[, ensembl_gene_id]
t1 <- posthoc.l.leku[ensembl_gene_id %in% genes]
setkey(t1, "ensembl_gene_id")
min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
setkey(min.pvals, "ensembl_gene_id")
t1 <- t1[min.pvals]
t1 <- t1[sig.pf == min.pval]
GO(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), all.genes, "L-LEKU-monosome-change", 0.05)
GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "L-LEKU-light-change", 0.05)
GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "L-LEKU-heavy-change", 0.05)
t1 <- t1[order(min.pval)]
plot_genes(head(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), 20), "L-LEKU-monosome-change")
plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "L-LEKU-light-change")
plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "L-LEKU-heavy-change")
plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "L-LEKU-heavy16-change")

t <- pol.all[!is.na(LE.LEKU.all.padj) & LE.LEKU.all.padj < 0.01 & !ensembl_gene_id %in% rownames(le.leku)]
genes <- t[, ensembl_gene_id]
t1 <- posthoc.le.leku[ensembl_gene_id %in% genes]
setkey(t1, "ensembl_gene_id")
min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
setkey(min.pvals, "ensembl_gene_id")
t1 <- t1[min.pvals]
t1 <- t1[sig.pf == min.pval]
GO(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), all.genes, "LE-LEKU-monosome-change", 0.05)
GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "LE-LEKU-light-change", 0.05)
GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "LE-LEKU-heavy-change", 0.05)
t1 <- t1[order(min.pval)]
plot_genes(head(unique(t1[pf < 8 & min.pval < 0.01, ensembl_gene_id]), 20), "LE-LEKU-monosome-change")
plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "LE-LEKU-light-change")
plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "LE-LEKU-heavy-change")
plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "LE-LEKU-heavy16-change")

# ##### VENN diagrams of gene sets based ANOVA analysis
# venn(list("LE-vs-L" = res.l.le$gene_id,
# 	"LEKU-vs-L" = res.l.leku$gene_id,
# 	"LE-vs-LEKU" = res.le.leku$gene_id), T, "polysome")
#
# venn(list("RNA-Seq" = rownames(le), "Polysome" = res.l.le$gene_id),
# 	T, "correlations-LE-L")
# venn(list("RNA-Seq" = rownames(leku), "Polysome" = res.l.leku$gene_id),
# 	T, "correlations-LEKU-L")
# venn(list("RNA-Seq" = rownames(le.leku), "Polysome" = res.le.leku$gene_id),
# 	T, "correlations-LE-LEKU")

go_correlations()
cluster_correlations()

# get the clusters
clusts.all <- get_clust()
# GO analysis of mutational gene sets and their clusters
go_clust(clusts.all)

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
