source("functions.R")

##### data preparation
process_data()

##### ANOVA analysis
anova_analysis()
# alternative way - using DESeq2
diff_expr_polysome()

##### differential expression analysis of mRNAs using DESeq2
diff_expr()

##### initialize gene sets
initialize_gene_sets()

##### VENN diagrams of gene sets based on RNA-Seq
venn(list("LE-vs-L" = rownames(le), "LEKU-vs-L" = rownames(leku),
	"LE-vs-LEKU" = rownames(le.leku)), T, "rna-seq")


# L-LE - changes in only monosome fractions
t <- pol.all[!is.na(L.LE.light.padj) & L.LE.light.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj > 0.5 &
                     !is.na(L.LE.monosome.padj) & L.LE.monosome.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le)]




t <- pol.all[!is.na(L.LE.all.padj) & L.LE.all.padj < 0.01 & !ensembl_gene_id %in% rownames(le)]

gene <- t$ensembl_gene_id[1]
cond1 <- "L"
cond2 <- "LE"


# cond1 and cond2 arguments are either "L", "LE", "LEKU"
# fraction argument is either "monosome", "light", "heavy" or "all"
d <- readRDS("files/data-table.rds")
d <- d[order(cond, pf)]

# TODO - can't return a data frame back to data.table - need to figure out how to do it differently
# also need to reduce a number of genes - the whole data set is too large
d[,list(sig.pf = posthoc_test_pf(data.frame(value = value, cond = cond, pf = pf), "L", "LE")), by = "ensembl_gene_id"]

posthoc_test_pf <- function(d, cond1, cond2) {
        d <- d[d$cond == cond1 | d$cond == cond2,]
        res <- pairwise.t.test(d$value, paste(d$cond, d$pf, sep = "_"), p.adjust.method="BH")
        res <- res$p.value
        p.vals <- NULL
        i <- 20
        j <- 8
        for(z in 1:6){
                p.vals[z] <- res[i, j]
                i <- i + 1
                j <- j + 1
        }
        i <- 13
        j <- 1
        for(z in 7:13){
                p.vals[z] <- res[i, j]
                i <- i + 1
                j <- j + 1
        }
        return(data.frame(pf = 4:16, p.adj = p.vals))
        # return(1)
}

plot_genes("ENSMUSG00000098243", "test")


# > dim(t)
# [1] 25 18
plot_genes(t[,ensembl_gene_id], "L-LE-monosome-with-de")
t <- pol.all[!is.na(L.LE.light.padj) & L.LE.light.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj > 0.5 &
                     !is.na(L.LE.monosome.padj) & L.LE.monosome.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 30 18
plot_genes(t[,ensembl_gene_id], "L-LE-monosome-without-de")

# L-LE - changes in only light fractions
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj > 0.5 &
                     !is.na(L.LE.light.padj) & L.LE.light.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 24 18
plot_genes(t[,ensembl_gene_id], "L-LE-light-with-de")
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj > 0.5 &
                     !is.na(L.LE.light.padj) & L.LE.light.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 26 18
plot_genes(t[,ensembl_gene_id], "L-LE-light-without-de")

# L-LE - changes in only heavy fractions
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.light.padj) & L.LE.light.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 113 18
plot_genes(head(t[order(L.LE.heavy.padj), ensembl_gene_id], 30), "L-LE-heavy-with-de")
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.light.padj) & L.LE.light.padj > 0.5 &
                     !is.na(L.LE.heavy.padj) & L.LE.heavy.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 130 18
plot_genes(head(t[order(L.LE.heavy.padj), ensembl_gene_id], 30), "L-LE-heavy-without-de")

# L-LE - changes in only merge fractions
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.merge.padj) & L.LE.merge.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 717 18
plot_genes(head(t[order(L.LE.merge.padj), ensembl_gene_id], 30), "L-LE-merge-with-de")
t <- pol.all[!is.na(L.LE.monosome.padj) & L.LE.monosome.padj > 0.5 &
                     !is.na(L.LE.merge.padj) & L.LE.merge.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le)]
# > dim(t)
# [1] 592 18
plot_genes(head(t[order(L.LE.merge.padj), ensembl_gene_id], 30), "L-LE-merge-without-de")





# LE-LEKU - changes in only monosome fractions
t <- pol.all[!is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj > 0.5 &
                     !is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 6 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-monosome-with-de")
t <- pol.all[!is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj > 0.5 &
                     !is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 9 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-monosome-without-de")

# LE-LEKU - changes in only light fractions
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj > 0.5 &
                     !is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 15 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-light-with-de")
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj > 0.5 &
                     !is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 37 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-light-without-de")

# LE-LEKU - changes in only heavy fractions
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 51 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-heavy-with-de")
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.light.padj) & LE.LEKU.light.padj > 0.5 &
                     !is.na(LE.LEKU.heavy.padj) & LE.LEKU.heavy.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 78 18
plot_genes(t[,ensembl_gene_id], "LE-LEKU-heavy-without-de")

# LE-LEKU - changes in only merge fractions
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.merge.padj) & LE.LEKU.merge.padj < 0.01 &
                     ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 234 18
plot_genes(head(t[order(LE.LEKU.merge.padj), ensembl_gene_id], 30), "LE-LEKU-merge-with-de")
t <- pol.all[!is.na(LE.LEKU.monosome.padj) & LE.LEKU.monosome.padj > 0.5 &
                     !is.na(LE.LEKU.merge.padj) & LE.LEKU.merge.padj < 0.01 &
                     !ensembl_gene_id %in% rownames(le.leku)]
# > dim(t)
# [1] 220 18
plot_genes(head(t[order(LE.LEKU.merge.padj), ensembl_gene_id], 30), "LE-LEKU-merge-without-de")




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
