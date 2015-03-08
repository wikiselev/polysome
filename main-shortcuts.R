process_data <- function() {
	d <- read.table("original-data/p28_expression_genes_edited.txt", sep = "\t", header = T)

	# remove Poly_L_2_17 sample
	d <- d[,!grepl("LUL16_e14", colnames(d))]

	# rename other samples in LUL16
	for(i in 15:40) {
		colnames(d) <- gsub(paste0("LUL16_e", i), paste0("LUL16_e", i-1), colnames(d))
	}
        colnames(d)[4] <- "ensembl_gene_id"
	d1 <- d[,c(4,13:dim(d)[2])]
	saveRDS(d1, "files/data-matrix.rds")

	# prepare a table for clustering
	tmp <- d1[,2:157]
	tmp <- as.data.frame(t(scale(t(tmp), center = T, scale = T)))
	tmp <- cbind(d1$ensembl_gene_id, tmp)
	colnames(tmp)[1] <- "ensembl_gene_id"
	saveRDS(tmp, "files/data-matrix-norm.rds")

	# prepare gene name/length annotation table
	ann <- read.table("annotations/genes-GRCm38.p3.txt", sep = "\t", header = T)
	colnames(ann) <- c("ensembl_gene_id", "gene_start", "gene_end", "gene_name")
	ann$gene_length <- abs(ann$gene_end - ann$gene_start)
	ann <- ann[,c(1,4,5)]
	saveRDS(as.data.table(ann), "files/ann-table.rds")

	# prepare a data table for anova analysis
	d2 <- rename_samples_plot(d1)
        d2 <- as.data.table(merge(as.data.frame(d2), ann))
        # 193 genes do not have gene names, so they are excluded-
        # but they never appear in the any list of significant genes
        # so we can safely ignore them
	saveRDS(d2, "files/data-table.rds")
	d3 <- as.data.frame(t(scale(t(d1[,2:157]), center = T, scale = T)))
	d3 <- cbind(d1$ensembl_gene_id, d3)
	colnames(d3)[1] <- "ensembl_gene_id"
	d4 <- rename_samples_plot(d3)
	d4 <- as.data.table(merge(as.data.frame(d4), ann))
	saveRDS(d4, "files/data-table-norm.rds")

	# prepare a data table for plotting
	# average replicates - find mean and sd
	d5 <- d2[,list(mean_value=mean(value), sd_value=sd(value)),by=c("ensembl_gene_id", "gene_name", "cond", "pf")]
	saveRDS(d5, "files/plot-table.rds")
	d6 <- d4[,list(mean_value=mean(value), sd_value=sd(value)),by=c("ensembl_gene_id", "gene_name", "cond", "pf")]
	saveRDS(d6, "files/plot-table-norm.rds")

	# principal component analysis + correlations of polysome data
	pca()

	# correlations of RNA-Seq data
	files <- list.files("htseq-count/", full.names = FALSE)
	count.matrix <- count_matrix_from_files(files, T)
	colnames(count.matrix) <- unlist(lapply(strsplit(colnames(count.matrix), "_"), function(x){return(x[2])}))
	count.matrix <- count.matrix[,order(colnames(count.matrix))]
	cor_heatmap(count.matrix, "rna-seq")
}

diff_expr_polysome <- function() {
        # taking into account all fractions
        diff_expr_two_time_courses_cond("L", "LE", "all")
        diff_expr_two_time_courses_cond("LE", "LEKU", "all")
        diff_expr_two_time_courses_cond("L", "LEKU", "all")
        # taking into account only monosome [4:7] fractions
        diff_expr_two_time_courses_cond("L", "LE", "monosome")
        diff_expr_two_time_courses_cond("LE", "LEKU", "monosome")
        diff_expr_two_time_courses_cond("L", "LEKU", "monosome")
        # taking into account only light [8:10] fractions
        diff_expr_two_time_courses_cond("L", "LE", "light")
        diff_expr_two_time_courses_cond("LE", "LEKU", "light")
        diff_expr_two_time_courses_cond("L", "LEKU", "light")
        # taking into account only heavy [11:16] fractions
        diff_expr_two_time_courses_cond("L", "LE", "heavy")
        diff_expr_two_time_courses_cond("LE", "LEKU", "heavy")
        diff_expr_two_time_courses_cond("L", "LEKU", "heavy")
        # taking into account only light and heavy [8:16] fractions
        diff_expr_two_time_courses_cond("L", "LE", "merge")
        diff_expr_two_time_courses_cond("LE", "LEKU", "merge")
        diff_expr_two_time_courses_cond("L", "LEKU", "merge")
}

diff_expr_rna <- function() {
	# differential expression analysis of RNA-seq data
	diff_expr_pairwise("L", "LE")
	diff_expr_pairwise("L", "LEKU")
	diff_expr_pairwise("LE", "LEKU")
}

initialize_gene_sets <- function() {
	# define rna-seq gene sets
	le <<- get_diff_expr("L.LE", 0.01)
	leku <<- get_diff_expr("L.LEKU", 0.01)
	le.leku <<- get_diff_expr("LE.LEKU", 0.01)

# 	# get significant genes from ANOVA analysis of polysome fractions
# 	# output: global variables res.l.le, res.le.leku, res.l.leku
# 	# switched to DESeq2 package and don't use ANOVA analysis anymore
# 	# now this function is commented:
# 	get_sig_genes_anova()

	# initialize polysome gene sets
	pol.all <- data.frame(ind = c(1:19600))
        coln <- c()
	for(cond in c("L.LE", "LE.LEKU", "L.LEKU")) {
                d <- readRDS(paste0("files/diff-expr-", cond, "-polysome-cond-all.rds"))
                pol.all <- cbind(pol.all, d$padj)
                coln <- c(coln, paste0(cond, ".all.padj"))
	}
        pol.all <- pol.all[,2:dim(pol.all)[2]]
        colnames(pol.all) <- coln
        pol.all$ensembl_gene_id <- rownames(d)
        # remove genes where padj is NA in all conditions and fractions
        pol.all <- pol.all[apply(pol.all[,1:3], 1, function(x) !all(is.na(x))),]

        # file gene_lengths_GRCh37.p13.txt was downloaded from Ensembl Biomart on 06/07/14
        gene.len <- readRDS("files/ann-table.rds")
        # add gene lengths
        pol.all <- merge(pol.all, gene.len)
        # note that 110 genes do not have gene names and are excluded from pol.all
        pol.all <<- as.data.table(pol.all)

	# old polysome differential expression gene sets - now we first split
	# polysome data by fraction and then do differential expression analysis
	# in each fraction category separately - see above
	res.l.le.deseq <- get_diff_expr("L.LE-polysome-cond-all", 0.01)
	res.l.leku.deseq <- get_diff_expr("L.LEKU-polysome-cond-all", 0.01)
	res.le.leku.deseq <- get_diff_expr("LE.LEKU-polysome-cond-all", 0.01)
	# add gene lengths
	res.l.le.deseq$ensembl_gene_id <- rownames(res.l.le.deseq)
	res.l.le.deseq <<- merge(as.data.frame(res.l.le.deseq), gene.len)
	res.l.leku.deseq$ensembl_gene_id <- rownames(res.l.leku.deseq)
	res.l.leku.deseq <<- merge(as.data.frame(res.l.leku.deseq), gene.len)
	res.le.leku.deseq$ensembl_gene_id <- rownames(res.le.leku.deseq)
	res.le.leku.deseq <<- merge(as.data.frame(res.le.leku.deseq), gene.len)

	posthoc.l.le <<- readRDS("files/posthoc-pf-sig-L-LE.rds")
	posthoc.l.leku <<- readRDS("files/posthoc-pf-sig-L-LEKU.rds")
	posthoc.le.leku <<- readRDS("files/posthoc-pf-sig-LE-LEKU.rds")

# 	# this gene set was required for GO analysis
# 	all.genes <<- unique(c(res.l.le.deseq$ensembl_gene_id,
# 	                       res.l.leku.deseq$ensembl_gene_id,
# 	                       res.le.leku.deseq$ensembl_gene_id,
# 	                       rownames(le),
# 	                       rownames(leku)))
}

anova_analysis <- function() {
	t <- readRDS("files/data-table.rds")
	# perform a two way anova test with repeated measures
	# d3 <- t[gene_id %in% t[,gene_id][1:10]]
	# res <- aov(value ~ cond * pf + Error(gene_id/(cond*pf)), data = d3)
	# simple two-way anova
	anova.l.le <- anova_results(t[cond == "L" | cond == "LE"])
	saveRDS(anova.l.le, "files/anova-result-L-LE.rds")
	anova.le.leku <- anova_results(t[cond == "LE" | cond == "LEKU"])
	saveRDS(anova.le.leku, "files/anova-result-LE-LEKU.rds")
	anova.l.leku <- anova_results(t[cond == "L" | cond == "LEKU"])
	saveRDS(anova.l.leku, "files/anova-result-L-LEKU.rds")
}

merge_anova_cond_int <- function(res) {
	# merge anova results - see description in the old report.pdf file in
	# report folder
	res.merged <- res$int
	res.merged$padj.cond <- res$cond$padj
	res.merged.sig <- res.merged[(padj < 0.01 & padj.cond < 0.05) | padj.cond < 0.01]
	return(res.merged.sig)
}

get_sig_genes_anova <- function() {
	res.l.le <<- merge_anova_cond_int(readRDS("files/anova-result-L-LE.rds"))
	res.le.leku <<- merge_anova_cond_int(readRDS("files/anova-result-LE-LEKU.rds"))
	res.l.leku <<- merge_anova_cond_int(res <- readRDS("files/anova-result-L-LEKU.rds"))
}

get_sig_genes_deseq <- function() {
	plot_genes(head(rownames(res.l.le.deseq[order(res.l.le.deseq$padj),]), 10),
		"L-LE-condition-pf-deseq")
	plot_genes(head(rownames(res.l.leku.deseq[order(res.l.leku.deseq$padj),]), 10),
		"L-LEKU-condition-pf-deseq")
	plot_genes(head(rownames(res.le.leku.deseq[order(res.le.leku.deseq$padj),]), 10),
		"LE-LEKU-condition-pf-deseq")
}

posthoc_analysis <- function() {
        d <- readRDS("files/data-table.rds")
        d <- d[order(cond, pf)]

        # also need to reduce a number of genes - the whole data set is too large
        res <- d[,list(sig.pf = posthoc_test_pf(data.frame(value = value, cond = cond, pf = pf), "L", "LE")), by = "ensembl_gene_id"]
        res$pf <- rep(4:16, length(unique(d[,ensembl_gene_id])))
        saveRDS(res, "files/posthoc-pf-sig-L-LE.rds")

        res <- d[,list(sig.pf = posthoc_test_pf(data.frame(value = value, cond = cond, pf = pf), "L", "LEKU")), by = "ensembl_gene_id"]
        res$pf <- rep(4:16, length(unique(d[,ensembl_gene_id])))
        saveRDS(res, "files/posthoc-pf-sig-L-LEKU.rds")

        res <- d[,list(sig.pf = posthoc_test_pf(data.frame(value = value, cond = cond, pf = pf), "LE", "LEKU")), by = "ensembl_gene_id"]
        res$pf <- rep(4:16, length(unique(d[,ensembl_gene_id])))
        saveRDS(res, "files/posthoc-pf-sig-LE-LEKU.rds")
}

go_correlations <- function() {
	GO(rownames(le)[rownames(le) %in% res.l.le$gene_id], all.genes, "corr-L-LE", 0.05)
	GO(rownames(leku)[rownames(leku) %in% res.l.leku$gene_id], all.genes, "corr-L-LEKU", 0.05)
	GO(rownames(le.leku)[rownames(le.leku) %in% res.le.leku$gene_id], all.genes, "corr-LE-LEKU", 0.05)
}

cluster_correlations <- function() {
	# clustering
	d <- readRDS("files/data-matrix-norm.rds")
	rownames(d) <- d[,1]
	d <- as.matrix(d[,2:dim(d)[2]])

	clust_boot(d[rownames(d) %in% rownames(le)[rownames(le) %in% res.l.le$gene_id],], 2, 6, "corr-L-LE")
	clust_boot(d[rownames(d) %in% rownames(leku)[rownames(leku) %in% res.l.leku$gene_id],], 2, 6, "corr-L-LEKU")
	clust_boot(d[rownames(d) %in% rownames(le.leku)[rownames(le.leku) %in% res.le.leku$gene_id],], 2, 6, "corr-LE-LEKU")

	plot_bootstrap_data("corr-L-LE")
	plot_bootstrap_data("corr-L-LEKU")
	plot_bootstrap_data("corr-LE-LEKU")
}

get_clust <- function() {
	clust1 <- get_clust_genes("corr-L-LE", 3)
	clust2 <- get_clust_genes("corr-L-LEKU", 3)
	clust3 <- get_clust_genes("corr-LE-LEKU", 3)

	clusts <- list(clust1$partition, clust2$partition, clust3$partition)

	plot_all_clusts(clusts, "corr", c(1, 2, 3))
	return(clusts)
}

go_clust <- function(clusts) {
	ind <- 1
	for(j in clusts) {
		for (i in 1:length(unique(j))) {
			GO(names(j[j == i]), all.genes, paste0("a66-", c(1, 2, 3)[ind], "-", i), 0.05)
		}
		# print(inds[ind])
		ind <- ind + 1
	}
}

plot_groups <- function() {
	# plot polysome data for the smallest gene sets
	plot_genes(setG.filt, "setG")
	plot_genes(setG.filt[setG.filt %in% res.l.le$gene_id], "setG-in-L-LE")
	plot_genes(setG.filt[setG.filt %in% res.le.leku$gene_id], "setG-in-LE-LEKU")
	plot_genes(setG.filt[setG.filt %in% res.l.leku$gene_id], "setG-in-L-LEKU")

	plot_genes(setH.filt, "setH")
	plot_genes(setH.filt[setH.filt %in% res.l.le$gene_id], "setH-in-L-LE")
	plot_genes(setH.filt[setH.filt %in% res.le.leku$gene_id], "setH-in-LE-LEKU")
	plot_genes(setH.filt[setH.filt %in% res.l.leku$gene_id], "setH-in-L-LEKU")
}

go_groups <- function() {
	GO(setA.filt[setA.filt %in% res.l.le$gene_id], all.genes, "A-L-LE", 0.05)
	GO(setA.filt[setA.filt %in% res.l.leku$gene_id], all.genes, "A-L-LEKU", 0.05)
	GO(setA.filt[setA.filt %in% res.le.leku$gene_id], all.genes, "A-LE-LEKU", 0.05)

	GO(setB.filt[setB.filt %in% res.l.le$gene_id], all.genes, "B-L-LE", 0.05)
	GO(setB.filt[setB.filt %in% res.l.leku$gene_id], all.genes, "B-L-LEKU", 0.05)
	GO(setB.filt[setB.filt %in% res.le.leku$gene_id], all.genes, "B-LE-LEKU", 0.05)

	GO(setC.filt[setC.filt %in% res.l.le$gene_id], all.genes, "C-L-LE", 0.05)
	GO(setC.filt[setC.filt %in% res.l.leku$gene_id], all.genes, "C-L-LEKU", 0.05)
	GO(setC.filt[setC.filt %in% res.le.leku$gene_id], all.genes, "C-LE-LEKU", 0.05)

	GO(setD.filt[setD.filt %in% res.l.le$gene_id], all.genes, "D-L-LE", 0.05)
	GO(setD.filt[setD.filt %in% res.l.leku$gene_id], all.genes, "D-L-LEKU", 0.05)
	GO(setD.filt[setD.filt %in% res.le.leku$gene_id], all.genes, "D-LE-LEKU", 0.05)

	GO(setE.filt[setE.filt %in% res.l.le$gene_id], all.genes, "E-L-LE", 0.05)
	GO(setE.filt[setE.filt %in% res.l.leku$gene_id], all.genes, "E-L-LEKU", 0.05)
	GO(setE.filt[setE.filt %in% res.le.leku$gene_id], all.genes, "E-LE-LEKU", 0.05)

	GO(setF.filt[setF.filt %in% res.l.le$gene_id], all.genes, "F-L-LE", 0.05)
	GO(setF.filt[setF.filt %in% res.l.leku$gene_id], all.genes, "F-L-LEKU", 0.05)
	GO(setF.filt[setF.filt %in% res.le.leku$gene_id], all.genes, "F-LE-LEKU", 0.05)

	GO(setG.filt[setG.filt %in% res.l.le$gene_id], all.genes, "G-L-LE", 0.05)
	GO(setG.filt[setG.filt %in% res.l.leku$gene_id], all.genes, "G-L-LEKU", 0.05)
	GO(setG.filt[setG.filt %in% res.le.leku$gene_id], all.genes, "G-LE-LEKU", 0.05)

	GO(setH.filt[setH.filt %in% res.l.le$gene_id], all.genes, "H-L-LE", 0.05)
	GO(setH.filt[setH.filt %in% res.l.leku$gene_id], all.genes, "H-L-LEKU", 0.05)
	GO(setH.filt[setH.filt %in% res.le.leku$gene_id], all.genes, "H-LE-LEKU", 0.05)
}

cluster_groups <- function() {
	# clustering
	d <- readRDS("files/data-matrix-norm.rds")
	rownames(d) <- d[,1]
	d <- as.matrix(d[,2:dim(d)[2]])

	clust_boot(d[rownames(d) %in% setA.filt[setA.filt %in% res.l.le$gene_id],], 2, 6, "clusts-A-L-LE")
	clust_boot(d[rownames(d) %in% setA.filt[setA.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-A-L-LEKU")
	clust_boot(d[rownames(d) %in% setA.filt[setA.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-A-LE-LEKU")

	clust_boot(d[rownames(d) %in% setB.filt[setB.filt %in% res.l.le$gene_id],], 2, 6, "clusts-B-L-LE")
	clust_boot(d[rownames(d) %in% setB.filt[setB.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-B-L-LEKU")
	clust_boot(d[rownames(d) %in% setB.filt[setB.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-B-LE-LEKU")

	clust_boot(d[rownames(d) %in% setC.filt[setC.filt %in% res.l.le$gene_id],], 2, 6, "clusts-C-L-LE")
	clust_boot(d[rownames(d) %in% setC.filt[setC.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-C-L-LEKU")
	clust_boot(d[rownames(d) %in% setC.filt[setC.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-C-LE-LEKU")

	clust_boot(d[rownames(d) %in% setD.filt[setD.filt %in% res.l.le$gene_id],], 2, 6, "clusts-D-L-LE")
	clust_boot(d[rownames(d) %in% setD.filt[setD.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-D-L-LEKU")
	clust_boot(d[rownames(d) %in% setD.filt[setD.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-D-LE-LEKU")

	clust_boot(d[rownames(d) %in% setE.filt[setE.filt %in% res.l.le$gene_id],], 2, 6, "clusts-E-L-LE")
	clust_boot(d[rownames(d) %in% setE.filt[setE.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-E-L-LEKU")
	clust_boot(d[rownames(d) %in% setE.filt[setE.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-E-LE-LEKU")

	clust_boot(d[rownames(d) %in% setF.filt[setF.filt %in% res.l.le$gene_id],], 2, 6, "clusts-F-L-LE")
	clust_boot(d[rownames(d) %in% setF.filt[setF.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-F-L-LEKU")
	clust_boot(d[rownames(d) %in% setF.filt[setF.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-F-LE-LEKU")

	clust_boot(d[rownames(d) %in% setG.filt[setG.filt %in% res.l.le$gene_id],], 2, 6, "clusts-G-L-LE")
	clust_boot(d[rownames(d) %in% setG.filt[setG.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-G-L-LEKU")
	clust_boot(d[rownames(d) %in% setG.filt[setG.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-G-LE-LEKU")

	clust_boot(d[rownames(d) %in% setH.filt[setH.filt %in% res.l.le$gene_id],], 2, 6, "clusts-H-L-LE")
	clust_boot(d[rownames(d) %in% setH.filt[setH.filt %in% res.l.leku$gene_id],], 2, 6, "clusts-H-L-LEKU")
	clust_boot(d[rownames(d) %in% setH.filt[setH.filt %in% res.le.leku$gene_id],], 2, 6, "clusts-H-LE-LEKU")

	plot_bootstrap_data("clusts-A-L-LE")
	plot_bootstrap_data("clusts-A-L-LEKU")
	plot_bootstrap_data("clusts-A-LE-LEKU")
	plot_bootstrap_data("clusts-B-L-LE")
	plot_bootstrap_data("clusts-B-L-LEKU")
	plot_bootstrap_data("clusts-B-LE-LEKU")
	plot_bootstrap_data("clusts-C-L-LE")
	plot_bootstrap_data("clusts-C-L-LEKU")
	plot_bootstrap_data("clusts-C-LE-LEKU")
	plot_bootstrap_data("clusts-D-L-LE")
	plot_bootstrap_data("clusts-D-L-LEKU")
	plot_bootstrap_data("clusts-D-LE-LEKU")
	plot_bootstrap_data("clusts-E-L-LE")
	plot_bootstrap_data("clusts-E-L-LEKU")
	plot_bootstrap_data("clusts-E-LE-LEKU")
	plot_bootstrap_data("clusts-F-L-LE")
	plot_bootstrap_data("clusts-F-L-LEKU")
	plot_bootstrap_data("clusts-F-LE-LEKU")
	plot_bootstrap_data("clusts-G-L-LE")
	plot_bootstrap_data("clusts-G-L-LEKU")
	plot_bootstrap_data("clusts-G-LE-LEKU")
	plot_bootstrap_data("clusts-H-L-LE")
	plot_bootstrap_data("clusts-H-L-LEKU")
	plot_bootstrap_data("clusts-H-LE-LEKU")
}
