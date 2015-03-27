process_data <- function() {
	d <- read.table("original-data/p28_expression_genes_edited.txt", sep = "\t", header = T)

	# remove Poly_L_2_17 sample
	d <- d[,!grepl("LUL16_e14", colnames(d))]

	# rename other samples in LUL16
	for(i in 15:40) {
		colnames(d) <- gsub(paste0("LUL16_e", i), paste0("LUL16_e", i-1), colnames(d))
	}
        colnames(d)[4] <- "ensembl_gene_id"
	d <- d[,c(4,13:dim(d)[2])]
	saveRDS(d, "files/data-matrix.rds")

	# normalize by sample library sizes (using DESeq2 library)
	colData <- data.frame(condition = colnames(d)[2:157])
	cds <- DESeqDataSetFromMatrix(countData = d[, c(2:157)], colData = data.frame(condition = colnames(d)[2:157]), design = ~ condition)
	cds <- estimateSizeFactors(cds)
	d1 <- as.data.frame(counts(cds, normalized = TRUE))
	colnames(d1) <- colnames(d)[2:157]
	d1 <- cbind(d$ensembl_gene_id, d1)
	colnames(d1)[1] <- "ensembl_gene_id"
	saveRDS(d1, "files/data-matrix-norm.rds")

	# prepare a table for clustering
	tmp <- as.data.frame(t(scale(t(d1[,2:157]), center = T, scale = T)))
	d2 <- cbind(d1[,1], tmp)
	colnames(d2)[1] <- "ensembl_gene_id"
	saveRDS(d2, "files/data-matrix-norm-scaled.rds")

	# prepare gene name/length annotation table
	ann <- read.table("annotations/genes-GRCm38.p3.txt", sep = "\t", header = T)
	colnames(ann) <- c("ensembl_gene_id", "gene_start", "gene_end", "gene_name")
	ann$gene_length <- abs(ann$gene_end - ann$gene_start)
	ann <- ann[,c(1,4,5)]
	saveRDS(as.data.table(ann), "files/ann-table.rds")

	# prepare a data table for plotting
	# average replicates - find mean and sd
	d3 <- rename_samples_plot(d1)
	saveRDS(d3, "files/data-table.rds")
	d4 <- d3[,list(mean_value=mean(value), sd_value=sd(value)),by=c("ensembl_gene_id", "cond", "pf")]
	d4 <- merge(as.data.frame(d4), ann)
	saveRDS(d4, "files/plot-table.rds")

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
        diff_expr_two_time_courses_cond("L", "LE")
        diff_expr_two_time_courses_cond("LE", "LEKU")
        diff_expr_two_time_courses_cond("L", "LEKU")
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

        # import gene lengths annotations
        gene.len <- readRDS("files/ann-table.rds")

	# old polysome differential expression gene sets - now we first split
	# polysome data by fraction and then do differential expression analysis
	# in each fraction category separately - see above
	res.l.le.deseq <- get_diff_expr("L.LE-polysome-cond", 0.01)
	res.l.leku.deseq <- get_diff_expr("L.LEKU-polysome-cond", 0.01)
	res.le.leku.deseq <- get_diff_expr("LE.LEKU-polysome-cond", 0.01)
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
	all.genes <<- unique(c(res.l.le.deseq$ensembl_gene_id,
	                       res.l.leku.deseq$ensembl_gene_id,
	                       res.le.leku.deseq$ensembl_gene_id,
	                       rownames(le),
	                       rownames(leku),
	                       rownames(le.leku)))
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

for_report <- function() {
        t <- res.l.le.deseq[!(res.l.le.deseq$ensembl_gene_id %in% rownames(le)),]
        genes <- t$ensembl_gene_id
        t1 <- posthoc.l.le[ensembl_gene_id %in% genes]
        setkey(t1, "ensembl_gene_id")
        min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
        setkey(min.pvals, "ensembl_gene_id")
        t1 <- t1[min.pvals]
        t1 <- t1[sig.pf == min.pval]
        t1 <- t1[min.pval < 0.01]
        GO(unique(t1[pf < 8, ensembl_gene_id]), all.genes, "L-LE-monosome-change", 0.05)
        GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "L-LE-light-change", 0.05)
        GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "L-LE-heavy-change", 0.05)
        GO(unique(t1[pf == 16, ensembl_gene_id]), all.genes, "L-LE-heavy16-change", 0.05)
        t1 <- t1[order(min.pval)]
        plot_genes(head(unique(t1[pf < 8, ensembl_gene_id]), 20), "L-LE-monosome-change")
        plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "L-LE-light-change")
        plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "L-LE-heavy-change")
        plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "L-LE-heavy16-change")

        t <- res.l.leku.deseq[!(res.l.leku.deseq$ensembl_gene_id %in% rownames(leku)),]
        genes <- t$ensembl_gene_id
        t1 <- posthoc.l.leku[ensembl_gene_id %in% genes]
        setkey(t1, "ensembl_gene_id")
        min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
        setkey(min.pvals, "ensembl_gene_id")
        t1 <- t1[min.pvals]
        t1 <- t1[sig.pf == min.pval]
        t1 <- t1[min.pval < 0.01]
        GO(unique(t1[pf < 8, ensembl_gene_id]), all.genes, "L-LEKU-monosome-change", 0.05)
        GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "L-LEKU-light-change", 0.05)
        GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "L-LEKU-heavy-change", 0.05)
        GO(unique(t1[pf == 16, ensembl_gene_id]), all.genes, "L-LEKU-heavy16-change", 0.05)
        t1 <- t1[order(min.pval)]
        plot_genes(head(unique(t1[pf < 8, ensembl_gene_id]), 20), "L-LEKU-monosome-change")
        plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "L-LEKU-light-change")
        plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "L-LEKU-heavy-change")
        plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "L-LEKU-heavy16-change")

        t <- res.le.leku.deseq[!(res.le.leku.deseq$ensembl_gene_id %in% rownames(le.leku)),]
        genes <- t$ensembl_gene_id
        t1 <- posthoc.le.leku[ensembl_gene_id %in% genes]
        setkey(t1, "ensembl_gene_id")
        min.pvals <- t1[,list(min.pval = min(sig.pf)),by = "ensembl_gene_id"]
        setkey(min.pvals, "ensembl_gene_id")
        t1 <- t1[min.pvals]
        t1 <- t1[sig.pf == min.pval]
        t1 <- t1[min.pval < 0.01]
        GO(unique(t1[pf < 8, ensembl_gene_id]), all.genes, "LE-LEKU-monosome-change", 0.05)
        GO(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), all.genes, "LE-LEKU-light-change", 0.05)
        GO(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), all.genes, "LE-LEKU-heavy-change", 0.05)
        GO(unique(t1[pf == 16, ensembl_gene_id]), all.genes, "LE-LEKU-heavy16-change", 0.05)
        t1 <- t1[order(min.pval)]
        plot_genes(head(unique(t1[pf < 8, ensembl_gene_id]), 20), "LE-LEKU-monosome-change")
        plot_genes(head(unique(t1[pf >= 8 & pf <= 10, ensembl_gene_id]), 20), "LE-LEKU-light-change")
        plot_genes(head(unique(t1[pf > 10 & pf != 16, ensembl_gene_id]), 20), "LE-LEKU-heavy-change")
        plot_genes(head(unique(t1[pf == 16, ensembl_gene_id]), 20), "LE-LEKU-heavy16-change")
}

cluster_correlations <- function() {
        # clustering
        d <- readRDS("files/data-matrix-norm-scaled.rds")
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
