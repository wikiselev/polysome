diff_expr_pairwise <- function(cond1, cond2) {
  files <- list.files("htseq-count/")
  files <- files[grepl(paste0("_", cond1, ".", "_"), files) |
                 grepl(paste0("_", cond2, ".", "_"), files)]
  condition <- sapply(strsplit(files, "\\_"), "[[", 2)
  condition <- substr(condition, 1, nchar(condition)-1)
  samples <- data.frame(sample.name = files, file.name = files,
    condition = condition)

  # change levels
  samples$condition <-
    factor(samples$condition, levels=unique(condition))

  # create a DESeqDataSet with interaction term in the design formula
  cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples,
    directory = "htseq-count",
    design = ~ condition)

  # perform the main computation
  cds <- DESeq(cds)

  saveRDS(results(cds), paste0("files/diff-expr-", cond1, "-", cond2, ".rds"))
}

norm_data_matrix <- function() {
        d <- readRDS("files/data-matrix.rds")
        colData <- data.frame(condition = colnames(d)[2:157])
        cds <- DESeqDataSetFromMatrix(countData = d[, c(2:157)], colData = data.frame(condition = colnames(d)[2:157]), design = ~ condition)
        cds <- estimateSizeFactors(cds)
        d1 <- counts(cds, normalized = TRUE)
        rownames(d1) <- d$ensembl_gene_id
        colnames(d1) <- colnames(d)[2:157]
        saveRDS(d1, "files/data-matrix-norm.rds")
}


diff_expr_two_time_courses_cond <- function(cond1, cond2) {
        # cond1 and cond2 arguments are either "L", "LE", "LEKU"
        d <- readRDS("files/data-matrix.rds")
        # rename the columns
        cols <- colnames(d[,2:157])
        j <- 4
        for(i in 1:13){
        cols[grep(paste0("LUL1._e", i, "_"), cols)] <- paste0("L_", j)
        j <- j + 1
        }
        j <- 4
        for(i in 14:26){
        cols[grep(paste0("LUL1._e", i, "_"), cols)] <- paste0("LE_", j)
        j <- j + 1
        }
        j <- 4
        for(i in 27:39){
        cols[grep(paste0("LUL1._e", i, "_"), cols)] <- paste0("LEKU_", j)
        j <- j + 1
        }
        colnames(d)[2:157] <- cols

        countData <- d[,grepl(paste0(cond1, "_"), colnames(d)) | grepl(paste0(cond2, "_"), colnames(d))]
        rownames(countData) <- d$ensembl_gene_id
        ann <- colnames(countData)
        ann <- sapply(strsplit(ann, "\\."), "[[", 1)
        colnames(countData) <- ann
        colData <- data.frame(condition = sapply(strsplit(ann, "\\_"), "[[", 1),
                        pf = sapply(strsplit(ann, "\\_"), "[[", 2))

        cds <- DESeqDataSetFromMatrix(countData = countData,
        colData = colData, design = ~ condition + pf + condition:pf)

        # perform the main computation
        cds <- DESeq(cds)
        resultsNames(cds)

        # likelihood ratio test
        cdsLRT <- nbinomLRT(cds, reduced = ~ pf)
        saveRDS(results(cdsLRT), paste0("files/diff-expr-", cond1, ".", cond2, "-polysome-cond.rds"))
}

get_diff_expr <- function(name, padj) {
  d <- readRDS(paste0("files/diff-expr-", name, ".rds"))
  return(d[!is.na(d$padj) & d$padj < padj,])
}

venn <- function(list, is.weights, name) {
  v <- Venn(list)
  pdf(file = paste0("plots/venn-", name, ".pdf"))
  plot(v, doWeights = is.weights)
  dev.off()
}

venn_ellipses <- function(list, is.weights, name) {
  v <- Venn(list)
  pdf(file = paste0("plots/venn-", name, ".pdf"))
  plot(v, doWeights = is.weights, type = "ellipses")
  dev.off()
}
