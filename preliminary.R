cor_heatmap <- function(count.matrix, name) {
  # plot a correlation matrix from a count matrix
  # calculate pearson's correlation coefficients
  cor.matrix <- cor(count.matrix, method = "pearson")

  # plot correlation matrix in a file with 'name'
  pdf(file = paste0("plots/", name, "-correlations.pdf"), w = 15, h = 15)
  heatmap.2(cor.matrix, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
    col=bluered(9), breaks = 10, trace = "none", keysize = 0.5)
  dev.off()
}

rename_samples_plot <- function(d) {
	d2 <- as.data.table(melt(d))
	j <- 4
	for(i in 1:13){
		d2[grep(paste0("LUL1._e", i, "_"), variable),cond:="L"]
		d2[grep(paste0("LUL1._e", i, "_"), variable),pf:=j]	
		j <- j + 1;
	}
	j <- 4
	for(i in 14:26){
		d2[grep(paste0("LUL1._e", i, "_"), variable),cond:="LE"]
		d2[grep(paste0("LUL1._e", i, "_"), variable),pf:=j]	
		j <- j + 1;
	}
	j <- 4
	for(i in 27:39){
		d2[grep(paste0("LUL1._e", i, "_"), variable),cond:="LEKU"]
		d2[grep(paste0("LUL1._e", i, "_"), variable),pf:=j]	
		j <- j + 1;
	}
	d2 <- d2[,list(gene_id, value, cond, pf)]
	d2$pf <- factor(d2$pf)
	d2$cond <- factor(d2$cond)
	d2$gene_id <- factor(d2$gene_id)
	return(d2)
}

pca <- function() {
  d <- readRDS("files/data-matrix-norm.rds")

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

  tmp <- as.matrix(d[,2:157])
  cor_heatmap(tmp[,order(colnames(tmp))], "polysome")

  res <- prcomp(d[,2:157])

  pdf(file = "plots/pca.pdf", w=7, h=6)
  print(screeplot(res))
  dev.off()

  data <- as.data.frame(res$rotation[,1:3])
  data$cond <- rownames(data)
  data$cond <- unlist(lapply(strsplit(data$cond, "\\."), function(x){return(x[1])}))
  data$pf <- as.numeric(unlist(lapply(strsplit(data$cond, "_"), function(x){return(x[2])})))
  data$cond <- unlist(lapply(strsplit(data$cond, "_"), function(x){return(x[1])}))
  p <- ggplot(data, aes(PC1,PC2)) +
      geom_point(aes(shape = cond, color = as.factor(pf))) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "plots/pca12.pdf", w=7, h=6)
  print(p)
  dev.off()

  p <- ggplot(data, aes(PC2,PC3)) +
      geom_point(aes(shape = cond, color = as.factor(pf))) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "plots/pca23.pdf", w=7, h=6)
  print(p)
  dev.off()

  p <- ggplot(data, aes(PC1,PC3)) +
      geom_point(aes(shape = cond, color = as.factor(pf))) +
      scale_size(range = c(3, 6)) +
      theme_bw()
  pdf(file = "plots/pca13.pdf", w=7, h=6)
  print(p)
  dev.off()

  # in 3D
  # open3d(windowRect=c(100,100,700,700))
  # plot3d(res$rotation,xlab="PC1",ylab="PC2",zlab="PC3")
  # spheres3d(res$rotation, radius=0.01,col=rainbow(length(res$rotation[,1])))
  # grid3d(side="z", at=list(z=0))
  # text3d(res$rotation, text=rownames(res$rotation), adj=1.3)
  # rgl.postscript(file="plots/pca-3d.pdf", fmt="pdf")
}

count_matrix_from_files <- function(files, normalized) {
  # this function create a count matrix from raw htseq read counts files
  # it uses DESeq2 library to construct the matrix

  # create a design matrix for DESeq2
  sample.name <- sapply(strsplit(files, "_trimmed"), "[[", 1)
  condition <- sapply(strsplit(sample.name, "_"), "[[", 1)
  samples <- data.frame(sample.name = sample.name,
    file.name = paste0("htseq-count/", files),
    condition = condition)
  samples$condition <- factor(samples$condition,
    levels = unique(samples$condition))

  # import data from files using 'samples' design matrix
  cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples,
    directory = ".", design = ~ condition)

  # if normalized is TRUE then the count matrix is normalized by library size
  # using estimateSizeFactors function of DESeq2
  if (normalized) {
    cds <- estimateSizeFactors(cds)
    return(counts(cds, normalized = TRUE))
  } # otherwise return the raw matrix
  else {return(counts(cds))}
}
