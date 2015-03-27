clust_boot <- function(matrix, min.clust.num, max.clust.num, name) {
	system(paste0("mkdir clusts/", name))
	for(i in min.clust.num:max.clust.num){
		c.boot <- clusterboot(matrix, B = 100, bootmethod="boot",
	                       clustermethod = pamkCBI,
				   		   # clustermethod = kmeansCBI,
	                       krange = i, seed = 155)
		saveRDS(c.boot, paste0("clusts/", name, "/", i, ".rds"))
	}
}

plot_bootstrap_data <- function(name) {
	files <- list.files(paste0("clusts/", name))
	jaccard.stats <- data.frame(numeric(), factor(), factor())
	for(file in files){
		clust.num <- as.numeric(strsplit(file, "\\.")[[1]][1])
		dat <- readRDS(paste0("clusts/", name, "/", file))
		for(i in 1:clust.num) {
			jaccard.stats <- rbind(jaccard.stats,
								   cbind(dat$bootresult[i,],
								   		 i,
								   		 clust.num)
								   )
		}

	}
	colnames(jaccard.stats) <- c("jacc.sim", "clust.ind", "clust.num")
	jaccard.stats$jacc.sim <- as.numeric(jaccard.stats$jacc.sim)
	jaccard.stats$clust.ind <- factor(jaccard.stats$clust.ind)
	jaccard.stats$clust.num <- factor(jaccard.stats$clust.num)

	p <- ggplot(jaccard.stats, aes(clust.ind, jacc.sim)) +
		geom_boxplot() +
		facet_grid(clust.num ~ .) +
		geom_hline(yintercept=0.75, color = "red") +
		theme_bw()
	pdf(file = paste0("plots/clust-boots-", name, ".pdf"),
	  width = 7,
	  height = 10)
	print(p)
	dev.off()
}

get_clust_genes <- function(name, clust.ind) {
	t <- readRDS(paste0("clusts/", name, "/", clust.ind, ".rds"))
	return(t)
}

plot_all_clusts <- function(clusts, name, inds) {
	plot.data <- readRDS("files/plot-table-norm.rds")
	plot.data[,clust:=-1]

	ind <- 1
	clust.size <- NULL

	for (j in clusts) {
		clust.size <- rbind(clust.size, cbind(as.data.frame(table(j)), inds[ind]))
		for (i in 1:length(unique(j))) {
			genes <- names(j[j == i])
			plot.data[gene_id %in% genes, clust:=i]
			plot.data[gene_id %in% genes, set:=inds[ind]]
		}
		ind <- ind + 1;
	}

	colnames(clust.size) <- c("clust", "gene.num", "set")
	clust.size$val <- 0
	clust.size$sd <- 0
	clust.size$cond <- "wt"
	clust.size$time <- 0

	plot.data <- plot.data[clust != -1]
	# plot.data[,gr:=paste(id, cond, sep = "_")]

	plot.data <-
		plot.data[, list(val = mean(mean_value), sd = sd(mean_value)), by = c("cond", "pf", "clust", "set")]

	limits <- aes(ymax = val + sd, ymin = val - sd)

	p <- ggplot(plot.data, aes(pf, val, group = cond, color = cond)) +
		# geom_rect(data = clust.size, aes(fill = log10(gene.num)), xmin = -Inf,xmax = Inf,
	            # ymin = -Inf,ymax = Inf, alpha = 0.2) +
	    geom_line() +
	    geom_point() +
		geom_text(aes(x = 12, y = 1, label = gene.num),
	         data = clust.size, colour = "black") +
	    # facet_grid(set ~ clust, scale = "free_y")
	    facet_grid(set ~ clust, scale = "free") +
	    theme_bw()
	    # geom_errorbar(limits, width = 0.25)
	pdf(file = paste0("plots/clusts-", name, "-av-all.pdf"), w = 8, h = 6)
	print(p)
	dev.off()
}

