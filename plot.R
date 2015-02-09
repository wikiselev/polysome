ensembl_id_to_hgnc_symbol <- function(genes) {
        mart <- readRDS("files/ann-table.rds")
        res <- unique(mart[mart$gene_id %in% genes])
        res <- res[res$gene_name != "", ]
        return(res)
}

hgnc_symbol_to_ensembl_id <- function(genes) {
        mart <- readRDS("files/ann-table.rds")
        res <- unique(mart[mart$gene_name %in% genes])
        # res <- res[res$ensembl_gene_id != "", ]
        return(res)
}

get_plot_data <- function(genes) {
        # if genes are represented by ensembl ids
        if(grepl("ENSMUSG000", genes[1])){
        res <- ensembl_id_to_hgnc_symbol(genes)
        } else {
        res <- hgnc_symbol_to_ensembl_id(genes)
        }
        # colnames(res) <- c("id", "name")
        
        plot.data <- readRDS("files/plot-table.rds")
        setkey(plot.data, "gene_id")
        setkey(res, "gene_id")
        plot.data <- plot.data[plot.data$gene_id %in% res$gene_id, ]
        plot.data <- merge(plot.data, res)
        return(plot.data)
}

plot_genes <- function(gene.names, file.name){
        plot.data <- get_plot_data(gene.names)
        
        limits <- aes(ymax = mean_value + sd_value, ymin = mean_value - sd_value)
        
        p <- ggplot(plot.data, aes(pf, mean_value, group = cond, color = cond)) +
                geom_line() + facet_wrap(~ gene_name, scale = "free_y") +
                geom_errorbar(limits, width = 0.25) +
        theme_bw()
        p
        
        pdf(file = paste0("plots/", file.name, ".pdf"),
          width = (sqrt(length(gene.names)) + 2)*3,
          height = sqrt(length(gene.names))*3)
        print(p)
        dev.off()
}
