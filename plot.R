get_plot_data <- function(genes) {
        plot.data <- readRDS("files/plot-table.rds")
        # if genes are represented by ensembl ids
        if(grepl("ENSMUSG000", genes[1])){
                return(plot.data[ensembl_gene_id %in% genes])
        } else {
                return(plot.data[gene_name %in% genes])
        }
}

plot_genes <- function(gene.names, file.name){
        plot.data <- get_plot_data(gene.names)
        plot.data$pf <- as.character(plot.data$pf)
        plot.data$pf <- as.numeric(plot.data$pf)
        limits <- aes(ymax = mean_value + sd_value, ymin = mean_value - sd_value)
        p <- ggplot(plot.data, aes(pf, mean_value, group = cond, color = cond)) +
                geom_line() + facet_wrap(~ gene_name, scale = "free_y") +
                geom_errorbar(limits, width = 0.25) +
                geom_vline(xintercept = c(7.5, 10.5)) +
                labs(x = "Polysome fraction", y = "Read counts") +
                theme_bw()
        pdf(file = paste0("plots/", file.name, ".pdf"),
          width = (sqrt(length(gene.names)) + 2)*3,
          height = sqrt(length(gene.names))*3)
        print(p)
        dev.off()
}
