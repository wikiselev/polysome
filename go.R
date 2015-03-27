GO <- function(selected.genes, all.genes, name, p.val) {
  system(paste0("rm -r GO/", name))
  system(paste0("mkdir GO/", name))

  geneList <- factor(as.integer(all.genes %in% selected.genes))
  names(geneList) <- all.genes

  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
    annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

  res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

  res.table <- GenTable(GOdata, classicFisher = res,
    topNodes = length(usedGO(GOdata)))
  res.table$classicFisher <- as.numeric(res.table$classicFisher)

  res.table <- res.table[res.table$classicFisher < p.val,]

  genes.ann <- genesInTerm(GOdata, whichGO=res.table$GO.ID)
  genes.ann <- lapply(genes.ann, function(x) {x[x %in% selected.genes]})
  saveRDS(genes.ann, paste0("GO/", name, "/genes-ann-BP.rds"))

  # write.table(res.table, file = paste0("../pip3-rna-seq-output/GO/", name, "/BP-Fisher.txt"),
  #   quote = F, row.names = F, col.names = F, sep = "\t")

  write.table(res.table[,c(1,6)], file = paste0("GO/", name, "/BP-Fisher-revigo.txt"),
    quote = F, row.names = F, col.names = F, sep = "\t")

  # printGraph(GOdata, res, firstSigNodes = 20,
  #   fn.prefix = paste0("../pip3-rna-seq-output/GO/", name, "/Fisher"),
  #   useInfo = "def", pdfSW = TRUE)

  system(paste0("python revigo.py BP-Fisher-revigo.txt ", name))
}
