two_way_anova <- function(dat, ret.val) {
  res <- anova(lm(value ~ cond + pf + cond * pf, dat))
  if(ret.val == "int") {
    return(res[[5]][3])
  }
  if(ret.val == "cond") {
    return(res[[5]][1])
  }
  if(ret.val == "pf") {
    return(res[[5]][2])
  }
}

anova_results <- function(dat) {
  anova_res_int <- dat[, list(pval =
      two_way_anova(data.frame(value = value, cond = cond, pf = pf), "int")),
    by = gene_id]
  anova_res_cond <- dat[,
    list(pval =
      two_way_anova(data.frame(value = value, cond = cond, pf = pf), "cond")),
    by = gene_id]
  anova_res_pf <- dat[,
    list(pval =
      two_way_anova(data.frame(value = value, cond = cond, pf = pf), "pf")),
    by = gene_id]

  anova_res_int$padj <- p.adjust(anova_res_int$pval, method = "BH")
  anova_res_cond$padj <- p.adjust(anova_res_cond$pval, method = "BH")
  anova_res_pf$padj <- p.adjust(anova_res_pf$pval, method = "BH")
  return(list(int = anova_res_int, cond = anova_res_cond, pf = anova_res_pf))
}
