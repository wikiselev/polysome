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
        return(p.vals)
}
