library(strataGdevel)
rm(list = ls())

data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x.g <- sequence2gtypes(x)
strata(x.g) <- c("A", "B")

phist <- sapply(locNames(x.g), function(n) {
  unname(statPhist(x.g[, n, ]))
})

phist
