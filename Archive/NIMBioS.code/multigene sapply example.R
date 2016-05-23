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

<<<<<<< HEAD
phist
=======

# A DNAbin example that behaves maybe
data(woodmouse)
genes <- list(gene1=woodmouse[,1:400], gene2=woodmouse[,2:402])
x <- new("multidna", genes)
x.g <- sequence2gtypes(x)
strata(x.g) <- c("A", "B")
results_gtype <- x.g


>>>>>>> 3252dce708b8b02ff116978934bd1c74ddbf76a5
