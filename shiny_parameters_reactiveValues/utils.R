#
# R functions to help with shiny skelesim interface
#
make.interactive.matrix <- function(mat,rows=dim(mat)[1],cols=dim(mat)[2])
    {
        if (!is.matrix(mat)) {stop("mat passed to make.interactive.matrix() is not an actual matrix")}
        retmat <- matrix("",rows,cols)
        for (row in 1:rows)
            for (col in 1:cols)
                {
                    intxt <- paste0("<input id='r",row,"c",col,
                                    "' class='input-tiny' type='number' value='",
                                    mat[row,col],"'>")
                    retmat[row,col] <-intxt
                }
        colnames(retmat) <- 1:cols
        rownames(retmat) <- 1:rows
        as.data.frame(retmat)
    }
