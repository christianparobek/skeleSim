#
# functions to provide an abstraction layer between the 'ssClass@analysis.results' object and dataframes that
# ggplot can use
#

#this one returns the global results
globalDF <- function(ar)
{
    if ("Global" %in% names(ar[[1]]))
    {
        gar <- ar[[1]][['Global']]
        gdf <- as.data.frame(as.table(gar))
        names(gdf) <- c("Locus","statistic","rep","value")
        gdf[complete.cases(gdf),] #get rid of statistics that are NAs
    } else {NULL}
}

locusDF <- function(ar)
{
    if ("Locus" %in% names(ar[[1]]))
    {
        lar <- ar[[1]][["Locus"]]
        ldf <- as.data.frame(as.table(lar))
        ploc <- strsplit(as.character(ldf$Var1)," ")
        population <- sapply(ploc,function(x){if(length(x)>1){x[2]} else {"overall"}})
        locus <- gsub("Locus_","",gsub("_Sample","",sapply(ploc,function(x){x[1]})))
        data.frame(locus=locus,pop=population,rep=ldf$Var3,statistic=ldf$Var2,value=ldf$Freq)
    } else {NULL}
}

pairwiseDF <- function(ar)
    {
        if ("Pairwise" %in% names(ar[[1]]))
            {
                pwar <- ar[[1]][["Pairwise"]]
                pwdf <- as.data.frame(as.table(pwar))
                pwdf$Var1 <- gsub("^ ","",gsub(" +"," ",gsub("_"," ",gsub("Sample|Locus","",as.character(pwdf$Var1)))))
                pop1 <- sapply(strsplit(pwdf$Var1," "),function(x){x[1]})
                pop2 <- sapply(strsplit(pwdf$Var1," "),function(x){x[2]})
                locus <- sapply(strsplit(pwdf$Var1," "),function(x){x[3]})
                df <- data.frame(locus=locus,pop1=pop1,pop2=pop2,rep=pwdf$Var3,statistic=pwdf$Var2,value=pwdf$Freq)
                df <- df[complete.cases(df),]
                df
            } else {NULL}
    }
