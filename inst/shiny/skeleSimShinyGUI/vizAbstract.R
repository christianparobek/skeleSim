#
# functions to provide an abstraction layer between the 'ssClass@analysis.results' object and dataframes that
# ggplot can use
#

#this one returns the global results
globalDF <- function(ssc)
{
    ar <- ssc@analysis.results
    if ("Global" %in% names(ar[[1]]))
    {
        ret <- do.call(rbind,lapply(1:length(ar),function(i)
        {
            gdf <- as.data.frame(as.table(ar[[i]][['Global']]));
            names(gdf) <- c("Locus","statistic","rep","value");
            gdf$scenario <- i
            gdf$Locus <- as.factor(gsub("Locus_","",as.character(gdf$Locus)))
            gdf
        }))
        ret[complete.cases(ret),] #get rid of statistics that are NAs
    } else {NULL}
}

####takes the analysis.results slot from a skelesim object and makes a data frame of locus information
locusDF <- function(ssc)
{
    ar <- ssc@analysis.results
    ret <- NULL
    if ("Locus" %in% names(ar[[1]]))
    {
#if (debug()) print("locusDF")
        ret <- do.call(rbind,lapply(1:length(ar),function(i)
        {
#if (debug()) print(i)
            if (ssc@scenarios[[i]]@locus.type %in% c("SNP","microsatellite"))
            {
                ldf <- as.data.frame(as.table(ar[[i]][['Locus']]))

#if (debug()) print(dim(ldf))
                
                ploc <- strsplit(as.character(ldf$Var1),"_")
                
                if (prod(sapply(ploc,function(x){x[1]=="Locus"})) == 1)
                    ploc <- lapply(ploc,function(x){x[-1]})
                
                population <- as.factor(gsub("Sample ","",sapply(ploc,function(x){if(length(x)>1){x[2]} else {NA}})))
                locus <- as.factor(gsub("L","",sapply(ploc,function(x){x[1]})))
               df <- data.frame(locus=locus,pop=population,rep=ldf$Var3,scenario=i,statistic=ldf$Var2,value=ldf$Freq)
            } else {  #data are sequence or other
                df <- NULL
            }
        }
        ))
        if (!is.null(ret))
            ret <- ret[complete.cases(ret),]
    }
    ret
}

####takes the analysis.results slot from a skelesim object and makes a data frame of pairwise information
pairwiseDF <- function(ssc)
    {
        ar <- ssc@analysis.results
        if ("Pairwise" %in% names(ar[[1]]))
        {
         ret <- do.call(rbind,lapply(1:length(ar),function(i)
             {
                pwdf <- as.data.frame(as.table(ar[[i]][["Pairwise"]]))
                pwdf$Var1 <- gsub("^ ","",gsub(" +"," ",gsub("_"," ",gsub("Sample|Locus","",as.character(pwdf$Var1)))))
                pop1 <- sapply(strsplit(pwdf$Var1," "),function(x){x[1]})
                pop2 <- sapply(strsplit(pwdf$Var1," "),function(x){x[2]})
                locus <- sapply(strsplit(pwdf$Var1," "),function(x){x[3]})
                df <- data.frame(locus=locus,pop1=pop1,pop2=pop2,rep=pwdf$Var3,scenario=i,statistic=pwdf$Var2,value=pwdf$Freq)
                df
             }))
         ret[complete.cases(ret),]
            } else {NULL}
    }
