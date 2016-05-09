#
# ggplot versions of of vizualization functions
#

gg.global <- function(ssc) #ssc is an ssClass (skelesim class) object
{
    gdf <- globalDF(ssc@analysis.results)

    p <- ggplot(data=gdf,aes(x=Locus,y=value)) +
        facet_wrap(~statistic, scales="free") +
        geom_violin()+geom_jitter(width=0.1)
    
    print(p)
}

gg.locus <- function(ssc)
{
    ldf <- locusDF(ssc@analysis.results)
    l <- ggplot(data=ldf,aes(x=locus,y=value)) + geom_violin()+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")
    print(l)
    p <- ggplot(data=ldf,aes(x=pop,y=value)) + geom_violin()+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")
    print(p)
    
}

gg.pairwise <- function(ssc)
    {
        pwdf <- pairwiseDF(ssc@analysis.results)
        pwdf.mn <- pwdf %>% group_by(pop1,pop2,rep,statistic) %>% summarise(value=mean(value))

        plts <- list()
        for (s in unique(pwdf$statistic))
            {
                df <- pwdf.mn[pwdf.mn$statistic==s,]
                p <- ggplot(data=df,aes(pop1,pop2)) 
                p <- p + ggtitle(s)
                p <- p + geom_tile(aes(fill=value), colour="white")
                p <- p + scale_fill_gradient(low="white",high="steelblue")
                plts[[length(plts)+1]] <- p
            }
        multiplot(plotlist=plts,cols=4)
    }


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
                                        # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
                                        # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
                                        # Make the panel
                                        # ncol: Number of columns of plots
                                        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
                                        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
                                        # Make each plot, in the correct location
        for (i in 1:numPlots) {
                                        # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
        }
    }
						  }	
