#
# ggplot versions of of vizualization functions
#


viz.scenarios <- function(ssc)
{
    if (!is.null(req(ssc@analysis.results)))
        length(ssc@analysis.results)
    else
        NULL
}

global.stats <- function(ssc)
{
    if (!is.null(ssc@analysis.results))
        unique(sort(as.character(globalDF(ssc)$statistic)))
}

gg.global <- function(ssc,stats=global.stats(ssc),scenario=1) #ssc is an ssClass (skelesim class) object
{
    gdf <- globalDF(ssc)
    gdf <- gdf[gdf$scenario%in%scenario,]
    gdf <- gdf[gdf$statistic%in%stats,]
    gdf <- gdf[gdf$Locus!="Overall",]

    p <- ggplot(data=gdf,aes(x=Locus,y=value)) +
        facet_wrap(~statistic, scales="free") +
        geom_violin()+geom_jitter(width=0.1)
    
    print(p)
}

gg.global.scmp <- function(ssc,stats=global.stats(ssc)) #ssc is an ssClass (skelesim class) object
{

    gdf <- globalDF(ssc)
    gdf <- gdf[gdf$Locus!="Overall",]
    gdf <- filter(gdf,statistic%in%stats) %>% filter(Locus!="Overall") %>%
        group_by(statistic,scenario,rep) %>% summarise(value=mean(value,na.rm=T))
    gdf$scenario <- as.factor(gdf$scenario)
    p <- ggplot(data=gdf,aes(x=scenario,y=value)) +
        facet_wrap(~statistic, scales="free") +
        geom_violin()+geom_jitter(width=0.1)
    
    print(p)
}



df.global <- function(ssc,stats=global.stats(ssc)) #ssc is an ssClass (skelesim class) object
{
    gdf <- globalDF(ssc)
    gdf[gdf$statistic%in%stats,]
}

locus.stats <- function(ssc)
{
    if (!is.null(ssc@analysis.results))
        unique(sort(as.character(locusDF(ssc)$statistic)))
}

gg.locus <- function(ssc,stats=locus.stats(ssc),scenario=1)
{
    ldf <- locusDF(ssc)
    ldf <- ldf[ldf$statistic%in%stats,]
    ldf <- ldf[ldf$scenario%in%scenario,]
    ldf <- ldf[ldf$pop!="overall",]
    l <- ggplot(data=ldf,aes(x=locus,y=value)) + geom_violin()+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")
    
    p <- ggplot(data=ldf,aes(x=pop,y=value)) + geom_violin()+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")
    multiplot(plotlist=list(l,p))
}


gg.locus.scmp <- function(ssc,stats=locus.stats(ssc))
{
    ldf <- locusDF(ssc)
    ldf <- ldf[ldf$statistic%in%stats,]
    ldf <- ldf[ldf$pop!="overall",]
    ldf$scenario <- as.factor(ldf$scenario)
    ldf <- group_by(ldf,pop,locus,scenario,statistic)%>%summarise(value=mean(value,na.rm=T))
    l <- ggplot(data=ldf,aes(x=scenario,y=value,group=locus,col=locus)) + geom_violin(aes(group=scenario))+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")+stat_summary(fun.y="mean",geom="line")

    p <- ggplot(data=ldf,aes(x=scenario,y=value,group=pop,col=pop)) + geom_violin(aes(group=scenario))+
        geom_jitter(width=0.15)+facet_wrap(~statistic,scales="free")+stat_summary(fun.y="mean",geom="line")

    p
    multiplot(plotlist=list(l,p))
}

df.locus <- function(ssc,stats=locus.stats(ssc))
{
    ldf <- locusDF(ssc)
    ldf <- ldf[ldf$statistic%in%stats,]
}

pairwise.stats <- function(ssc)
{
        if (!is.null(ssc@analysis.results))
            unique(sort(as.character(pairwiseDF(ssc)$statistic)))
}

gg.pairwise <- function(ssc,stats=pairwise.stats(ssc),scenario=1)
    {
        pwdf <- pairwiseDF(ssc)
        pwdf <- pwdf[pwdf$statistic%in%stats,]
        pwdf <- pwdf[pwdf$scenario%in%scenario,]
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

gg.pairwise.scmp <- function(ssc,stats=pairwise.stats(ssc))
    {
        pwdf <- pairwiseDF(ssc)
        pwdf <- pwdf[pwdf$statistic%in%stats,]
        pwdf$scenario <- as.factor(pwdf$scenario)
        
        pwdf.mn <- pwdf %>% group_by(pop1,pop2,statistic,scenario) %>% summarise(value=mean(value))

        plts <- list()
        for (s in unique(pwdf$statistic))
#            for (scen in unique(pwdf.mn$scenario))
            {
                df <- pwdf.mn[pwdf.mn$statistic==s,]
                p <- ggplot(data=df,aes(pop1,pop2)) 
                p <- p + ggtitle(s)
                p <- p + geom_tile(aes(fill=value), colour="white")
                p <- p + scale_fill_gradient(low="white",high="steelblue")
                p <- p + facet_wrap(~scenario,scales="free")
                plts[[length(plts)+1]] <- p
            }
        multiplot(plotlist=plts,cols=ifelse(length(plts)>=4,4,length(plts)))
    }

df.pairwise <- function(ssc,stats=pairwise.stats(ssc))
{
        pwdf <- pairwiseDF(ssc)
        pwdf <- pwdf[pwdf$statistic%in%stats,]
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

vizMegaPlot <- function(ssc)
{

    for (s in 1:length(ssc@scenarios)) gg.global(ssc=ssc,scenario=s)
    gg.global.scmp(ssc=ssc)
    for (s in 1:length(ssc@scenarios))  gg.locus(ssc=ssc,scenario=s)
    gg.locus.scmp(ssc=ssc)
    for (s in 1:length(ssc@scenarios))  gg.pairwise(ssc=ssc,scenario=s)
    gg.pairwise.scmp(ssc=ssc)

}
