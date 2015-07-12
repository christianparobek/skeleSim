#AES 6-28-15
#
# these functions are designed as helpers and interfaces for creating simcoal histories.
# the basic data structure is a data frame that has the same columns as the elements in a
# simcoal history entry:
#  time, source, sink, migrants, new size, growth rate, migr. matrix
#

create.new.history <- function(npop=3,
                               deepestPopcoal=10000,
                               method=c("allroot","exponential")[2]
                               )
    {
        if (npop<2) {stop("function needs at least 2 populations")}
        history <- data.frame(
            time=NA, source=NA, sink=NA, migrants=NA, new.size=NA, growth.rate=NA, migr.matrix=NA
            )

        history <- history[rep(1,npop-1),]
        history$time <- deepestPopcoal
        history$migr.matrix <- 0
        history$growth.rate <- 1
        history$migrants <- 1
        history$growth.rate <- 0
        history$new.size <- 1
        history$sink <- 0
        history$source <- 1:dim(history)[1]

        if (method=="exponential")
            {
                history$time <- round(rexp(npop-1,rate=1/(deepestPopcoal/2)))
            }
        
        history
    }



simcoal.history.plot <- function(history
                         )
    {
        npop <- max(c(history$source,history$sink))+1
        pops <- 0:(npop-1)

        plot(1~1,type="n",
             ylim=c(0,max(history$time)*1.1),
             xlim=c(-1,npop),axes=F,xlab="population",ylab="time")
        axis(1,labels=pops+1,at=0:(npop-1))
        axis(2)
        for (p in 0:(npop-1))
            {
                x.source <- p
                p.hist <- history[history$source==p,]
                if (dim(p.hist)[1]>0)
                    for (e in 1:dim(p.hist)[1])
                    {
                        x.sink <- p.hist$sink[e]
                        y <-  p.hist$time[e]
                        points(x=c(x.source,x.source),
                               y=c(0,y), type="l")
                        arrows(x.source,y,x.sink,y)
                    }
            }
        ###now deal with the sink-only situations (should be one at most)
        spops <- unique(history$sink[!history$sink%in%history$source])
        for (i in spops)
            {
                arrows(i,0,i,max(history$time)*1.1)
            }
    }

simcoal.history.change <- function(history,
                           click,
                           dblclick
                           )
{
    oldhist <- history
    src <- round(click$x)
    if (src %in% c(history$source,history$sink))
        {
            history <- history[history$source!=src,]
            history <- history[c(1:dim(history)[1],1),]
            row <- dim(history)[1]
            history[row,] <- c(round(dblclick$y),
                               src,
                               round(dblclick$x),
                               1,
                               1,
                               0,
                               0)
            history
        }
    if (is.history(history)) history else oldhist
}


is.history <- function(history)
    {
        err <- F
        msg <- NULL
        ##check for names
        if (prod(names(history) %in% c("time", "source", "sink", "migrants", "new.size", "growth.rate", "migr.matrix"))!=1) err <- T
        if (prod(c("time", "source", "sink", "migrants", "new.size", "growth.rate", "migr.matrix")%in%names(history))!=1) err <- T
        if (err) msg <- c(msg,"problem with names")

        ######## this code is essentially worthless, but the problem
        ######## needs solving
        ##check to make sure that everybody coalesces
        ##first run through the sources
        last.sink <- history$sink[history$time==max(history$time)]
        if (length(unique(last.sink))>1) #multiple sinks at the "root"
            {
                err <- T
                msg <- c(msg,"multiple sinks at deepest time")
            }
        last.sink <- last.sink[1]
        
        history$coal <- FALSE
        for  (s in 1:dim(history))
            {
                snk <- history[s,"sink"]
                
            }
        
        if (err)
            {
                FALSE
            } else {
                TRUE
            }
    }
