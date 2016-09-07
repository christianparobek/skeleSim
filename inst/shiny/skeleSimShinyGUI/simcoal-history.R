#AES 6-28-15
#
# these functions are designed as helpers and interfaces for creating simcoal histories.
# the basic data structure is a data frame that has the same columns as the elements in a
# simcoal history entry:
#  time, source, sink, migrants, new size, growth rate, migr. matrix
#

create.new.history <- function(npop=3,
                               deepestPopcoal=5000,
                               method=c("allroot","exponential")[2]
                               )
    {
        if (npop>1) {
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
        } else {history <- NULL}
        history
    }



simcoal.history.plot <- function(history
                         )
    {
#                        print("in shp")
        if (!is.null(history))
            {

                history <- as.data.frame(history)
#                print(history)
                npop <- max(c(history[,2],history[,3]))+1
                pops <- 0:(npop-1)
                plot(x=0,y=0,type="n",
                     ylim=c(0,max(history[,1])*1.2),
                     xlim=c(-1,npop),axes=F,xlab="population",ylab="time")
                axis(1,labels=pops,at=0:(npop-1))
                axis(2)
                for (p in 0:(npop-1))
                    {
                        x.source <- p
                        p.hist <- history[history[,2]==p,]
                        if (dim(p.hist)[1]>0)
                            for (e in 1:dim(p.hist)[1])
                                {
                                    x.sink <- p.hist[e,3]
                                    y <-  p.hist[e,1]
                                    points(x=c(x.source,x.source),
                                           y=c(0,y), type="l")
                                    arrows(x.source,y,x.sink,y,length=0)
                                }
                    }
        ###now deal with the sink-only situations (should be one)
                spops <- history[history[,2]!=history[,3],]
                spopvec <- unique(spops[,3][!spops[,3]%in%spops[,2]])
                for (i in spopvec)
                    {
                        arrows(i,0,i,max(history[,1])*1.2)
                    }
            }
    }

simcoal.history.change <- function(history,
                           click,
                           dblclick
                           )
{
    if (!is.null(history))
        {
            if ((!is.null(click)) & (!is.null(dblclick)))
                {
                    oldhist <- history
                    src <- round(click$x)
                    if (src %in% c(history[,2],history[,3]))
                        {
                            history <- history[history[,2]!=src,]
                            history <- history[c(1:dim(history)[1],1),]
                            row <- dim(history)[1]
                            history[row,] <- c(abs(round(dblclick$y)),
                                               abs(src),
                                               abs(round(dblclick$x)),
                                               1,
                                               1,
                                               0,
                                               0)
                            ##            history
                        }
                    ##check and make sure that the oldest sink does not sink into another deme
                    history <- history[order(-history[,1]),]
                    if (history[1,3] %in% history[,2]) {history <- history[-1,]}
                }

if (debug())              print(paste("test of is.history",is.history(history)))
            if (is.history(history)) history else oldhist
        } else {NULL}
}


all.coalesce <- function(history)
{
    ac = TRUE
    pops <- unique(unlist(c(history[,2:3])))
    if (length(pops)!=(max(pops)+1))
    {
        ac=FALSE
    }
    history <- history[order(history[,1]),]
    popcoal <- rep(FALSE,length(pops))
    for (i in 1:dim(history)[1])
    {
        if (history[i,2]!=history[i,3]) #migration event occured
        {
            popcoal[pops==history[i,2]] <- TRUE
        }
    }
    if (sum(!popcoal)>0)
    {
        if (pops[!popcoal][1]%in%history[,3])
            popcoal[!popcoal] <- TRUE
    }
    
    if ((ac) & (sum(!popcoal)==0)) 
        ac=TRUE
    else
        ac=FALSE
#    print(ac)
    ac
}

is.history.new <- function(hist,ps,gr,mmats)
{
    err <- fsc.histEvCheck(hist,ps,gr,length(mmats))
    if (err)
    {
        if (sum(!(unique(hist[,7])%in%((1:length(mmats))-1)))>0) #asking for matrices that dont exist 
        {
            err <- FALSE
        }
    }
    
    err
}

is.history <- function(history)
    {
        
        err <- F
        msg <- NULL
        ##check for names
#        if (prod(names(history) %in% c("time", "source", "sink", "migrants", "new.size", "growth.rate", "migr.matrix"))!=1) err <- T
#        if (prod(c("time", "source", "sink", "migrants", "new.size", "growth.rate", "migr.matrix")%in%names(history))!=1) err <- T
#        if (err) msg <- c(msg,"problem with names")

        ######## this code is essentially worthless, but the problem
        ######## needs solving
        ##check to make sure that everybody coalesces
        ##first run through the sources
        last.sink <- history[history[,1]==max(history[,1]),3]
        if (length(unique(last.sink))>1) #multiple sinks at the "root"
            {
                err <- T
                msg <- c(msg,"multiple sinks at deepest time")
            }
        last.sink <- last.sink[1]
        
#        history$coal <- FALSE
        for  (s in 1:dim(history)[1])
            {
                snk <- history[s,3]
            }
        
        if (err)
        {
            print(paste("in simcoal is.history, value of err:",err) )
            print(msg)
            print(history)
            
            FALSE
        } else {
            TRUE
        }
    }


historiesEqual <- function(h1,h2)
{
    eq <- T
    if (dim(h1)[1]!=dim(h2)[1]) eq <- F
    if ((eq))
    {
        h1 <- h1[order(h1[,1]),]
        h2 <- h2[order(h2[,1]),]
        eq <- identical(h1,h2)
    }
if (debug())      print(paste("histories equal?",eq))
    eq
}
