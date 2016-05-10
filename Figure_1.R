# set working directory
setwd("/Users/Jan/Documents/UGent/Data/Experiment 1/Supporting information")

# Figure 1
    load("data_biov.RData")
    letter <- c("A","B","C","D")    

    time    <- c(7,14,21,28)
    conc    <- c(0,25,250)
    div     <- c(1,2,4,6,8)
    col     <- c("salmon4","darkorange2","grey25")

    # open window and define graphical parameters
    quartz(width=10,height=3.5)
    layout(matrix(c(1:4),byrow=T,nrow=1),heights=c(0.5,0.5),width=c(0.31,0.23,0.23,0.23))
    par(oma=c(5,1,1,1),mai=c(0.25,0,0.5,0.05),mar=c(0.1,1,1,1))
    
    for(t in time){ temp <- data_biov[which(data_biov$day==t),]
                    for(c in conc){ temp2 <- temp[which(temp$concentration==c),]
        
                    if(t==7){ par(mai=c(0.25,0.8,0.5,0.05))
                              if(c==0){plot(0,0,col="white",xlim=c(-0.1,1),ylim=c(6,10.6),ylab=expression(paste("Log biovolume (",paste(paste("mm"^3," L"^-1),")"))),xaxt="n",yaxt="n",xlab="",main="",cex.lab=1.5,cex.main=1.3,cex.axis=1.4)
                                       rect(-1,0,2,7,border="white",col="grey90")
                                       rect(-1,7,2,12,border="white",col="grey99")
                                       legend("topleft",c(expression(paste("0",paste(mu,"g L"^-1))),expression(paste("25",paste(mu,"g L"^-1))),expression(paste("250",paste(mu,"g L"^-1)))),pch=16,col=col,bty="n",horiz=F,cex=1.2)
                                       axis(1,lab=c("","0","0.2","0.4","0.6","0.8","1",""),at=c(-10,0,0.2,0.4,0.6,0.8,1,10),cex.axis=1.4)
                                       axis(2,lab=c("","","","","","",""),at=c(0,6,7,8,9,10,12),hadj=-.2)
                                       mtext("A: Day 7",side=3,line=1,cex=1.2)
                                       mtext("6",1,at=-0.24,cex=.9,line=-.9)
                                       mtext("7",1,at=-0.24,cex=.9,line=-3.8)
                                       mtext("8",1,at=-0.24,cex=.9,line=-6.8)
                                       mtext("9",1,at=-0.24,cex=.9,line=-9.8)
                                       mtext("10",1,at=-0.24,cex=.9,line=-12.83)}
                    }else{
                      par(mai=c(0.25,0,0.5,0.05))
                      if(c==0){plot(0,0,col="white",xlim=c(-0.1,1),ylim=c(6,10.6),yaxt="n",ylab="",xlab="",main="",cex.lab=1.4,cex.main=1.3,cex.axis=1.4,xaxt="n")
                               rect(-1,0,2,7,border="white",col="grey90")
                               rect(-1,7,2,12,border="white",col="grey99")
                               axis(1,lab=c("","0","0.2","0.4","0.6","0.8","1",""),at=c(-10,0,0.2,0.4,0.6,0.8,1,10),cex.axis=1.4)
                               mtext(paste(paste(letter[which(time==t)],": Day ",sep=""),t),side=3,line=1,cex=1.2)}}
        
        sapply(div,function(x) points(rep(log10(x)+(which(conc==c)-2)*0.025,length(which(temp2$diversity==x))),log10(temp2$total[which(temp2$diversity==x)]),pch=16,col=col[which(conc==c)],cex=1.2))
        if(c==0){abline(7.798528+0.041006*t,0.475281-0.046298*t,col=col[which(conc==c)],lty=2)}
        if(c==25){abline(7.798528-0.016252+0.041006*t+0.001073*t,0.475281+0.140380-0.046298*t-0.007792*t,col=col[which(conc==c)],lty=2)}
        if(c==250){abline(7.798528-0.307387+0.041006*t-0.043171*t,0.475281-0.785644-0.046298*t+0.097789*t,col=col[which(conc==c)],lty=2)}    }}
mtext("Log species richness",side=1,line=2,cex=1.2,outer=T,at=0.55)

