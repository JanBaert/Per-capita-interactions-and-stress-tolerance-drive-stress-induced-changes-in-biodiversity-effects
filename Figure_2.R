# set working directory
setwd("/Users/Jan/Documents/UGent/Data/Experiment 1/Supporting information")

# Figure 2
load("change_effects.RData")
biodiversityeffects <- c("A: Dominance ","B: Trait-dep. complementarity ","C: Trait-indep. complementarity ")

div_range <- c(2,4,6,8)
quartz(width=8,height=3.5)
layout(matrix(c(1:3),nrow=1),width=c(0.4,0.3,0.3))
par(oma=c(5,1,1,1))
colours <- c("salmon4","darkorange2","grey25","lightsteelblue4")
points  <- rep(16,4)
lim_plot <- c(-0.02,0.02)
for(i in c(1:length(biodiversityeffects)))
        { 
                  if(i==1){par(mai=c(0.25,0.8,0.5,0.05))}else{par(mai=c(0.25,0,0.5,0.05))}
                  plot(1,1,col="white",xlab="",ylab="",main="",xlim=c(0.25,1),ylim=c(-0.02,0.02),xaxt="n",yaxt="n")
                  rect(-10,-10,45,20,col="grey99",border=NA)
                  rect(-10,-20,45,0,col="grey90",border=NA)
                  axis(1,lab=c("","0","0.2","0.4","0.6","0.8","1",""),at=c(-10,0,0.2,0.4,0.6,0.8,1,10),cex.axis=1.4)
                  if(i==1){
                            axis(side=2,at=c(-10,-0.02,-0.01,0,0.01,0.02,10),labels=c("",round(lim_plot[1],digits=2),round(lim_plot[1]/2,digits=2),0,round(lim_plot[2]/2,digits=2),round(lim_plot[2],digits=2),""),cex.axis=1.25)
                            mtext(expression(paste(Delta,paste("biodiveristy effect (",paste(paste("mm"^3,paste(mu,"g"^-1)),")")))),side=2,line=2.5,cex=0.95)
                            legend('topleft',c("day 7","day 14","day 21","day 28"),col=colours,pch=16,bty='n',cex=1.1,ncol=2)}
                          
                  mtext(biodiversityeffects[i],side=3,line=1,cex=.9)
                  for(t in c(7,14,21,28))
                            { for(div in div_range)
                                { temp <- change_effects[which(change_effects$day==t&change_effects$diversity==div),4+(i-1)*3]
                                  points(rep(log10(div)+(which(time==t)-2)*0.015,length(temp)),temp,pch=points[which(time==t)],col=colours[which(time==t)],cex=1.2)}
                                if(i==1){abline(0.001425687-0.000173067*t,-0.003680934+0.000481080*t,col=colours[which(time==t)],lty=2)}
                                if(i==2){ abline(-0.00083902+0.0000945862*t,0,col=colours[which(time==t)],lty=2)}
                                if(i==3){abline(-0.0012799,0,col=colours[which(time==t)],lty=2)}}}

mtext("Log species richness",side=1,line=1.6,outer=T,at=0.55)


