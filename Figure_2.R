# set working directory
setwd("/Users/Jan/Documents/UGent/Papers/Per capita interactions and stress tolerance drive stress-induced changes in biodiversity effects on ecosystem functions/Supplementary information")

# Figure 2
load("change_effects.RData")
biodiversityeffects <- c("A: Dominance ","B: Trait-dep. complementarity ","C: Trait-indep. complementarity ")

div_range <- c(2,4,6,8)
quartz(width=8,height=3.5)
layout(matrix(c(1:3),nrow=1),width=c(0.4,0.3,0.3))
par(oma=c(5,1,1,1))
colours <- c("grey50","darkorange2","firebrick2","dodgerblue3")
points  <- rep(16,4)
lim_plot <- c(-0.02,0.02)
time     <- c(7,14,21,28)
pch      <- c(21,21,21,21)
for(i in c(1:length(biodiversityeffects)))
        { 
                  if(i==1){par(mai=c(0.25,0.8,0.5,0.05))}else{par(mai=c(0.25,0,0.5,0.05))}
                  plot(1,1,col="white",xlab="",ylab="",main="",xlim=c(0.25,1),ylim=c(-0.02,0.02),xaxt="n",yaxt="n",bty="n")
                  rect(-10,-10,45,20,col="white",border=NA)
                  rect(-10,-20,45,0,col="grey90",border=NA)
                  axis(1,lab=c("","0","0.2","0.4","0.6","0.8","1",""),at=c(-10,0,0.2,0.4,0.6,0.8,1,10),cex.axis=1.4)
                  if(i==1){
                            axis(side=2,at=c(-10,-0.02,-0.01,0,0.01,0.02,10),labels=c("",round(lim_plot[1],digits=2),round(lim_plot[1]/2,digits=2),0,round(lim_plot[2]/2,digits=2),round(lim_plot[2],digits=2),""),cex.axis=1.25)
                            mtext(expression(paste(Delta,paste("biodiveristy effect (",paste(paste("mm"^3,paste(mu,"g"^-1)),")")))),side=2,line=3,cex=.95)
                            legend(x=0.22,y=0.024,c("day 7","day 14","day 21","day 28"),col="grey10",pt.bg=colours,pch=pch,bty='n',cex=1.2,ncol=2)
                  }else{    axis(side=2,at=c(-100,100),c("",""))}
                          
                  mtext(biodiversityeffects[i],side=3,line=2,cex=.9)
                  for(t in c(7,14,21,28))
                            { for(div in div_range)
                                { temp <- change_effects[which(change_effects$day==t&change_effects$diversity==div),4+(i-1)*3]
                                  points(rep(log10(div)+(which(time==t)-2)*0.0175,length(temp)),temp,pch=pch[which(time==t)],col="grey10",bg=colours[which(time==t)],cex=1.22)}
                                  if(i==1){segments(log10(1.8),(0.001425687-0.000173067*t)+(-0.003680934+0.000481080*t)*log10(1.8),log10(9.5),(0.001425687-0.000173067*t)+(-0.003680934+0.000481080*t)*log10(9.5),col=colours[which(time==t)],lty=2,lwd=1.35)}
                                  if(i==2){segments(log10(1.8),(-0.00083902+0.0000945862*t),log10(9.5),(-0.00083902+0.0000945862*t),col=colours[which(time==t)],lty=2,lwd=1.35)}  
                                  if(i==3){segments(log10(1.8),-0.0012799,log10(9.5),-0.0012799,col=colours[which(time==t)],lty=2,lwd=1.35)}
                  }}

mtext("Log species richness",side=1,line=1.6,outer=T,at=0.55)


