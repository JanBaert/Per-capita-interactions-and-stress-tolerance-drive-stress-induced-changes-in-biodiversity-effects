# set working directory
setwd("/Users/Jan/Documents/UGent/Papers/Per capita interactions and stress tolerance drive stress-induced changes in biodiversity effects on ecosystem functions/Supplementary information")

require(nlme)
    
quartz(width=6.5,height=3.75)
layout(matrix(c(1,2,3,3),nrow=2,byrow=T),height=c(.9,.1))
par(oma=c(0,1,1,1))
load("change_effects.RData")
load("data_traits.RData")
par(mar=c(4.5,3.5,2.75,1.75))
data <- change_effects
plot(0,0,xlim=c(-0.02,0.02),ylim=c(-0.02,0.02),col="white",main="",xlab="",ylab="",cex.lab=0.95,cex.axis=1.05,bty="n")
axis(side=2,c(-10,10),c("",""))
axis(side=1,c(-10,10),c("",""))
abline(0,1,lty=2,col="grey10",lwd=1.5)
model <- lme(betadom~log10(diversity)*day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              
points(data$betadom,predict(model),pch=21,col="grey10",bg="grey50",cex=1.1)
model <- lme(betatdc~+day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              
points(data$betatdc,predict(model),pch=21,col="grey10",bg="darkorange3",cex=1.1)
model <- lme(betatic~1,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              
points(data$betatic,predict(model),pch=21,col="grey10",bg="dodgerblue4",cex=1.1)

plot(0,0,xlim=c(-0.02,0.02),ylim=c(-0.02,0.02),col="white",main="",xlab="",ylab="",cex.lab=0.95,cex.axis=1.05,bty="n")
axis(side=2,c(-10,10),c("",""))
axis(side=1,c(-10,10),c("",""))
abline(0,1,lty=2,col="grey10",lwd=1.5)
model  <- lme(betadom~sens_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              
points(data_traits$betadom[which(data_traits$day>14)],predict(model),pch=21,col="grey10",bg="grey50",cex=1.1)
model  <- lme(betatdc~sens_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              
points(data_traits$betatdc[which(data_traits$day>14)],predict(model),pch=21,col="grey10",bg="darkorange3",cex=1.1)
model  <- lme(betatic~int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              
points(data_traits$betatic[which(data_traits$day>14)],predict(model),pch=21,col="grey10",bg="dodgerblue4",cex=1.1)
mtext(expression(paste("Observed ",paste(Delta,"biodiversity effect"))),side=1,outer=T,line=-3.8,cex=1)
mtext(expression(paste("Predicted ",paste(Delta,"biodiversity effect"))),side=2,outer=T,line=-.85,at=.57,cex=1)
mtext("A: Model 1",3,at=.275,outer=T,line=-1.25,cex=1.2)
mtext("B: Model 2",3,at=.775,outer=T,line=-1.25,cex=1.2)

par(mar=c(0,3.2,0,1))
plot.new()
legend(0,1,c("Dominance"),pch=16,col=c("grey30"),bty="n",cex=1,horiz=T)
legend(.22,1,c("Trait-dep. complementarity","Trait-indep. complementarity"),pch=21,col="grey10",pt.bg=c("darkorange3","dodgerblue4"),bty="n",cex=1,horiz=T)
mtext(expression(paste("R"^2,"=0.052")),3,at=.16,outer=T,line=-4.1,cex=.85)
mtext(expression(paste("R"^2,"=0.459")),3,at=.66,outer=T,line=-4.1,cex=.85)
