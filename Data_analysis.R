#  custom functions 


  # convert factor variable to numeric variable
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}                      
  
  # plot correlations in upper triangle of pair-plot
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){                  
                                                                    usr <- par("usr"); on.exit(par(usr))
                                                                    par(usr = c(0, 1, 0, 1))
                                                                    r <- abs(cor(x, y))
                                                                    txt <- format(c(r, 0.123456789), digits=digits)[1]
                                                                    txt <- paste(prefix, txt, sep="")
                                                                    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
                                                                    text(0.5, 0.5, txt, cex = 2)}


#######################################################################
#                                                                     #
#                      Trait-independent analysis                     #
#                                                                     #
#######################################################################

# set working directory
setwd("/Users/Jan/Documents/UGent/Papers/Per capita interactions and stress tolerance drive stress-induced changes in biodiversity effects on ecosystem functions/Supplementary information")
  
# calculate biodiversity effects
source("Calculate_biodiversity_effects.R")
  
# load nlme package
require(nlme)


data <- change_effects

  # 1. Dominance effect
  
        # initial model
        model <- lme(betadom~log10(diversity)*day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # optimal model
        model <- lme(betadom~log10(diversity)*day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # print residual diagnostics
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$diversity,ylab="Model residuals",xlab="Diversity",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$day,ylab="Model residuals",xlab="Day (d)",cex.lab=1.2,cex.axis=1.2)
        plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta,"dominace effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data$betadom,predict(model),ylab=expression(paste("Predicted ",paste(Delta,"dominace effect"))),xlab=expression(paste("Oserved ",paste(Delta,"dominace effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-25,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-25,cex=1.25)
        
  
        quartz(width=9,height=6)
        par(oma=c(11,11,1,1))
        resmat <- cbind(residuals(model)[which(data$day==7)],residuals(model)[which(data$day==14)],residuals(model)[which(data$day==21)],residuals(model)[which(data$day==28)])
        pairs(resmat,labels=c("Day 7","Day 14","Day 21","Day 28"),upper.panel=panel.cor,cex.lab=1.2,cex.axis=1.2)
        mtext("Model residuals (-)",side=1,line=9.5,at=.36,outer=T)
        mtext("Model residuals (-)",side=2,line=9.5,at=.30,outer=T)
  
  # 2. Trait-dependent complementarity effect
  
        # initial model 
        model <- lme(betatdc~log10(diversity)*day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # optimal model
        model <- lme(betatdc~+day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # print residual diagnostics
        summary(model)
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$diversity,ylab="Model residuals",xlab="Diversity",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$day,ylab="Model residuals",xlab="Day (d)",cex.lab=1.2,cex.axis=1.2)
        plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta," trait-dep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data$betatdc,predict(model),ylab=expression(paste("Predicted ",paste(Delta,"trait-dep. comp. effect"))),xlab=expression(paste("Observed ",paste(Delta,"trait-dep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-24,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-24,cex=1.25)
        
        quartz(width=9,height=6)
        par(oma=c(11,11,1,1))
        resmat <- cbind(residuals(model)[which(data$day==7)],residuals(model)[which(data$day==14)],residuals(model)[which(data$day==21)],residuals(model)[which(data$day==28)])
        pairs(resmat,labels=c("Day 7","Day 14","Day 21","Day 28"),upper.panel=panel.cor,cex.lab=1.2,cex.axis=1.2)
        mtext("Model residuals (-)",side=1,line=9.5,at=.36,outer=T)
        mtext("Model residuals (-)",side=2,line=9.5,at=.30,outer=T)
  
  # 3. Trait-independent complementarity effect  
  
        # initial model
        model <- lme(betatic~log10(diversity)*day,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # optimal model
        model <- lme(betatic~1,correlation=corCAR1(form=~day|compositions),random=(~1|compositions),data=data)              # initial model
        summary(model)
        
        # print residual diagnostics
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$diversity,ylab="Model residuals",xlab="Diversity",cex.lab=1.2,cex.axis=1.2)
        boxplot(residuals(model)~data$day,ylab="Model residuals",xlab="Day (d)",cex.lab=1.2,cex.axis=1.2)
        plot(round(predict(model),4),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta,"trait-indep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data$betatic,round(predict(model),4),ylab=expression(paste("Predicted ",paste(Delta,"trait-indep. comp. effect"))),xlab=expression(paste("Observed ",paste(Delta,"trait-indep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-24,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-24,cex=1.25)
        
        quartz(width=9,height=6)
        par(oma=c(11,11,1,1))
        resmat <- cbind(residuals(model)[which(data$day==7)],residuals(model)[which(data$day==14)],residuals(model)[which(data$day==21)],residuals(model)[which(data$day==28)])
        pairs(resmat,labels=c("Day 7","Day 14","Day 21","Day 28"),upper.panel=panel.cor,cex.lab=1.2,cex.axis=1.2)
        mtext("Model residuals (-)",side=1,line=9.5,at=.36,outer=T)
        mtext("Model residuals (-)",side=2,line=9.5,at=.30,outer=T)
  
  
#######################################################################
#                                                                     #
#                         Trait-based analysis                        #
#                                                                     #
#######################################################################
  
# load trait-independent datasets
source("Calculate_traits.R")  
  
  # 1. Dominance effect
  
        # initial model
        model  <- lme(betadom~sens_w*int_w+sens_w*d_int_w+d_int_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # optimal model
        model  <- lme(betadom~sens_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # plot residual diagnostics
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        plot(data_traits$sens_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted species sensitivity",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted per capita interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$d_int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted change in  interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta,"dominace effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$betadom[which(data_traits$day>14)],predict(model),ylab=expression(paste("Predicted ",paste(Delta,"dominace effect"))),xlab=expression(paste("Observed ",paste(Delta,"dominace effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-25,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-25,cex=1.25)
        mtext("F",3,at=.67,outer=T,line=-25,cex=1.25)
        
  # 2. Trait-dependent complementarity effect
  
        # initial model
        model  <- lme(betatdc~sens_w*int_w+sens_w*d_int_w+d_int_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # optimal model
        model  <- lme(betatdc~sens_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # plot residual diagnostics
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        plot(data_traits$sens_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted species sensitivity",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted per capita interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$d_int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted change in  interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta,"trait-dep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$betatdc[which(data_traits$day>14)],predict(model),ylab=expression(paste("Predicted ",paste(Delta,"trait-dep. comp. effect"))),xlab=expression(paste("Observed ",paste(Delta,"trait-dep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-24,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-24,cex=1.25)
        mtext("F",3,at=.67,outer=T,line=-24,cex=1.25)
        
  # 3. Trait-dependent complementarity effect
  
        # initial model
        model  <- lme(betatic~sens_w*int_w+sens_w*d_int_w+d_int_w*int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # optimal model
        model  <- lme(betatic~int_w,random=(~1|comp),data=data_traits[which(data_traits$day>14),])              # initial model
        summary(model)
        
        # plot residual diagnostics
        quartz(width=9,height=6)
        par(mfrow=c(2,3),oma=rep(1,4))
        qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
        plot(data_traits$sens_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted species sensitivity",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted per capita interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$d_int_w[which(data_traits$day>14)],residuals(model),ylab="Model residuals",xlab="Weighted change in  interaction strength",cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste("Predicted ",paste(Delta,"trait-indep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(h=0,lty=2)
        plot(data_traits$betatic[which(data_traits$day>14)],predict(model),ylab=expression(paste("Predicted ",paste(Delta,"trait-indep. comp. effect"))),xlab=expression(paste("Observed ",paste(Delta,"trait-indep. comp. effect"))),cex.lab=1.2,cex.axis=1.2)
        abline(0,1,lty=2)
        mtext("A",3,at=.01,outer=T,line=-3,cex=1.25)
        mtext("B",3,at=.34,outer=T,line=-3,cex=1.25)
        mtext("C",3,at=.67,outer=T,line=-3,cex=1.25)
        mtext("D",3,at=.01,outer=T,line=-24,cex=1.25)
        mtext("E",3,at=.34,outer=T,line=-24,cex=1.25)
        mtext("F",3,at=.67,outer=T,line=-24,cex=1.25)
        
#######################################################################
#                                                                     #
#                           Log10 Biovolume                           #
#                                                                     #
#######################################################################
  
data  <- data_biov[-which(data_biov$day<7|data_biov$day>28|data_biov$day==9),]
  
  

  # optimize model
  model <- lme(log10(total)~as.factor(concentration)*day*log10(diversity),random=(~1|composition),data=data)        # initial model
  
  # test for random intercept
  model2 <- lm(log10(total)~as.factor(concentration)*day*log10(diversity),data=data)   
  anova(model,model2)
  
  summary(model)
  quartz(width=9,height=6)
  par(mfrow=c(2,3),oma=c(1,1,1,1),mar=c(5,5,3,2))
  qqnorm(residuals(model),main="",cex.lab=1.2,cex.axis=1.2)
  boxplot(residuals(model)~data$concentration,ylab="Model residuals",xlab=expression(paste("Concentration (",paste(mu,paste("g L"^-1,")")))),cex.lab=1.2,cex.axis=1.2)
  boxplot(residuals(model)~data$diversity,ylab="Model residuals",xlab="Diversity",cex.lab=1.2,cex.axis=1.2)
  boxplot(residuals(model)~data$day,ylab="Model residuals",xlab="Day (d)",cex.lab=1.2,cex.axis=1.2)
  plot(predict(model),residuals(model),ylab="Model residuals",xlab=expression(paste(paste("Predicted log"[10]," biovolume ("),paste(mu,paste("m"^3,")")))),cex.lab=1.2,cex.axis=1.2)
  abline(h=0,lty=2)
  plot(log10(data$total),predict(model),ylab=expression(paste(paste("Predicted log"[10]," biovolume ("),paste(mu,paste("m"^3,")")))),xlab=expression(paste(paste("Observed log"[10]," biovolume ("),paste(mu,paste("m"^3,")")))),cex.lab=1.15,cex.axis=1.2)
  abline(0,1,lty=2)
  mtext("A",3,at=.02,outer=T,line=-2,cex=1.25)
  mtext("B",3,at=.36,outer=T,line=-2,cex=1.25)
  mtext("C",3,at=.69,outer=T,line=-2,cex=1.25)
  mtext("D",3,at=.02,outer=T,line=-24,cex=1.25)
  mtext("E",3,at=.36,outer=T,line=-24,cex=1.25)
  mtext("F",3,at=.69,outer=T,line=-24,cex=1.25)
  
  quartz(width=9,height=6)
  par(oma=c(10,10,1,1))
  resmat <- cbind(residuals(model)[which(data$day==7)],residuals(model)[which(data$day==14)],residuals(model)[which(data$day==21)],residuals(model)[which(data$day==28)])
  pairs(resmat,labels=c("Day 7","Day 14","Day 21","Day 28"),upper.panel=panel.cor,cex.lab=1.2,cex.axis=1.2)
  mtext("Model residuals (-)",side=1,line=8.5,at=.39,outer=T)
  mtext("Model residuals (-)",side=2,line=8.5,at=.33,outer=T)
  
    