# set working directory
setwd("/Users/Jan/Documents/UGent/Data/Experiment 1/Supporting information")

# load packages and data
require(deSolve)     
load("data_biov.RData")    

# load vector with the order in which parameters were estimated
load("compositions_list.RData")

# load model precitions
load("MC_out.RData")
load("M_out.RData")


# calculate predicted and observed relative species yields
time <- c(7,14,21,28)
Y_pred_obs     <-  array(dim=c(length(compositions)*8*2*4,4,5))
RY_pred_obs    <-  array(dim=c(length(compositions)*8*2*4,4,5))

for(scen in 1:5){
    for(comp in compositions[1:31]){
      spec            <- strsplit(toString(comp),"/",fixed=TRUE)[[1]]
      n_spec          <- length(spec)
      
      # obeserved yields
      data_0          <- data_biov[which(data_biov$composition==comp&data_biov$concentration==0),c(2,13:20)]
      ind             <- which(data_0[1,2:9]!=0)
      data_0          <- data_0[,c(1,1+which(data_0[1,c(2:9)]!=0))]
      data_0          <- t(sapply(time,function(x) apply(data_0[which(data_0[,1]==x),2:(n_spec+1)],2,mean)))
      data_250        <- data_biov[which(data_biov$composition==comp&data_biov$concentration==250),c(2,13:20)]
      data_250        <- data_250[,c(1,1+which(data_250[1,c(2:9)]!=0))]
      data_250        <- t(sapply(time,function(x) apply(data_250[which(data_250[,1]==x),2:(n_spec+1)],2,mean)))
      Y_obs_0         <- data_0
      Y_obs_250       <- data_250
      RY_obs_0        <- t(sapply(1:nrow(data_0),function(x) data_0[x,]/sum(data_0[x,])))
      RY_obs_250      <- t(sapply(1:nrow(data_250),function(x) data_250[x,]/sum(data_250[x,])))
      
      #predicted realtive yields
      pred_0          <-  MC_out[which(compositions==comp),ind,,which(c(0,25,250)==0),,scen]
      pred_0          <- t(apply(pred_0,c(1,2),median))
      pred_250        <- MC_out[which(compositions==comp),ind,,which(c(0,25,250)==250),,scen]
      pred_250        <- t(apply(pred_250,c(1,2),median))
      Y_pred_0        <- t(sapply(1:length(time),function(x) pred_0[x,]))
      Y_pred_250      <- t(sapply(1:length(time),function(x) pred_250[x,]))
      RY_pred_0       <- t(sapply(1:length(time),function(x) pred_0[x,]/sum(pred_0[x,])))
      RY_pred_250     <- t(sapply(1:length(time),function(x) pred_250[x,]/sum(pred_250[x,])))
      BC_0            <- sapply(1:nrow(Y_pred_0),function(x) sum(abs(Y_obs_0[x,]-Y_pred_0[x,]))/(sum(Y_obs_0[x,]+Y_pred_0[x,])))
      BC_250          <- sapply(1:nrow(Y_pred_0),function(x) sum(abs(Y_obs_250[x,]-Y_pred_250[x,]))/(sum(Y_obs_250[x,]+Y_pred_250[x,])))
      BC              <- rbind(sum(BC_0)/4,sum(BC_250)/4)
      pred_0          <- rowSums(pred_0)
      pred_250        <- rowSums(pred_250)
      resp_pred       <- ((pred_0-pred_250)/pred_0)
      Y_pred_obs[(which(is.na(Y_pred_obs[,1,scen]))[1]+c(0:(2*length(spec)*length(time)-1))),,scen] <- cbind(as.vector(c(Y_obs_0,Y_obs_250)),as.vector(c(Y_pred_0,Y_pred_250)),rep(length(spec),2*length(spec)*length(time)),c(rep(0,length(spec)*length(time)),rep(250,length(spec)*length(time))))
      RY_pred_obs[(which(is.na(RY_pred_obs[,1,scen]))[1]+c(0:(2*length(spec)*length(time)-1))),,scen] <- cbind(as.vector(c(RY_obs_0,RY_obs_250)),as.vector(c(RY_pred_0,RY_pred_250)),rep(length(spec),2*length(spec)*length(time)),c(rep(0,length(spec)*length(time)),rep(250,length(spec)*length(time))))
      }}

# Supplementray Figure 2
colours=c("black","grey")
quartz(width=13,height=6)
layout(matrix(c(1:10,11,11,11,11,11),byrow=T,nrow=3),height=c(.45,.45,.06))
par(oma=rep(1,4),mai=c(.55,.75,.5,.0))
for(i in 1:5){
  plot(0,0,xlim=c(1,10),ylim=c(1,10),xlab=expression(paste("Observed density (mm"^3,paste(" L"^-1,")"))),ylab=expression(paste("Predicted density (mm"^3,paste(" L"^-1,")"))),cex.lab=1.1,cex.axis=1.2,main=paste("Scenario ",i),cex.main=1.5)
  for(c in c(0,250)){
  points(log10(Y_pred_obs[which(Y_pred_obs[,4,i]==c),1,i]),log10(Y_pred_obs[which(Y_pred_obs[,4,i]==c),2,i]),col=colours[which(c(0,250)==c)],cex=.3)}
  abline(0,1,lty=2)}
par(mai=c(.55,.75,.5,.0))
for(i in 1:5){
  plot(0,0,xlim=c(0,1),ylim=c(0,1),main="",xlab="Observed relative abundance (-)",ylab="Predicted relative abundance (-)",cex.lab=1.1,cex.axis=1.2)
  for(c in c(0,250)){
    points(RY_pred_obs[which(RY_pred_obs[,4,i]==c),1,i],RY_pred_obs[which(RY_pred_obs[,4,i]==c),2,i],col=colours[which(c(0,250)==c)],cex=.3)}
  abline(0,1,lty=2)}
par(mai=rep(0,4))
plot.new()
legend("bottom",c(expression(paste("0 ",paste(mu,"g L"^-1))),expression(paste("250 ",paste(mu,"g L"^-1)))),pch=16,col=colours,horiz=T,bty="n",cex=1.2)

