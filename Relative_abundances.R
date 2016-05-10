# set working directory
setwd("/Users/Jan/Documents/UGent/Data/Experiment 1/Supporting information")

# load data
load("data_biov.RData")

# calculate relative abundances
data   <- data_biov[-which(data_biov$diversity==1),]
comp   <- unique(data$composition)

rel_ab <- NULL

for(i in 1:length(comp)){
  
      for(t in c(7,14,21,28)){
          
          for(rep in 1:3){
              
                temp    <- data[which(data$composition==comp[i]&data$day==t&data$replicate==rep),]
                if(nrow(temp)==3){
                      rel_ab  <- rbind(rel_ab,c(comp[i],t,rep,temp[1,13:20]/temp[1,21],temp[2,13:20]/temp[2,21],temp[3,13:20]/temp[3,21]))}}}}


colours <- c("black","blue","red","green","yellow","orange","cyan","purple")

quartz(width=12,height=7)
layout(matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow=3,byrow=T),height=c(.45,.45,.1))
par(oma=c(1,3,1,1))
for(t in c(7,14,21,28)){
  temp <- rel_ab[which(rel_ab[,2]==t),]
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="",ylab="",main=paste("Day ",t),cex.axis=1.2,cex.main=1.2)
  sapply(1:8,function(x) points(temp[,3+x],temp[,11+x],col=colours[x],pch=16))
  abline(0,1,lty=2)}
mtext(expression(paste("Relative abundance 25",paste(mu,"g L"^-1))),2,outer=T,cex=.85,at=.775)

for(t in c(7,14,21,28)){
  temp <- rel_ab[which(rel_ab[,2]==t),]
  plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="",ylab="",main=paste("Day ",t),cex.axis=1.2,cex.main=1.2)
  sapply(1:8,function(x) points(temp[,3+x],temp[,19+x],col=colours[x],pch=16))
  abline(0,1,lty=2)}
mtext(expression(paste("Relative abundance 250",paste(mu,"g L"^-1))),2,outer=T,cex=.85,at=.325)
mtext(expression(paste("Relative abundance 0",paste(mu,"g L"^-1))),1,outer=T,cex=.9,at=.5,line=-6)
mtext(expression(paste("Relative abundance 0",paste(mu,"g L"^-1))),1,outer=T,cex=.9,at=.5,line=-29)


par(mar=rep(0,4))
plot.new()  
legend("bottom",c("Coscinodiscus sp.","Ditylum sp.","Bacillaria sp.","Odontella sp.","Thalassiosira sp.1","Gyrosigma sp.","Guinardia sp.","Thalassiosira sp.2"),col=colours,cex=1.2,pch=16,bty="n",horiz=T)
