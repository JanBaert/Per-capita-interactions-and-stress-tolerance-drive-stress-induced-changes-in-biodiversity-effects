# set working directory
setwd("/Users/Jan/Documents/UGent/Data/Experiment 1/Supporting information")

# MCMC predictions, data and compositions order
load("MC_out.RData")
load("M_out.RData")
load("data_biov.RData")
load("change_effects.RData")
load("compositions_list.RData")

# simulation settings
n_scen  <- 5
n_it    <- 1000
time  <- c(7,14,21,28)


# Calculate likelihoods and Bray-Curtis similarity
    lik_out_resp <- array(dim=c(n_scen,n_it))
    BC_sim <- array(dim=c(n_scen,n_it))
    
    for(scen in 1:n_scen){
      for(it in 1:n_it){
        BC_all         <- NULL
        res_resp_all   <- NULL
        for(comp in compositions){
          spec            <- strsplit(toString(comp),"/",fixed=TRUE)[[1]]
          n_spec          <- length(spec)
          
          # observed data
          data_0          <- data_biov[which(data_biov$composition==comp&data_biov$concentration==0),c(2,13:20)]
          ind             <- which(data_0[1,2:9]!=0)
          data_0          <- data_0[,c(1,1+which(data_0[1,c(2:9)]!=0))]
          data_0          <- t(sapply(time,function(x) apply(data_0[which(data_0[,1]==x),2:(n_spec+1)],2,median)))
          data_250        <- data_biov[which(data_biov$composition==comp&data_biov$concentration==250),c(2,13:20)]
          data_250        <- data_250[,c(1,1+which(data_250[1,c(2:9)]!=0))]
          data_250        <- t(sapply(time,function(x) apply(data_250[which(data_250[,1]==x),2:(n_spec+1)],2,median)))
          Y_obs_0         <- t(sapply(1:nrow(data_0),function(x) data_0[x,]))
          Y_obs_250       <- t(sapply(1:nrow(data_250),function(x) data_250[x,]))
          RY_obs_0        <- t(sapply(1:nrow(data_0),function(x) data_0[x,]/sum(data_0[x,])))
          RY_obs_250      <- t(sapply(1:nrow(data_250),function(x) data_250[x,]/sum(data_250[x,])))
          data_0          <- rowSums(data_0)
          data_250        <- rowSums(data_250)    
          resp_obs        <- (data_0-data_250)/data_0
          
          # predicted data
          pred_0          <- t(MC_out[which(compositions==comp),ind,,which(c(0,25,250)==0),it,scen])
          pred_250        <- t(MC_out[which(compositions==comp),ind,,which(c(0,25,250)==250),it,scen])
          Y_pred_0        <- t(sapply(1:length(time),function(x) pred_0[x,]))
          Y_pred_250      <- t(sapply(1:length(time),function(x) pred_250[x,]))
          RY_pred_0       <- t(sapply(1:length(time),function(x) pred_0[x,]/sum(pred_0[x,])))
          RY_pred_250     <- t(sapply(1:length(time),function(x) pred_250[x,]/sum(pred_250[x,])))
          BC_0            <- sapply(1:nrow(RY_pred_0),function(x) sum(abs(Y_obs_0[x,]-Y_pred_0[x,]))/(sum(Y_obs_0[x,]+Y_pred_0[x,])))
          BC_250          <- sapply(1:nrow(RY_pred_0),function(x) sum(abs(Y_obs_250[x,]-Y_pred_250[x,]))/(sum(Y_obs_250[x,]+Y_pred_250[x,])))
          BC              <- rbind(sum(BC_0),sum(BC_250))
          pred_0          <- rowSums(pred_0)
          pred_250        <- rowSums(pred_250)
          resp_pred       <- ((pred_0-pred_250)/pred_0)
    
          # calculate residuals
          res_resp        <- resp_pred-resp_obs
    
          # store residuals and BC similarities
          res_resp_all    <- c(res_resp_all,res_resp)
          BC_all          <- rbind(BC_all,BC)}
        
        # remove any NAs
        res_resp_all    <- res_resp_all[-which(is.na(res_resp_all))]
        BC_all          <- BC_all[-which(is.na(BC_all[,1])),]
        
        # calculate log likelihood
        log_lik_resp    <- sum(res_resp_all^2/var(res_resp_all))
        
        # store log likelihood and Bray-Curtis similarity
        lik_out_resp[scen,it] <- log_lik_resp
        BC_sim[scen,it]      <- sum(BC_all)}}


# Figure 3
    quartz(width=8,height=4)
    par(oma=rep(0.5,4),mfrow=c(1,2),mai=c(.8,0.8,0.5,0.3))
    
    # Log Likelihood
    boxplot(-lik_out_resp[1,],-lik_out_resp[2,],-lik_out_resp[3,],-lik_out_resp[4,],-lik_out_resp[5,],outline=F,xlab="",ylab="-Log likelihood",ylim=c(-150,-116),cex.lab=1.2)
    axis(1,at=c(-1,1,2,3,4,5,10),labels=c("","1","2","3","4","5",""))
    text(c(1:5),rep(-116.5,5),c("a","b","c","d","d"))
    mtext("Scenario",1,cex=1.2,line=2.25)
    mtext("A: Function",3,line=0.6,at=3,cex=1.45)
    
    # Bray-Curtis similarity
    boxplot(-BC_sim[1,]/(4*29),-BC_sim[2,]/(4*29),-BC_sim[3,]/(4*29),-BC_sim[4,]/(4*29),-BC_sim[5,]/(4*29),outline=F,xlab="",ylab="- Bray-Curtis index",ylim=c(-1,-0.58),cex.lab=1.2)
    axis(1,at=c(-1,1,2,3,4,5,10),labels=c("","1","2","3","4","5",""))
    text(c(1:5),rep(-0.59,5),c("a","b","c","d","e"))
    mtext("B: Composition",3,line=0.6,at=3,cex=1.45)
    mtext("Scenario",1,cex=1.2,line=2.25)

# test significance

  # Likelihood
  p12 <- wilcox.test(lik_out_resp[1,],lik_out_resp[2,])$p.value
  p13 <- wilcox.test(lik_out_resp[1,],lik_out_resp[3,])$p.value
  p14 <- wilcox.test(lik_out_resp[1,],lik_out_resp[4,])$p.value
  p15 <- wilcox.test(lik_out_resp[1,],lik_out_resp[5,])$p.value
  p23 <- wilcox.test(lik_out_resp[2,],lik_out_resp[3,])$p.value
  p24 <- wilcox.test(lik_out_resp[2,],lik_out_resp[4,])$p.value
  p25 <- wilcox.test(lik_out_resp[2,],lik_out_resp[5,])$p.value
  p34 <- wilcox.test(lik_out_resp[3,],lik_out_resp[4,])$p.value
  p35 <- wilcox.test(lik_out_resp[3,],lik_out_resp[5,])$p.value
  p45 <- wilcox.test(lik_out_resp[4,],lik_out_resp[5,])$p.value
  
  p.adjust(c(p12,p13,p14,p15,p23,p24,p25,p34,p35,p45),method="bonferroni")

  # Bray-Curtis similarity
  p12 <- wilcox.test(BC_sim[1,],BC_sim[2,])$p.value
  p13 <- wilcox.test(BC_sim[1,],BC_sim[3,])$p.value
  p14 <- wilcox.test(BC_sim[1,],BC_sim[4,])$p.value
  p15 <- wilcox.test(BC_sim[1,],BC_sim[5,])$p.value
  p23 <- wilcox.test(BC_sim[2,],BC_sim[3,])$p.value
  p24 <- wilcox.test(BC_sim[2,],BC_sim[4,])$p.value
  p25 <- wilcox.test(BC_sim[2,],BC_sim[5,])$p.value
  p34 <- wilcox.test(BC_sim[3,],BC_sim[4,])$p.value
  p35 <- wilcox.test(BC_sim[3,],BC_sim[5,])$p.value
  p45 <- wilcox.test(BC_sim[4,],BC_sim[5,])$p.value
  
  p.adjust(c(p12,p13,p14,p15,p23,p24,p25,p34,p35,p45),method="bonferroni")

