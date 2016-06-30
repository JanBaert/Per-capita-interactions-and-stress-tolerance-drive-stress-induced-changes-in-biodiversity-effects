# custom function to convert factor variable in a numeric variable
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# 1. Load data and create dataframes
      
      # load data
      data_dens         <- read.delim("data.txt", header = TRUE)
      data_dens         <- as.data.frame(data_dens)
      names(data_dens)  <- c("composition","day","replicate","concentration","8","10","17","18","E18","23","31","32","nitrate","phosphate","silicate","atrazine")

      # remove NA from dataframe
      data_dens         <- data_dens[which(is.na(data_dens[,5])==FALSE),]                                         # dataframe with cell densities (cells/ml)

      # convert cell densities to bio-volumes
      cellvolumes       <- c(367008,24757,6448,752767,25727,3940,14001,12556)                                     # cell volumes ?m?
      densities         <- as.matrix(data_dens[,5:12])
      biovolumes        <- as.matrix(data_dens[,5:12])%*%diag(cellvolumes)
      total_biovolume   <- as.matrix(data_dens[,5:12])%*%as.matrix(cellvolumes)
      temp <- as.matrix(biovolumes==0)
      temp[temp==TRUE]  <- 1
      diversity <- t(t(8-rowSums(temp)))

      # construct dataframe with biovolumes (cubic micrometer/ml)
      data_biov         <- data.frame(data_dens[,1:12],biovolumes,total_biovolume,diversity,data_dens[13:16])       # dataframe with volumetric densities (?m?/ml)
      names(data_biov)  <- c("composition","day","replicate","concentration","n8","n10","n17","n18","nE18","n23","n31","n32","8","10","17","18","E18","23","31","32","total","diversity","nitrate","phosphate","silicate","atrazine")

      # convert composition to character string
      for(i in 1:nrow(data_biov)){data_biov$composition[i] <- as.character(data_biov$composition[i])}

      # save data_frame
      save(data_biov,file="data_biov.RData")

# 2. Calculate biodiversity effects (Fox 2005)

      # load biovolumes
      load("data_biov.RData")

      # vector with reference numbers of species
      spec_range <- c("8","10","17","18","E18","23","31","32")

      # calculate monoculture yiels
      monoculture_yields <- NULL
      
      for(t in c(0,7,14,21,28)){   m_yield_temp     <- array(dim=c(9,11))
                                      m_yield_temp[,1] <- rep(t,9)
                                      m_yield_temp[,2] <- c(1,2,3,1,2,3,1,2,3)
                                      m_yield_temp[,3] <- c(0,0,0,25,25,25,250,250,250)
                                      for(i in c(1:length(spec_range))){ m_yield_temp[,3+i] <- data_biov[which(data_biov$composition==spec_range[i]&data_biov$day==t),12+i]}
                                      monoculture_yields <- rbind(monoculture_yields,m_yield_temp)}


      # calculate biodiveristy effects
      effects_analysis <- NULL

      for(i in c(2,4,6,8)){ temp         <- data_biov[which(data_biov$diversity==i),]
                            compositions <- unique(temp$composition)
                            
                            for(j in compositions){ temp2      <- temp[which(temp$composition==j),]
                                                    spec       <- strsplit(toString(j),"/",fixed=TRUE)[[1]]                              
                                                    spec_ind   <- which(spec_range%in%spec)
                            
                                                    for(t in c(7,14,21,28,49)){ temp3    <- temp2[which(temp2$day==t),c(3,4,13:21)]
           
                                                                                for(conc in unique(temp3$concentration)){ M             <- as.vector(apply(monoculture_yields[which(monoculture_yields[,1]==t&monoculture_yields[,3]==conc),3+spec_ind],2,median))
                                                                                                                         for(rep in c(1:3)){ RY               <- as.vector(as.matrix(temp3[which(temp3$concentration==conc&temp3$replicate==rep),2+spec_ind]/M))
                                                                                                                                             RYT              <- sum(RY)
                                                                                                                                             dev_exp_freq     <- RY/RYT-(1/i)
                                                                                                                                             dev_zerosum      <- RY - RY/RYT
                                                                                                                                             dY               <- sum(apply(temp3[which(temp3$concentration==conc&temp3$replicate==rep),2+spec_ind],2,mean))-sum((1/i)*M)
                                                                                                                                             dom              <- cov(dev_exp_freq,M)*(i-1)
                                                                                                                                             tdc              <- cov(dev_zerosum,M)*(i-1)
                                                                                                                                             tic              <- mean(RY-(1/i))*mean(M)*i
                                                                                                                                             total            <- temp3$total[which(temp3$concentration==conc&temp3$replicate==rep)]
                                                                                                                                             info             <- c(toString(j),i,t,rep,conc,cov(dev_exp_freq,M),cor(dev_exp_freq,M),sd(dev_exp_freq),cov(dev_zerosum,M),cor(dev_zerosum,M),sd(dev_zerosum),sd(M),mean(M),mean(RY-(1/i)),dY,dom,tdc,tic,total)
                                                                                                                                             effects_analysis <- rbind(effects_analysis,as.vector(info))}}}}}
        # save data frame          
        effects_analysis        <- data.frame(effects_analysis)
        names(effects_analysis) <- c("composition","diversity","day","replicate","concentration","cov(RY/RYT-RYe,M)","cor(RY/RYT-RYe,M)","sd(RY/RYT-RYe)","cov(RY-RY/RYT,M)","cor(RY-RY/RYT,M)","sd(RY-RY/RYT)","sd(M)","mean(M)","mean(deltaRY)","dY","dom","tdc","tic","total")
        for(i in c(2:dim(effects_analysis)[2])){effects_analysis[,i] <- as.numeric.factor(effects_analysis[,i])}
        save(effects_analysis,file="effects_analysis.RData")
                  
# 3. Calculate change in biodiversity effects

        # calculate changes in biodiveristy effects
        change_effects <- NULL
        for(i in c(2,4,6,8)){ temp         <- effects_analysis[which(effects_analysis$diversity==i),]
                              compositions <- unique(temp$composition)
                              for(j in compositions){ temp2      <- temp[which(temp$composition==j),]
                                                      for(t in c(7,14,21,28)){ temp3    <- temp2[which(temp2$day==t),]
                                                                                  reg_dom  <- lm(dom/10^9~concentration,data=temp3)
                                                                                  reg_tdc  <- lm(tdc/10^9~concentration,data=temp3)
                                                                                  reg_tic  <- lm(tic/10^9~concentration,data=temp3)
                                                                                  info     <- c(j,i,t,summary(reg_dom)$coefficients[2,1],summary(reg_dom)$coefficients[2,2],summary(reg_dom)$coefficients[2,4],summary(reg_tdc)$coefficients[2,1],summary(reg_tdc)$coefficients[2,2],summary(reg_tdc)$coefficients[2,4],summary(reg_tic)$coefficients[2,1],summary(reg_tic)$coefficients[2,2],summary(reg_tic)$coefficients[2,4])
                                                                                  change_effects <- rbind(change_effects,info)}}}


        # save data frame
        change_effects        <- data.frame(change_effects)
        names(change_effects) <- c("compositions","diversity","day","betadom","sddom","pdom","betatdc","sdtdc","ptdc","betatic","sdtic","ptic")
        for(i in c(2:12)){change_effects[,i] <- as.numeric.factor(change_effects[,i])}
        save(change_effects,file="change_effects.RData")


