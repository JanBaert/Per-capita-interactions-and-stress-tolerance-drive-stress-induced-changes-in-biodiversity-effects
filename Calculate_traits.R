# load average parameters estimates for scenarios 4 and 5
load("estimates_4.RData")
load("estimates_5.RData")

# load data biovolumes and calculated changes biodiveristy effects
load("data_biov.RData")
load("change_effects.RData")

# load vector with the corresponding compositions
load("compositions_list.RData")

# load average monoculture yields at day 7, 14, 21 and 28 (i.e. colums) at 0 ppb and 250 ppb
species     <- c("8","10","17","18","E18","23","31","32")
load("M_yield.RData")

# make trait-based dataset
  data        <- change_effects
  data_traits <- NULL
  times       <- c(7,14,21,28)
  for(k in 1:124){  
                    # estimate weighted sensitivity and per-capita interaciton strength from scenario 4
                    temp   <- data[k,]
                    spec   <- strsplit(toString(unlist(temp[1])),"/",fixed=TRUE)[[1]]
                    n_spec <- length(spec)
                    est    <- estimates_4[which(compositions==toString(unlist(temp[1]))),1:(n_spec*(n_spec+3))]   
                    x      <- est[1:(n_spec*(n_spec-1))]
                    K      <- est[((n_spec*(n_spec+1))+1):(n_spec*(n_spec+2))]
                    int    <- mean(x)
                    biov   <- data_biov[which(data_biov$composition==toString(unlist(temp[1]))&data_biov$day==unlist(temp[3])&data_biov$concentration==0),13:20]     
                    if(n_spec!=8){biov <- biov[,which(biov[1,]!=0)]}
                    biov   <- sapply(1:nrow(biov),function(x) biov[x,]/sum(biov[x,]))
                    RY     <- sapply(1:nrow(biov),function(x) mean(as.numeric(biov[x,])))
                    interact  <- matrix(1,ncol=n_spec,nrow=n_spec)
                    count  <- 1
                    for(i in c(1:(n_spec-1)))
                      { for(j in c((i+1):n_spec))
                          { interact[i,j] <- x[count]
                            count         <- count+1
                            interact[j,i] <- x[count]
                            count         <- count+1 }}
                    for(i in 1:n_spec){interact[i,-i] <- interact[i,-i]/K[i]} 
                    interact <- sapply(1:n_spec,function(x) mean(interact[x,-x]))
                    int_w_4  <- sum(RY*interact)
                    ind      <- which(species%in%spec)
                    M_0      <- M_yield[ind,which(times==28),1]
                    M_250    <- M_yield[ind,which(times==28),2]
                    sens     <- M_250/M_0
                    sens_w   <- sum(RY*sens)
                    sens     <- mean(sens)
                    
                    # estimate changes in per-capita interaction strenght from scenario 5
                    est    <- estimates_5[which(compositions==toString(unlist(temp[1]))),]   
                    x      <- est[1:(n_spec*(n_spec-1))]
                    x_250  <- est[(n_spec*(n_spec-1))+1:(n_spec*(n_spec-1))]
                    K      <- est[(2*(n_spec*(n_spec-1))+2*n_spec)+1:n_spec]
                    K_250  <- est[(2*(n_spec*(n_spec-1))+3*n_spec)+1:n_spec]
                    int    <- mean(x)
                    interact  <- matrix(1,ncol=n_spec,nrow=n_spec)
                    count  <- 1
                    for(i in c(1:(n_spec-1)))
                      { for(j in c((i+1):n_spec))
                          { interact[i,j] <- x[count]
                            count         <- count+1
                            interact[j,i] <- x[count]
                            count         <- count+1 }}
                    for(i in 1:n_spec)
                    {interact[i,-i] <- interact[i,-i]/K[i]} 
                    interact <- sapply(1:n_spec,function(x) mean(interact[x,-x]))
                    int_w    <- sum(RY*interact)
                    interact  <- matrix(1,ncol=n_spec,nrow=n_spec)
                    count  <- 1
                    for(i in c(1:(n_spec-1)))
                      { for(j in c((i+1):n_spec))
                        { interact[i,j] <- x_250[count]
                          count         <- count+1
                          interact[j,i] <- x_250[count]
                          count         <- count+1 }}
                    for(i in 1:n_spec)
                    {interact[i,-i] <- interact[i,-i]/K_250[i]} 
                    interact    <- sapply(1:n_spec,function(x) mean(interact[x,-x]))
                    int_w_250   <- sum(RY*interact)
                    d_int_w     <- int_w_250-int_w
                    d_int_w     <- sign(d_int_w)*log10(abs(d_int_w))
                    data_traits <- rbind(data_traits,c(toString(unlist(temp[1])),unlist(temp[c(2,3,4,7,10)]),log10(int_w_4),sens_w,d_int_w))}
                  
# save traits data frame
data_traits <- data.frame(comp=data_traits[,1],diversity=as.numeric(data_traits[,2]),day=as.numeric(data_traits[,3]),betadom=as.numeric(data_traits[,4]),betatdc=as.numeric(data_traits[,5]),betatic=as.numeric(data_traits[,6]),int_w=as.numeric(data_traits[,7]),sens_w=as.numeric(data_traits[,8]),d_int_w=as.numeric(data_traits[,9]))
data_traits <- data_traits[-which(is.na(data_traits$int)),]
save(data_traits,file="data_traits.RData")


