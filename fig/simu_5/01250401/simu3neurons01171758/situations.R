
source("spwave_template.R")
source("init.R")
source("spikegen.R")
source("spsorterWavelet.R")
##set random seed
set.seed(101)

##spike-sorting simulation in 3 neurons of 10 different SNRs 
filenames <- vector(length=10)
for( i in 1:10)
  filenames[i] <- sprintf("./data/SNR%s",i)

SNR.vec <- c(100,90,80,70,60,50,40,30,20,10)
eval.vec <- vector(length=10)

res.spikegen.count.overlapped <- vector(length=10) 

res.count <- array(dim=c(10,3,3))

for(i in 1:10)
  {
    res.init <- init(filenames[i],SNR.vec[i])
    res.spikegen <- spikegen(res.init)

    res <- spsorterWavelet(filenames[i],res.init)

    res.spikegen.count.overlapped[i] <- res.spikegen$count.overlapped
    res.count[i,,] <- res$count
    
    ##browser()
  }

i = 1
  {
    res.init <- init(filenames[i],SNR.vec[i])
    res.spikegen <- spikegen(res.init)

    res <- spsorterWavelet(filenames[i],res.init)

    res.spikegen.count.overlapped[i] <- res.spikegen$count.overlapped
    res.count[i,,] <- res$count
    
    ##browser()
  }

