
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

##for(i in 1:10)
##  {
##    res.init <- init(filenames[i],SNR.vec[i])
##    res.spikegen <- spikegen(res.init)

##    res <- spsorterWavelet(filenames[i],res.init)

##    res.spikegen.count.overlapped[i] <- res.spikegen$count.overlapped
##    res.count[i,,] <- res$count
    
    ##browser()
##  }

i = 1
  {
    res.init <- init(filenames[i],SNR.vec[i])
    res.spikegen <- spikegen(res.init)

    times <- system.time(res <- spsorterWavelet(filenames[i],res.init))

    res.spikegen.count.overlapped[i] <- res.spikegen$count.overlapped
    res.count[i,,] <- res$count
    
    ##browser()

    data.samples <- 1:64
    cex.val = 1.5
    
    pdf("spikes1.pdf")
    par(cex=cex.val)
    plot(data.samples,res$interp.spikes[1,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()
    
    jpeg("spikes1.jpg")
    par(cex=cex.val)
     par(lwd=5)
    plot(data.samples,res$interp.spikes[1,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()
    
    pdf("dwt_spikes1.pdf")
    par(cex=cex.val)
    par(lwd=5)
    plot(data.samples,res$dwt.spikes[1,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coeff Number")
    dev.off()
    
    jpeg("dwt_spikes1.jpg")
    par(cex=cex.val)
    plot(data.samples,res$dwt.spikes[1,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coef Number")
    dev.off()

    ##spike 2
    pdf("spikes2.pdf")
    par(cex=cex.val)
    plot(data.samples,res$interp.spikes[2,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()
    
    jpeg("spikes2.jpg")
    par(cex=cex.val)
    plot(data.samples,res$interp.spikes[2,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()

    pdf("dwt_spikes2.pdf")
    par(cex=cex.val)
     par(lwd=5)
    plot(data.samples,res$dwt.spikes[2,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coeff Number")
    dev.off()
    
    jpeg("dwt_spikes2.jpg")
    par(cex=cex.val)
    par(lwd=5)
    plot(data.samples,res$dwt.spikes[2,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coeff Number")
    dev.off()


    ##spike 3
    pdf("spikes3.pdf")
    par(cex=cex.val)
    par(lwd=5)
    plot(data.samples,res$interp.spikes[3,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()
    
    jpeg("spikes3.jpg")
    par(cex=cex.val)
    plot(data.samples,res$interp.spikes[3,],lwd=5,type="l",
         ylab="Amplitude",xlab="Data Samples")
    dev.off()

 
    pdf("dwt_spikes3.pdf")
    par(lwd=5)
    par(cex=cex.val)
    plot(data.samples,res$dwt.spikes[3,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coeff Number")
    dev.off()
    
    jpeg("dwt_spikes3.jpg")
    par(lwd=5)
    par(cex=cex.val)
    plot(data.samples,res$dwt.spikes[3,],lwd=5,type="h",
         ylab="Amplitude",xlab="Coeff Number")
    dev.off()


    
    
    
    pdf("all.spikes.pdf")
    
    plot(data.samples,res$interp.spikes[1,],type="l",ylim=c(-1100,800),
         ylab="Amplitude",xlab="Data Samples")

    for(j in 2:100)
      points(data.samples,res$interp.spikes[j,],col=res$answer[j],type="l")
    dev.off()
  }

