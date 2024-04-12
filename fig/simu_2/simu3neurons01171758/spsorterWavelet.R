library(waveslim)
library(scatterplot3d)
source("read_simudatas.R")

spsorterWavelet <- function(filename,res.init)
  {
    
    num.data.points=40
    num.interp.data.points=64
    num.spikes = sum(res.init$num.spikes)
    num.neurons = res.init$num.temps
    
    ##spikes is a matrix
    res.read.simudatas <- read.simudatas(filename)
    spikes <- res.read.simudatas$spwave
    answer <- res.read.simudatas$answer#true spike label
    
    
    ##to interporate spikes, we need 64 datapoints to compute DWT
    interp.spikes <- array(dim=c(num.spikes,num.interp.data.points))
    dwt.spikes <- array(dim=c(num.spikes,num.interp.data.points))


    
    for(i in 1:num.spikes){
      interp.spikes[i,]<-spline(spikes[i,],n=num.interp.data.points)$y  
      dwt<- dwt(interp.spikes[i,],wf="la8",n.levels=6)
      dwt.spikes[i,] <- c(dwt$d1,dwt$d2,dwt$d3,dwt$d4,dwt$d5,dwt$d6,dwt$s6)
    }

    ##record each axes(coefficients) score, I use Standard Deviation of eace coefficients. in this version 
    axes.score <- c(1:64)

    for(i in 1:num.interp.data.points){
      ## All Coefs set like a pdf to compare each SD in same condittion
      min.val = min(dwt.spikes[,i])
      coefs <- dwt.spikes[,i] - min.val
      mu = mean(coefs)
      coefs <- coefs/mu

      axes.score[i] = sd(coefs)
      ##  Next Version I will use ks-test p-value for axes.score (How far from Gaussian Distribution)
      ##  y<-amplitude*rnorm(numInterpDataPoints,mean=mu,sd=sigma)
      ## axes.score[i] <- ks.test(ks.spike, y )$p.value
      
    }

    most.sig.axes <- c(1:64)
    axes.order <- order(axes.score)

    
    jpeg(sprintf("%s%s",filename,"Axes_Scores.jpeg"))
    plot(1:64,axes.score,type="h",main="Each Coefficients SD (sorted)",xlab="Coefs",ylab="SD")
    dev.off()
   
    jpeg(sprintf("%s%s",filename,"Axes_Score_sorted.jpeg"))
    plot(rev(sort(axes.score)),type="h",main="Each Coefficients SD (sorted)",xlab="Coefs",ylab="SD")
    dev.off()

    pdf(sprintf("%s%s",filename,"Axes_Scores.pdf"))
    plot(1:64,axes.score,type="h",main="Each Coefficients SD ",xlab="Coefs",ylab="SD")
    dev.off()
   
    pdf(sprintf("%s%s",filename,"Axes_Score_sorted.pdf"))
    plot(rev(sort(axes.score)),type="h",main="Each Coefficients SD ",xlab="Coefs",ylab="SD")
    dev.off()
    
    

    ##significant axes are orderd with SD 
    for(i in 1:64){
      most.sig.axes[i] = which(order(axes.score) == i)
    }
  
    sig.wavelet.space <- cbind(dwt.spikes[, most.sig.axes[1]],
                               dwt.spikes[ ,most.sig.axes[2]],
                               dwt.spikes[ ,most.sig.axes[3]],
                               dwt.spikes[ ,most.sig.axes[4]],
                               dwt.spikes[ ,most.sig.axes[5]],
                               dwt.spikes[ ,most.sig.axes[6]],
                               dwt.spikes[ ,most.sig.axes[7]],
                               dwt.spikes[ ,most.sig.axes[8]],
                               dwt.spikes[ ,most.sig.axes[9]],
                               dwt.spikes[ ,most.sig.axes[10]]
                               )
                          
    w.cl<-kmeans(sig.wavelet.space,num.neurons)
    write(w.cl$cluster,file=paste(filename,"10coefsKmean.dat"))

    Angle=200
    cex.val=1.1
    ##plot in wavelet.space
    jpeg(sprintf("%s%s",filename,"wavelet_trueans.jpeg"))
    par(pch=20)
    par(cex=cex.val)
    p3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       color=answer,
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="plot spikes in  3D space (colored by answer)")

    dev.off()

    jpeg(sprintf("%s%s",filename,"wavelet_kmeans10.jpeg"))
    par(pch=20)
    w3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       color=w.cl$cluster,
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="wavelet Tr & colored by kmeans in 10D space")
    dev.off()
    ##no colored,  display for master thesis 
    jpeg(sprintf("%s%s",filename,"no_colored_wavelet_kmeans10.jpeg"))
    par(pch=20)
    par(cex=cex.val)
    w3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="plot spikes in 3D space")
    dev.off()

##no colored,  display for master thesis 
    pdf(sprintf("%s%s",filename,"no_colored_wavelet_kmeans10.pdf"))
    par(pch=20)
    par(cex=cex.val)
    w3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="plot spikes in 3D space")
    dev.off()
    
    pdf(sprintf("%s%s",filename,"wavelet_trueans.pdf"))
    par(pch=20)
    par(cex=cex.val)
    p3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       color=answer,
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="plot spikes in  3D space (colored by answer)")

    dev.off()
    pdf(sprintf("%s%s",filename,"wavelet_kmeans10.pdf"))
    par(pch=20)
    par(cex=cex.val)
    w3d<-scatterplot3d(sig.wavelet.space[1:num.spikes,1],
                       sig.wavelet.space[1:num.spikes,2],
                       sig.wavelet.space[1:num.spikes,3],
                       color=w.cl$cluster,
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="plot spikes in  3D space (colored by k-means)")
    dev.off()

    ##check the wavelet answers
    w.sp1.1=length(which(w.cl$cluster[1:num.spikes] == 1))
    w.sp1.2=length(which(w.cl$cluster[1:num.spikes] == 2))
    w.sp1.3=length(which(w.cl$cluster[1:num.spikes] == 3))
    ##in spite of reach the answer,w.cl$cluster may be assigned
    ##as a number which are different with answer number.
    ##num.read = 100
    ##num.spikes.ans <- c(0,0,0)
    ##num.match <- matrix(rep(0,3*3),nrow=3,ncol=3)


    ##for(i in 1:num.read)
    ##  num.match[answer[i],w.cl$cluster[i]] = num.match[answer[i],w.cl$cluster[i]] + 1
    



    ##I compute kmeans with Top 10 Axes(coeffs),but result is bad. 
    ##So, I compute kmeans-clustering using by only Top 3 axes.
    top3.wavelet.space <- cbind(dwt.spikes[,most.sig.axes[1]],
                                dwt.spikes[,most.sig.axes[2]],
                                dwt.spikes[,most.sig.axes[3]]
                                )

    top3.cl<-kmeans(top3.wavelet.space,num.neurons)


    p3d<-scatterplot3d(dwt.spikes[1:100,most.sig.axes[1]],
                       dwt.spikes[1:100,most.sig.axes[2]],
                       dwt.spikes[1:100,most.sig.axes[3]],
                       color=1,
                       angle=Angle,
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3",
                       main="wavelet Tr")


                                        #png("wavelet_kmeans.png")
    jpeg(sprintf("%s%s",file=filename,"wavelet_kmeans.jpeg"))
    w3d<-scatterplot3d(dwt.spikes[1:100,most.sig.axes[1]],
                       dwt.spikes[1:100,most.sig.axes[2]],
                       dwt.spikes[1:100,most.sig.axes[3]],
                       color=top3.cl$cluster,
                       angle=Angle,
                       main="wavelet Tr & kmeans in top3 wavelet space",
                       xlab="Coef 1",ylab="Coef 2",zlab="Coef 3"
                       )
    dev.off()

    pdf(sprintf("%s%s",file=filename,"wavelet_kmeans.pdf"))
    w3d<-scatterplot3d(dwt.spikes[1:100,most.sig.axes[1]],
                       dwt.spikes[1:100,most.sig.axes[2]],
                       dwt.spikes[1:100,most.sig.axes[3]],
                       color=top3.cl$cluster,
                       angle=Angle,
                       main="wavelet Tr & kmeans in top3 wavelet space",
                        xlab="Coef 1",ylab="Coef 2",zlab="Coef 3"
                       )
    dev.off()
    
    top3.sp1.1=length(which(top3.cl$cluster[1:100] == 1))
    top3.sp1.2=length(which(top3.cl$cluster[1:100] == 2))
    top3.sp1.3=length(which(top3.cl$cluster[1:100] == 3))

    
    write(answer,file=sprintf("%s%s",filename,"true_answer.dat"))
    write(w.cl$cluster,file=sprintf("%s%s",filename,"top10Dkmeans_spikelabel.dat"))
    write(top3.cl$cluster,file=sprintf("%s%s",filename,"top3Dkmeans_spikelabel.dat"))

    
    check.answer <- function(answer.vec,cluster.vec,num.vec,num.label=2)
      {
        ##count increases count[j,k] +1 when answer.vec[i] == j and cluster.vec[i] == k,
        ##which means i-th spike is labeled k.
        count <- matrix(rep(0,num.label*num.label),ncol=num.label,nrow=num.label)
    
        for(i in 1:num.vec){

          for(j in 1:num.label){
            if(j == answer.vec[i]){
          
              for(k in 1:num.label)
                if(k == cluster.vec[i]){               
                  count[j,k] = count[j,k] + 1
                }
            }
          }
        }
              
        return(list("count"=count))
      }
 
    
   

    
    
    return(list("count"=check.answer(answer,w.cl$cluster,100,3)$count,
                "interp.spikes"=interp.spikes,
                "dwt.spikes"=dwt.spikes,
                "answer"=answer))
    
  }
