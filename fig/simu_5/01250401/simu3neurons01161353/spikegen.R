

                                  
# sample per sec = 1/25000(s) 1poits=1/25000(s)
# template = 40 points
#

spikegen <- function(res.init)
  {
    #res.init <- init()

    filename <- res.init$filename
    templates <-  res.init$templates
    
    options(digits=20)
    set.seed(101)
    
    SAMPLE.PER.REC = res.init$SAMPLE.PER.REC
    
    ##one spike waveform has 40 data points
    DATAPOINTS = res.init$DATAPOINTS
    
    ##number of spikes
    num.temps = res.init$num.temps#number of neurons
    num.spikes <- res.init$num.spikes#number of spikes in each neurons

    
    ##firing pattern
    ##spikes are generated from poisson distributions
    ##lambda is a mean of a poisson distribution in each neurons
    lambda <- res.init$pois.lambda.neurons
    ##rpois() function return integer values(too large for isi), so that I convert from integer
    ##to decimal value by using isi.time.part.  
    isi.time.part = res.init$convert.time


    ##Constant value for spike waveform Amplitude
    A.const = res.init$Amp.const

    ##Noise, I add the noise to each spike waveform. 
    SNR = res.init$SNR #SN ratio
     
    ## filenme for save datas
    
    ##attenuation
    ##I suppose neurons's spike amplitude  are attenuated with this expression '1-\delta * (exp( - \lambda * isi)' .
    atten.delta <- res.init$atten.delta
    atten.lambda <- res.init$atten.lambda
    atten.SNR <- res.init$atten.SNR
    
    attenuation <- function(lapsed.time,
                            atten.delta,atten.lambda,atten.SNR,
                            neuron.number,datapoints)
      {
        
        delta = atten.delta[neuron.number]
        lambda = atten.lambda[neuron.number] 
        val = 1 - delta*exp(-(lambda * lapsed.time))
        noise.sd = 1/atten.SNR
        
        return( rep( val +  rnorm(1,mean=0,sd=noise.sd*val) ,datapoints))
        #return(rep(1,datapoints))
      }
    
    ## observed.time <-  seq(0,50,by = 1/SAMPLE.PER.REC)
    count.overlapped = 0
    
    ##generate spike's firing times from Poisson distributions 
    #real time
    
    isi <- array(dim=c(num.temps,max(num.spikes)))
    for(i in 1:num.temps)
      for(j in 1:num.spikes[i]){
      isi[i,j] <-  isi.time.part * rpois(1,lambda[i]) 
     }
  
    ##real time on data point
    isi.point <- array(dim=c(num.temps,max(num.spikes))) 
    for(i in 1:num.temps){
      v <-  round(isi[i,]*SAMPLE.PER.REC,digits=100000)
      isi.point[i,] <- v
    }
   
    ##firing time
    firing.time <-  array(dim=c(num.temps,max(num.spikes)))
    ##convert time data to data point data
    firing.time.point <- array(dim=c(num.temps,max(num.spikes)))
      
          
    for(i in 1:num.temps){
      firing.time[i,1] = isi[i,1]
      firing.time.point[i,1] = isi.point[i,1]
  
      for(j in 2:num.spikes[i]){
        firing.time[i,j] = firing.time[i,j-1] + isi[i,j]
        firing.time.point[i,j] = firing.time.point[i,j-1] + isi.point[i,j]
      }  
    }
          
  
    
    
          
    ##In  experiments , we observe spikes in order of their firing-time
    all.firing.time <- firing.time[1,1:num.spikes[1]]
    all.firing.time.point <- firing.time.point[1,1:num.spikes[1]]
    if(num.temps > 1){
      for(i in 2:num.temps){
        x <- firing.time[i,1:num.spikes[i]]
        all.firing.time <- c(all.firing.time,x)
        y <- firing.time.point[i,1:num.spikes[i]]
        all.firing.time.point <- c(all.firing.time.point,y)
      }
    }
    all.firing.time.order <- order(all.firing.time)
    all.firing.time <- sort(all.firing.time)
    all.firing.time.point <- sort(all.firing.time.point)
        
    ## Generate spike waveforms from templates    
    templates <- A.const*templates
    
    x <- mean(templates[1,])
    mean.amp <- x
    y <- sd(templates[1,])
    sd.amp <- y
    if(num.temps > 1){
      for(i in 2:num.temps){
        x <- mean(templates[i,])
        mean.amp <- c(mean.amp,x)
        y <- sd(templates[i,])
        sd.amp <- c(sd.amp,y)
      }
    }

    
    ## spike waveform
    spike.waveform <- array(dim=c(num.temps,max(num.spikes),DATAPOINTS))
    spike.waveform.time.order <- array(dim=c(sum(num.spikes),DATAPOINTS))
    spike.label <-  array(dim=c(num.temps,max(num.spikes)))
    for(i in 1:num.temps)
      spike.label[i,] <- rep(i,max(num.spikes))
    
    ##Plot spike waveforms
    ## waveforms are atenuated with elapsed time and  added with noise 
    for(i in 1:num.temps){
      spike.waveform[i,1,] <- templates[i,]
      for(j in 2:num.spikes[i]){
        #Attenuation <-  attenuation(isi[[i]][j],isi.time.part*lambda[i],DATAPOINTS)
        Attenuation <-  attenuation(isi[i,j],atten.delta,atten.lambda,atten.SNR,i,DATAPOINTS)
        v <- rep(0,DATAPOINTS)
        if(isi[i,j] > 1.0){ v <- templates[i,] + rnorm(DATAPOINTS,0,sd.amp[i]/SNR)}
        else{ v<- spike.waveform[i,j-1,]}
        ## Spike Amp = Template * Atten  +  Noise
        spike.waveform[i,j,] <- Attenuation * v + rnorm(DATAPOINTS,0,sd.amp[i]/SNR)
        #spike.waveform[i,j,] <- Attenuation * v 
       
     }
   }


    buf.spike.waveform <-  spike.waveform

  
    
   #A spike is overlapped with other spike ,
   #if other spike fires within  a spike waveform's data points(40 points) from the time a spike fires.

    #a spike fires: firing.time.point[[i]][k]
    #other spike fires: firing.time.point[[j]][l]
    #
    #I wonder you image the situation
    #with makeing a spike firing time diagram.
    #
   
   for(i in 1:num.temps){
     for(j in 1:num.temps){
         for(k in 1:num.spikes[i]){
           for(l in 1:num.spikes[j]){            
              ##a spike fires after the other spike fires, as a result spikes are overlaped. 
             if((firing.time.point[i,k] > firing.time.point[j,l]) &&
                (firing.time.point[i,k] < firing.time.point[j,l] + DATAPOINTS))
               {
                 count.overlapped = count.overlapped + 1
                 ##if overlped,spike.label tells you other spike's label 
                 #spike.label[i,k] = j
                 ##overlapped                                                                       
                 for(x in 1:(firing.time.point[j,l] + DATAPOINTS - firing.time.point[i,k])){
                   end=firing.time.point[j,l] + DATAPOINTS - firing.time.point[i,k]
                   y=firing.time.point[i,k] - firing.time.point[j,l] + x
                   ##draw a spike overlapped with the othe spike
                   spike.waveform[i,k,x] = max( buf.spike.waveform[i,k,x],buf.spike.waveform[j,l,y])
                   ##draw the othe spike overlapped with a spike
                   spike.waveform[j,l,y] = max( buf.spike.waveform[i,k,x],buf.spike.waveform[j,l,y])
                 }
                     
                    
                 for(x in (end+1):(firing.time.point[j,l] + DATAPOINTS - firing.time.point[i,k]+1):DATAPOINTS){
                   ##draw a spike not overlapped
                   spike.waveform[i,k,x] = buf.spike.waveform[i,k,x]
                   ##draw the other spike not overlapped
                  # y=firing.time.point[i,k] - firing.time.point[j,l] + x
                  # spike.waveform[j,l,y] = buf.spike.waveform[j,l,y]
                 }
               }


           }
         }
       }
   }
         
                  
   
  
 
  
    ##Max amplitude
    max.amp <- array(dim=c(num.temps,max(num.spikes)))
    
    
    for(i in 1:num.temps)
      for(j in 1:num.spikes[i])
        max.amp[i,j] <- max(spike.waveform[i,j,])

    max.amp.neuron <- max.amp#max amplitude in one spike waveform for each neuron
    max.amp.vec <- c(max.amp[1,1:num.spikes[1]])
    if(num.temps>1)
      for(i in 2:num.temps)
        max.amp.vec <- c(max.amp.vec,max.amp[i,1:num.spikes[i]])
    
    ##ordered max.amp with all.firing.time.order
    max.amp <- max.amp.vec[all.firing.time.order]
    

    ##Spike label(this is the answer for a spike sorter
    answer <- array(dim=c(num.temps,max(num.spikes)))
    
    for(i in 1:num.temps)
      for(j in 1:num.spikes[i])
        answer[i,j] <- spike.label[i,j]

    answer.vec <- c(answer[1,1:num.spikes[1]])
    if(num.temps>1)
      for(i in 2:num.temps)
        answer.vec <- c(answer.vec,answer[i,1:num.spikes[i]])
      
    
    answer <- answer.vec
    ##ordered max.amp with all.firing.time.order
    answer <- answer[all.firing.time.order]

    amp.ratio <- array(dim=c(num.temps,max(num.spikes)))

    for(i in 1:num.temps){
      amp.ratio[i,1] = 1
      if(num.spikes[i]>1)
        for(j in 2:num.spikes)
          amp.ratio[i,j] = max.amp.neuron[i,j-1]/max.amp.neuron[i,j]
    }
    
    ##save datas
    save.datas(filename,
               firing.time.point,
               firing.time,
               all.firing.time,
               max.amp,
               max.amp.neuron,
               isi.point,
               isi,
               spike.waveform,
               spike.label,
               count.overlaped,
               answer
               )

    
    
    return(list("time.pts"=firing.time.point,
                "time"=firing.time,
                "all.time.pts"=all.firing.time.point,
                "all.time"=all.firing.time,
                "max.amp"=max.amp,
                "max.amp.neuron"=max.amp.neuron,
                "isi.pts"=isi.point,
                "isi"=isi,
               "wf"=spike.waveform,
                "label"=spike.label,
                "answer"=answer,
                "count.overlapped" = count.overlapped
                ))

  }

save.datas <- function(filename,
                       firing.time.point,#the vector
                       firing.time,#the vectort
                       all.firing.time,#the vector
                       max.amp,# the vector
                       max.amp.chan,#the list
                       isi.point,
                       isi,
                       spike.waveform,
                       spike.label,
                       count.overlapped,
                       answer
                       )
{
  
  
  write(all.firing.time,file=sprintf("%s%s",filename,"all_firing_time.dat"))
  write(max.amp,file=sprintf("%s%s",filename,"max_amp.dat"))

  
  write(firing.time.point,file=sprintf("%s%s",filename,"firing_time_point_1.dat"))
 

  write(firing.time,file=sprintf("%s%s",filename,"firing_time_1.dat"))
 

  write(isi.point,file=sprintf("%s%s",filename,"isi_point_1.dat"))
 
 
  write(isi,file=sprintf("%s%s",filename,"isi_1.dat"))
 
  spike.waveform.1 <- matrix(spike.waveform[1,,],ncol=40)
  spike.waveform.2 <- matrix(spike.waveform[2,,],ncol=40)
  spike.waveform.3 <- matrix(spike.waveform[3,,],ncol=40)
  
  
  write(spike.waveform.1,file=sprintf("%s%s",filename,"spike_waveform_1.dat"))
  write(spike.waveform.2,file=sprintf("%s%s",filename,"spike_waveform_2.dat"))
  write(spike.waveform.3,file=sprintf("%s%s",filename,"spike_waveform_3.dat"))

  write(max.amp.chan,file=sprintf("%s%s",filename,"max_amp_1.dat"))
 

  write(answer,file=sprintf("%s%s",filename,"answer.dat"))

  
}
