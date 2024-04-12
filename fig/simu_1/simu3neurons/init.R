init <- function(init.filename,init.SNR)
  {
    res.temp <-  spwave.template()
    templates <-  rbind(res.temp$sp1,res.temp$sp2,res.temp$sp3)
    
    options(digits=20)

    SAMPLE.PER.REC = 25000.0
    
    ##one spike waveform has 40 data points
    DATAPOINTS = 40
    
    ##number of spikes
    num.temps = 3#number of neurons
    num.spikes <- c(35,35,30)#number of spikes in each neurons

    
    ##firing pattern
    ##spikes are generated from poisson distributions
    ##lambda is a mean of a poisson distribution in each neurons
    
    lambda <- c(1500,1700,2000)
    ##rpois() function return integer values(too large for isi), so that I convert from integer
    ##to decimal value by using isi.time.part.  
    isi.time.part = 0.0001
    ## neuron-[i]'s mean of isi(firing.rate) is isi.time.part*lambda[i]s (s^-1)

    ##Constant value for spike waveform Amplitude
    Amp.const = 1

    ##Noise, I add the noise to each spike waveform. 
    SNR = init.SNR #SN ratio
     
    ## filenme for save datas
    filename =init.filename
    
    ##attenuation
    ##I suppose neurons's spike amplitude  are attenuated with
    ##this expression '1-\delta * (exp( - \lambda * isi)' .
    atten.delta <- c(0.8,0.8,0.8)
    #atten.lambda <- c(0.00002,0.000015)#referenced to Yoshikawa's master thesis
    atten.lambda <- c(100,100,150)#referenced to Yoshikawa's master thesis
    ## SNR  noise  
    atten.SNR = 100
    ## observed.time <-  seq(0,50,by = 1/SAMPLE.PER.REC)
    count.overlaped = 0
    
    return(list("filename"=filename,
                "SAMPLE.PER.REC"=SAMPLE.PER.REC,
                "DATAPOINTS"=DATAPOINTS,
                "templates"=templates,
                "num.temps"=num.temps,
                "num.spikes"=num.spikes,
                "pois.lambda.neurons"=lambda,
                "convert.time"=isi.time.part,
                "Amp.const"=Amp.const,
                "SNR"=SNR,
                "atten.delta"=atten.delta,
                "atten.lambda"=atten.lambda,
                "atten.SNR"=atten.SNR
                ))

  }
