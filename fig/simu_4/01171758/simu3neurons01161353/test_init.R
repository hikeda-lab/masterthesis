init <- function()
  {
    res.temp <-  template()
    templates <-  rbind(res.temp$A,res.temp$B)
    
    options(digits=20)

    SAMPLE.PER.REC = 25000.0
    
    ##one spike waveform has 40 data points
    DATAPOINTS = 40
    
    ##number of spikes
    num.temps = 2#number of neurons
    num.spikes <- c(50,50)#number of spikes in each neurons

    
    ##firing pattern
    ##spikes are generated from poisson distributions
    ##lambda is a mean of a poisson distribution in each neurons
    lambda <- c(20000,15000)
    ##rpois() function return integer values(too large for isi), so that I convert from integer
    ##to decimal value by using isi.time.part.  
    isi.time.part = 0.0001


    ##Constant value for spike waveform Amplitude
    Amp.const = 50

    ##Noise, I add the noise to each spike waveform. 
    SNR = 0.1 #SN ratio
     
    ## filenme for save datas
    filename ="test.dat"
    
    ##attenuation
    ##I suppose neurons's spike amplitude  are attenuated with
    ##this expression '1-\delta * (exp( - \lambda * isi)' .
    atten.delta <- c(1.0,0.5)
    atten.lambda <- c(100,200)

    
    ## observed.time <-  seq(0,50,by = 1/SAMPLE.PER.REC)
    count.overlaped = 0
    
    return(list("SAMPLE.PER.REC"=SAMPLE.PER.REC,
                "DATAPOINTS"=DATAPOINTS,
                "templates"=templates,
                "num.temps"=num.temps,
                "num.spikes"=num.spikes,
                "pois.lambda.neurons"=lambda,
                "convert.time"=isi.time.part,
                "Amp.const"=Amp.const,
                "SNR"=SNR,
                "atten.delta"=atten.delta,
                "atten.lambda"=atten.lambda
                ))

  }
