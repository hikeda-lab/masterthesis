read.simudatas <- function(filename)
  {
    ##Read 3 neurons spike waveforms datas  as 3 vector, 
    ##and rewrite it as an array.
    ##
    SAMPLES = 40
    num.read = 100 #I read 100 spikes
    sp1 <- matrix(scan(file=sprintf("%s%s",filename,"spike_waveform_1.dat")),ncol=SAMPLES)
    sp2 <- matrix(scan(file=sprintf("%s%s",filename,"spike_waveform_2.dat")),ncol=SAMPLES)
    sp3 <- matrix(scan(file=sprintf("%s%s",filename,"spike_waveform_3.dat")),ncol=SAMPLES)
    ##spike labeles , it is arraied with time order
    answer <- scan(file=sprintf("%s%s",filename,"answer.dat"))

    ##make the list
    spwaveforms <- list(sp1,sp2,sp3)

    ##rearranged each spike's spike waveforms with a time order
    total = length(sp1)+length(sp2)+length(sp3)# total number of  spikes
    spwave.time.order <- matrix(nrow=num.read,ncol=SAMPLES)

    
    count.sp <- c(1,1,1)
    
      for(i in 1:num.read)
        {
          spwave.time.order[i,] = spwaveforms[[ answer[i] ]][ count.sp[answer[i]], ]
          count.sp[answer[i]] = count.sp[answer[i]] + 1
        }

    
    #browser()
    ##return:
    ##spwave.time.order: matrix
    return(list("spwave"=spwave.time.order,
                "answer"=answer[1:num.read]))
  }
