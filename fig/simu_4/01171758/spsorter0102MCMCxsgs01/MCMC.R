library(adapt)
 

spikesort.MCMC <- function()
  {
    ##set variables
    init.MCMC <- init()
    num.neurons = init.MCMC$neurons
    num.spikes = init.MCMC$spikes
    MCMC.steps = init.MCMC$MCMC.steps

    
    ##observed spike datas
    firing.time <- init.MCMC$firing.time
    amp.time <- init.MCMC$amp.time
    spike.label <- init.MCMC$spike.label
    spike.label.history <- matrix(nrow=MCMC.steps,ncol=num.spikes)

    ##Amplitude parameters
    lambda = init.MCMC$lambda;lambda.min = init.MCMC$lambda.min;lambda.max = init.MCMC$lambda.max
    ##lambda history for each neuron
    lambda.history <- matrix(nrow=MCMC.steps*num.spikes,ncol=num.neurons)
    delta = init.MCMC$delta;delta.min = init.MCMC$delta.min;delta.max = init.MCMC$delta.max   
    ##delta.history for each neuron 
    delta.history <- matrix(nrow=MCMC.steps*num.spikes,ncol=num.neurons)

    
    
	
    ##Energy
    energy = 0
    energy.history <- matrix(nrow=MCMC.steps*num.spikes,nco=1)
    
    ##for chooose proper parameter(sd) in generate.random.walk
    ##a proper sd generates an accepted delta about 40%. 
    count.accept.delta = 0;count.reject.delta =0
    count.accept.lambda = 0;count.reject.lambda =0

    
    ##Markov Chain Monte Calro subroutine    
    for(now.step in 1:MCMC.steps)
      {
        cat("steps=",now.step,"\n")
        ##One transition ends when all spike's label and all parameters' trasnsision are comletely done 
        cat("spike.label.transition \n")
        ##all spike's label transition(one step for spike label transition)
        for(now.label.number in 1:num.spikes)
          {
             
            cat(now.label.number," ")
            before.spike.label <- spike.label
            proposed.spike.label <- spike.label
            proposed.spike.label[now.label.number] <-
              generate.random.walk.spike.label(spike.label, now.label.number , num.neurons)
            
            ##accept or reject the transision
            u = log(runif(1)) #to use log lilklihood,we need a log value by generated from unif distribution
           
            
            Accept.proposed.label = Accept.prob.spike.label(proposed.spike.label,before.spike.label,now.label.number,
                                    firing.time,amp.time,
                                    delta,lambda,
                                    delta.min,delta.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
            cat(Accept.proposed.label,"  ")
            
            if( Accept.proposed.label > u) # accept the prposed spike label
              {
                spike.label <- proposed.spike.label
                spike.label.history[now.step,] <- spike.label
                cat("spike label changed",before.spike.label[now.label.number],"->",proposed.spike.label[now.label.number],"\n")
              }
            else
              {
                spike.label.history[now.step,] <- spike.label
                cat("spike label not changed",before.spike.label[now.label.number],"->",spike.label[now.label.number],"\n")
              }
            ##When all spike's label are inquired (accept or reject),add spikelabel to history
            spike.label.history[now.step,] <- spike.label
          

            #cat("delta transition",(now.step-1)*num.spikes+now.label.number,"\n")
            ##delta for each neuron
            for(label in 1:num.neurons){
              u.delta = log(runif(1)) #use a log likelihood  
              before.delta <- delta
              before.delta[label] <- delta[label]
              proposed.delta <- delta
              proposed.delta[label] =
                generate.random.walk.delta(label,spike.label,delta,lambda,delta.min,delta.max,now.step)
              ##generate.piece.wise.delta(label,firing.time,amp.time,
              ##spike.label,delta,lambda,delta.min,delta.max,now.step)

              Accept = Accept.prob.delta(label,proposed.delta,before.delta,
                firing.time,amp.time,spike.label,
                delta,delta.min,delta.max,lambda,lambda.min,lambda.max,
                num.neurons,num.spikes,
                now.step)

          
              if( Accept > u.delta)
                {
                  delta <- proposed.delta
                  delta.history[(now.step-1)*num.spikes + now.label.number,label] <- delta[label]

                  #cat("delta changed, neuron=",label," ",before.delta[label],"->",delta[label],"\n")
                  
                  count.accept.delta = count.accept.delta + 1
                }
              else
                {
                  delta.history[(now.step-1)*num.spikes + now.label.number,label] <- delta[label]
                  #cat("delta not changed, neuron=",label," ",before.delta[label],"->",delta[label],"\n")
                  count.reject.delta = count.reject.delta + 1
                }
            }
            
            #cat("lambda transition ",(now.step-1)*num.spikes+now.label.number,"\n")
        
            ##amplitude parameter lambda for each neuron transision
            for(label in 1:num.neurons){
              ##accept or reject the proposed lambda
              u.lambda = log(runif(1))
              before.lambda <- lambda
              before.lambda[label] <- lambda[label]
              proposed.lambda <- lambda
              ##proposed.lambda =generate.piece.wise.lambda(label,firing.time,amp.time,spike.label,delta,lambda,lambda.min,lambda.max,now.step)
              proposed.lambda[label] =
                generate.random.walk.lambda(label,spike.label,delta,lambda,lambda.min,lambda.max,now.step)
          
        
              Accept = Accept.prob.lambda(label,proposed.lambda,before.lambda,
                firing.time,amp.time,spike.label,
                delta,delta.min,delta.max,
                lambda,lambda.min,lambda.max,
                num.neurons,num.spikes,
                now.step)
          
          
              if( Accept > u.lambda)
                {
                  lambda = proposed.lambda
                  lambda.history[(now.step-1)*num.spikes+now.label.number,label] <- lambda[label]
                  #cat("lambda changed, neuron=",label," ",before.lambda[label],"->",lambda[label],"\n")
                  count.accept.lambda = count.accept.lambda + 1
                }
              else
                {
                  lambda.history[(now.step-1)*num.spikes+now.label.number,label] <- lambda[label]
                  #cat("lambda not changed, neuron=",label," ",before.lambda[label],"->",lambda[label],"\n")
                  count.reject.lambda = count.reject.lambda + 1
                }
            }
          
      
            ##calculate Energy
            energy = calculate.log.Likelihood(firing.time,amp.time,spike.label,delta,lambda,
              delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)
            
            ##energy = Energy(firing.time,amp.time,spike.label)
            energy.history[(now.step-1)*num.spikes+now.label.number] = energy
          }
      }
    
    write(spike.label.history,file="./result/spike_label_history.txt")
    write(delta.history,file="./result/delta_history.txt")
    write(lambda.history,file="./result/lambda_history.txt")
    write(energy.history,file="./result/energy_history.txt")
    
    
    return(list("spike.label.history"=spike.label.history,
           "delta.history"=delta.history,
           "lambda.history"=lambda.history,
           "energy.history"=energy.history))
  }

#
#Y: observed datas(firing.time, amp.time)
#C: spike label
#params:  parameters(delta,lambda for amplitude attenuation)
#


#
#Y: observed datas(firing.time, amp.time)
#C: spike label
#params:  parameters(delta,lambda for amplitude attenuation)
#

init <- function()
  {


   set.seed(101) 
   #users need to edit these variables on each spike file
    
   ##number of neurons 
   num.neurons=3

   ##number of spikes in one file
   num.spikes = 100
   #To huge
   #num.spikes = 1000
   
   ##spike label
   ##all spikes are labeled with randomized number 1~3
   #spike.label <- round(runif(num.spikes)*2)+1
   ##spikes are labeled by some analysys
   #spike.label <- scan("./simu_data1205/answer.dat")[1:num.spikes]
   #spike.label <- scan("./simu_data1224/SNR3 top3Dkmeans_spikelabel.dat")[1:num.spikes]
   spike.label <- scan("./simu_data01171758/SNR1top10Dkmeans_spikelabel.dat")[1:num.spikes]
   ##Read datas
   ##firing time
   #firing.time <- scan("./simu_data1205/all_firing_time.dat")[1:num.spikes]
   #firing.time <- scan("./simu_data1224/SNR3all_firing_time.dat")[1:num.spikes]
   firing.time <- scan("./simu_data01171758/SNR1all_firing_time.dat")[1:num.spikes]
   #amp.time <- scan("./simu_data1205/max_amp.dat")[1:num.spikes] * 0.01
   #amp.time <- scan("./simu_data1224/SNR3max_amp.dat")[1:num.spikes] * 0.00001
   amp.time <- scan("./simu_data01171758/SNR1max_amp.dat")[1:num.spikes] * 0.00001
   ##firing.time <- scan("./data/all_firing_time.dat")
   #amp.time <- scan("./data/max_amp.dat")

   ##amplitude parameters
   ##amplitude or datas which are attenuated with interspike-interval time
   lambda.min=1
   lambda.max=200#s^-1
   delta.min=0.1
   delta.max=0.9
   ## first lambda settings for each neuron
   first.lambda = (lambda.min + lambda.max)/2
   lambda <- c(first.lambda,first.lambda,first.lambda)
   ## first lambda settings for each neuron
   first.delta = (delta.min + delta.max)/2
   delta <- c(first.delta,first.delta,first.delta)
   

   MCMC.steps = 100

   
   return(list("neurons"=num.neurons,"spikes"=num.spikes,
               "spike.label"=spike.label,
               "firing.time"=firing.time,
               "amp.time"=amp.time,
               "lambda"=lambda,
               "lambda.min"=lambda.min,
               "lambda.max"=lambda.max,
               "delta"=delta,
               "delta.min"=delta.min,
               "delta.max"=delta.max,
               "MCMC.steps"=MCMC.steps
               ))
   
  }

model.amplitude <- function(lapsed.time,amp.value,delta.value,lambda.value)
{
  return(amp.value - amp.value*(1-delta.value*exp(-lambda.value*lapsed.time)))
}

prob.amplitude <- function(lapsed.time,before.amp.value,now.amp.value,delta,lambda)
  {
    ##prior probability density for amplitude
    ##P(amp | isi, delta,lambda)
    ##
    ##spike amplitude attenuates with lapsed time from which former spike fired
    ##
    mean = before.amp.value*(1-delta*exp(-lambda*lapsed.time))
    #sd.amp = abs(now.amp.value - mean)
    sd.amp = abs(mean * 0.1)
    result = dnorm(now.amp.value - mean,mean,sd.amp)
    
    return(result)              
  }
prob.marginal.amplitude.delta <- function(lapsed.time,before.amp.value,now.amp.value,delta,
                                           lambda.min,lambda.max)
  {
    ##prob.amplitude integrated on lambda
    lapse <<- lapsed.time
    before.amp <<- before.amp.value
    now.amp <<- now.amp.value
    delta.value <-  delta
    prob <- function(v)
      {    return(prob.amplitude(lapse,before.amp,now.amp, delta.value, v)) }
    
    result = integrate(prob,lambda.min,lambda.max)
    return(result$value)
  }

prob.marginal.amplitude.lambda <- function(lapsed.time,before.amp.value,now.amp.value,lambda,
                                           delta.min,delta.max)
  {
    ##prob.amplitude integrated on delta 
    lapse <<- lapsed.time
    before.amp <<- before.amp.value
    now.amp <<- now.amp.value
    lambda.value <-  lambda
    prob <- function(v)
      {    return(prob.amplitude(lapse,before.amp,now.amp, v ,lambda.value)) }
    
    result = integrate(prob,delta.min,delta.max)
    return(result$value)
  }

prob.marginal.amplitude.spike.label <- function(lapsed.time,before.amp.value,now.amp.value,
                                                delta.min,delta.max,lambda.min,lambda.max)
  {
    
    lapse <<-lapsed.time
    before.amp <<- before.amp.value
    now.amp <<- now.amp.value
    
    prob <- function(v)
      {    return(prob.amplitude(lapse,before.amp,now.amp, v[1], v[2])) }

    v.lower <-c(delta.min,lambda.min)
    v.upper <- c(delta.max,lambda.max)
    
    result = adapt(2,v.lower,v.upper,functn=prob)
    
  
    
    return(result$value)
  }


prior.prob.delta <- function(delta,delta.min,delta.max)
  {
    ##amplitude parameter
    ##P(delta)
    ##we suppose uniform distribution
  
    return(dunif(delta,min=delta.min,max=delta.max))
    
  }
           
prior.prob.lambda <- function(lambda,lambda.min,lambda.max)
  {
    ##amplitude parameter
    ##P(lambda)
    ##we suppose uniform ditribution

    return(dunif(lambda,min=lambda.min,max=lambda.max))
    
  }
           
prior.prob <- function(delta,delta.min,delta.max,lambda,lambda.min,lambda.max)
  {
    return(prior.prob.delta(delta,delta.min,delta.max)*prior.prob.lambda(lambda,lambda.min,lambda.max))
  }

         
prob.model <- function(label.number,firing.time,before.amp.value,now.amp.value,
                       spike.label,
                       delta,lambda,delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,
                       case=1)
  {
    ##This function is used when calculate likelihood 
    
    ##P(isi,amp | delta,lambda)
    ##my model: PriorProb(isi,amp|delta,lambda,mean,var) == P(amp | delta,lambda)
    ##
    ##my model have no assumption on a interspike-interval histogram shape.
    ##
    
    isi <- interspike.interval(firing.time,spike.label,num.neurons,num.spikes)
    lapsed.time = isi[label.number] 

    
    ##prior.prob = prior.prob.attenuation(isi,amp,delta,lambda) * prob.amplitude(delta,lambda)
    prob.value =
      switch(case,
             prob.amplitude(lapsed.time,before.amp.value,now.amp.value,delta,lambda),
             ##case 2:integrated on lambda
             prob.marginal.amplitude.delta(lapsed.time,before.amp.value,now.amp.value,delta,lambda.min,lambda.max),
             ##case 3:integrated on delta
             prob.marginal.amplitude.lambda(lapsed.time,before.amp.value,now.amp.value,lambda,delta.min,delta.max),
             ##case 4:integrated on delta,lambda
             prob.marginal.amplitude.spike.label(lapsed.time,before.amp.value,now.amp.value,delta.min,delta.max,lambda.min,lambda.max),
             )

    
    return(prob.value)
  }

interspike.interval <- function(firing.time,spike.label,num.neurons,num.spikes)
  {
    ##return  interspike-interval time vector for each neuron with firing.time and spike.label

    ##interspike-interval vector's element isi[i] is  corresponding to  same number of spike.label[i], firing.time[i]
  

    isi <- rep(0,num.spikes)
    
     for(now.neuron in 1:num.neurons)
      {
        for(j in 1:num.spikes)
          {
            before.time=0
            isi.time=0
           
          if(spike.label[j] == now.neuron)
            {
            isi.time = firing.time[j] - before.time
            isi[j] = isi.time
            before.time = firing.time[j]
            }
          }
      }
    ##return as vector
    return(isi)
  }
  
#L(Y,C|params)             
calculate.Likelihood <- function(firing.time,amp.time,spike.label,delta,lambda,num.neurons,num.spikes)
  {
    ##delta,lambda are given as vector
    
    ##Likelihood for all neuron's spikes
    Likelihood = 1
    ##Likelihood for each neuron's spikes
    likelihood <- rep(1,num.neurons) 
    
    for(i in 1:num.neurons)
      {
        before.amp.value = mean(amp.time)##I define the first spike's previous spike's amp 
        for(j in 1:num.spikes)
          {
            before.time=0
            isi.time=0
           
          if(spike.label[j] == num.neurons)
            {
            isi.time = firing.time[j] - before.time
            now.amp.value = amp.time[j]

            
            
            likelihood[i] = likelihood[i] *
              prob.model(i,firing.time,before.amp.value,now.amp.value,spike.label,delta[spike.label[j]],
                         lambda[spike.label[j]],num.neurons,num.spikes)
            before.time = firing.time[j]
            before.amp.value = now.amp.value
            }
       }

        Likelihood = Likelihood * likelihood[i]
      }

    return(Likelihood)
  }
calculate.log.Likelihood <- function(firing.time,amp.time,spike.label,delta,lambda,
                                     delta.min=0.1,delta.max=0.9,lambda.min=10,lambda.max=200,
                                     num.neurons,num.spikes,case=1)
  {
    ##delta,lambda are given as vector
    ##return log.likelihood value
    
    ##log.Likelihood for all neuron's spikes
    log.Likelihood = 0
    ##log.Likelihood for each neuron's spikes
    log.likelihood <- rep(1,num.neurons) 

    ##pre-spike amp value,I suppose first spike's pre-spike 's amp is 1
    before.amp.value <- rep(mean(amp.time),num.neurons)
    
    amp.value <- rep(1,num.neurons)

    
    for(i in 1:num.neurons)
      {
        
        for(j in 1:num.spikes)
          {
            before.time=0
            isi.time=0
            
           
            if(spike.label[j] == i)
              {
                isi.time = firing.time[j] - before.time
            
                amp.value[i] = amp.time[j]
            
                
                
                value = 
                  switch(case,
                         log(prob.model(i,firing.time,before.amp.value[i],amp.value[i],
                                        spike.label,
                                        delta[i],lambda[i],
                                        delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)),
                         ##integrate on lambda
                         log(prob.model(i,firing.time,before.amp.value[i],amp.value[i],
                                        spike.label,
                                        delta[i],lambda[i],
                                        delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=2)),
                         ##integrate on delta
                         log(prob.model(i,firing.time,before.amp.value[i],amp.value[i],
                                        spike.label,
                                        delta[i],lambda[i],
                                        delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=3)),
                         ##integrate on delta,lambda
                         log(prob.model(i,firing.time,before.amp.value[i],amp.value[i],
                                        spike.label,
                                        delta[i],lambda[i],
                                        delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=4))
                         )

                log.likelihood[i] = log.likelihood[i] + value
                    #if(case==1){browser()}
                before.time = firing.time[j]
                before.amp.value[i] = amp.value[i]
              }
          }
        #cat(log.Likelihood)
        log.Likelihood = log.Likelihood + log.likelihood[i]
      }


    return(log.Likelihood[1])
  }


#L(Y,C|params)             
calculate.Likelihood.amp <- function(firing.time,amp.time,spike.label,delta,lambda,num.neurons,num.spikes)
  {
    ##delta,lambda are given as vector
    
    ##Likelihood for all neuron's spikes
    Likelihood = 1
    ##Likelihood for each neuron's spikes
    likelihood <- rep(1,num.neurons) 
    
    for(i in 1:num.neurons)
      {
        for(j in 1:num.spikes)
          {
            before.time=0
            isi.time=0
           
          if(spike.label[j] == num.neurons)
            {
            isi.time = firing.time[j] - before.time
            likelihood[i] = likelihood[i] +
              prob.amplitude(isi.time,amp.time[i],delta[spike.label[j]],lambda[spike.label[j]]) 
            before.time = firing.time[j]
            }
       }

        Likelihood = Likelihood * log.likelihood[i]
      }
    
    return(Likelihood)
  }

#L(Y,C|params)             
log.Likelihood.amp <- function(firing.time,amp.time,spike.label,delta,lambda,num.neurons,num.spikes)
  {
    ##Likelihoo is too small,so we compute log likelihood
    
    ##delta,lambda are given as vector
    
    ##Likelihood for all neuron's spikes
    log.Likelihood = 0
    ##Likelihood for each neuron's spikes
    log.likelihood <- rep(1,num.neurons) 
    
    for(i in 1:num.neurons)
      {
        for(j in 1:num.spikes)
          {
            before.time=0
            isi.time=0
           
          if(spike.label[j] == num.neurons)
            {
            isi.time = firing.time[j] - before.time
            log.likelihood[i] = log.likelihood[i] +
              log(prob.amplitude(isi.time,amp.time[i],delta[spike.label[j]],lambda[spike.label[j]]))
            before.time = firing.time[j]
            }
       }

        log.Likelihood = log.Likelihood + log.likelihood[i]
      }
    return(log.Likelihood)
  }

post.prob <- function(firing.time,amp.time,spike.label,delta,lambda,
                      delta.min,delta.max,lambda.min,lambda.max,
                      num.neurons,num.spikes
                      )
  {
    ##delta,lambda are given as vector
    
    ##P(spike.label,params | firing.time,amp.time) not normolized
    ##posterior probability are only used when we calc the Accept probability
    ##we do not need to calc Z(normolization constant)
     
    ##posterior prob density is proportial to Likelihood * prior.prob    
    Likelihood = calculate.Likelihood(firing.time,amp.time,spike.label,delta,lambda,num.neurons,num.spikes)
    
    return(Likelihood * prior.prob(delta,delta.min,delta.max,lambda,lambda.min,lambda.max))
  }


log.post.prob <- function(firing.time,amp.time,spike.label,delta,lambda,
                          delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=1)
{
  
  log.likelihood = switch(case,
    calculate.log.Likelihood(firing.time,amp.time,spike.label,delta,lambda,
                             delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes),
    calculate.log.Likelihood(firing.time,amp.time,spike.label,delta,lambda,
                             delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=2),
    calculate.log.Likelihood(firing.time,amp.time,spike.label,delta,lambda,
                             delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=3),
     calculate.log.Likelihood(firing.time,amp.time,spike.label,delta,lambda,
                             delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes,case=4)
    )

  log.likelihood = log.likelihood + log( prior.prob(delta,delta.min,delta.max,lambda,lambda.min,lambda.max) )
  
  return(log.likelihood[1])
}

likelihood.amplitude <- function(label.number,firing.time,amp.time,spike.label,delta,lambda,
                                 num.neurons,num.spikes)
  {
    ## log-likelihood function from a single neuron
    ## See the 'Technique for spikesorting' p11,(11)
    ##
    ##L(y,C| delta.lambda)
    likelihood = 0
    isi <- interspike.interval(firing.time,spike.label,num.neurons,num.spikes)    
    likelihood=0
    for(i in 1:num.spikes){
      if(spike.label[i] == label.number){
        amp.value = amp.time[i]
        lapsed.time = isi[i]
        delta.value = delta[label.number]
        lambda.value = lambda[label.number]
        likelihood = likelihood +
              model.amplitude(lapsed.time,amp.value,delta.value,lambda.value)*model.amplitude(lapsed.time,amp.value,delta.value,lambda.value) 
    #    likelihood = model.amplitude(lapsed.time,amp.value,delta.value,lambda.value) 
      }
    }
   
    return(exp(-0.5*likelihood))
  }

generate.piece.wise.lambda <-function(label.number,firing.time,amp.time,spike.label,
                                      delta,lambda,lambda.min,lambda.max,MCMC.steps)
  {
    ##proposed a new lambda from piece-wice approximate Likelihood
    ##
    ##
    partitions = 100
    
    if(MCMC.steps > 100)
      partitions = 15
    lambda.number = spike.label[label.number]
    d.lambda <- seq(lambda.min,lambda.max,by=(lambda.max-lambda.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    lambda.value = lambda[label.number]
    
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] =
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,delta,d.lambda[i],num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.lambda[i] <= lambda.value && d.lambda[i+1] >= lambda.value)
          nearest.lambda = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
    

    i = nearest.lambda
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.lambda[i+1] -d.lambda[i])*(lambda.value-d.lambda[i]) 
      
    N=0 # Normoarized factor
    S <- rep(0,partitions+1) # cmf value in delta[i]
    for(i in 1:partitions)
      {
        N = N + (likelihood.table[i+1] + likelihood.table[i])*(d.lambda[i+1] -d.lambda[i])*0.5
        S[i] = N
      }
    gen = gen/N
    S <- S/N
    u <- runif(1)
    
    for(i in 1:partitions)
      {
        if((S[i]<u) && (u<S[i+1]))
           new.lambda = (d.lambda[i+1] + d.lambda[i])*0.5
      }
    
    return(new.lambda)
  }

propopsed.prob.lambda <- function(label.number,firing.time,amp.time,
                                  spike.label,delta,lambda,lambda.min,lambda.max,MCMC.steps)
  {
    ##proposed a new lambda from piece-wice approximate Likelihood
    ##
    ##
    partitions = 100
    
    if(MCMC.steps > 100)
      partitions = 15
    lambda.number = spike.label[label.number]
    d.lambda <- seq(lambda.min,lambda.max,by=(lambda.max-lambda.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    lambda.value = lambda[label.number]
    
    
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] =
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,delta,d.lambda[i],num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.lambda[i] <= lambda.value && d.lambda[i+1] >= lambda.value)
          nearest.lambda = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
   
     i = nearest.lambda
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.lambda[i+1] -d.lambda[i])*(lambda.value-d.lambda[i]) 
   
    N=0 # Normoarized factor
    S <- rep(0,partitions+1) # cmf value in delta[i]
    for(i in 1:partitions)
      {
        N = N + (likelihood.table[i+1] + likelihood.table[i])*(d.lambda[i+1] -d.lambda[i])*0.5
        S[i] = N
      }
    gen = gen/N
    
    
    return(gen)
  }
  
generate.piece.wise.delta <-
  function(label.number,firing.time,amp.time,spike.label,delta,lambda,delta.min,delta.max,now.steps)
  {
    ##proposed a new delta from piece-wice approximate Likelihood
    ##
    ##delta:vector
    ##
    partitions = 100
    
    if(now.steps > 100)
      partitions = 15

    delta.number = spike.label[label.number]
    d.delta <- seq(delta.min,delta.max,by=(delta.max-delta.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    delta.value = delta[delta.number]
    
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] = 
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,d.delta[i],lambda,num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.delta[i] <= delta.value && d.delta[i+1] >= delta.value)
          nearest.delta = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
    
    
    i = nearest.delta
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.delta[i+1] -d.delta[i])*(delta.value-d.delta[i])
       
    N=0 # Normoarized factor
    S <- rep(0,partitions+1) # cmf value in delta[i]
    for(i in 1:partitions)
      {
        N = N + (likelihood.table[i+1] + likelihood.table[i])*(d.delta[i+1] -d.delta[i])*0.5
        S[i] = N
      }
    gen = gen/N
    S <- S/N
    u <- runif(1)
   
    for(i in 1:partitions)
      {
        if((S[i]<u) && (u<S[i+1]))
           new.delta = (d.delta[i+1] + d.delta[i])*0.5
      }
    
    return(new.delta)
  }

proposed.prob.delta <-
  function(label.number,firing.time,amp.time,spike.label,delta,lambda,delta.min,delta.max,now.steps)
  {
    
    ##proposed a new delta from piece-wice approximate Likelihood
    ##
    ##delta:vector
    ##
    partitions = 100
    
    if(now.steps > 100)
      partitions = 15

    delta.number = spike.label[label.number]
    d.delta <- seq(delta.min,delta.max,by=(delta.max-delta.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    delta.value = delta[delta.number]
   
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] = 
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,d.delta[i],lambda,num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.delta[i] <= delta.value && d.delta[i+1] >= delta.value)
          nearest.delta = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
  
    i=nearest.delta
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.delta[i+1] -d.delta[i])*(delta.value-d.delta[i])
  
    
    N=0 # Normoarized factor
    Integral=0
    S <- rep(0,partitions+1) # cmf value in delta[i]
    for(i in 1:partitions)
      {
        Integral = Integral + (likelihood.table[i+1] + likelihood.table[i])*(d.delta[i+1] -d.delta[i])*0.5
        S[i] = Integral
      }
    gen = gen/Integral
   
    
    
   
    return(gen)
  }

generate.random.walk.delta <-
  function(label.number,spike.label,delta,lambda,delta.min,delta.max,now.steps)
  {
    ##generate new delta from random walk process
    old.delta = delta[spike.label[label.number]]
    ##you choose (search for) a proper value
    delta.sd = (delta.max - delta.min) * 0.1
    
    new.delta = old.delta + rnorm(1,mean=0,sd=delta.sd)
    while(new.delta < delta.min || new.delta > delta.max)
      {
      new.delta = old.delta + rnorm(1,mean=0,sd=delta.sd)
    }

    return(new.delta)
  }

generate.random.walk.lambda <-
    function(label.number,spike.label,delta,lambda,lambda.min,lambda.max,now.steps)
  {
    ##generate new lambda from random walk process
    old.lambda = lambda[spike.label[label.number]]
    ##you choose (search for) a proper value
    lambda.sd = (lambda.max - lambda.min) * 0.1

    
    new.lambda = old.lambda + rnorm(1,mean=0,sd=lambda.sd)
    while((new.lambda < lambda.min)||(new.lambda > lambda.max))
      {
        new.lambda = old.lambda + rnorm(1,mean=0,sd=lambda.sd)
      }
    
    return(new.lambda)
  }

log.proposed.prob.random.walk.delta <-function(label.number,spike.label,
                                                new.delta,old.delta,lambda,delta.min,delta.max,now.steps)
  {
    
    delta.sd = (delta.max - delta.min) * 0.4  
    log.prob = log(pnorm(new.delta - old.delta,mean = 0, sd = delta.sd))
    return(log.prob)
  }

log.proposed.prob.randaom.walk.lambda <- function(label.number,spike.label,delta,
                                                  new.lambda,old.lambda,lambda.min,lambda.max,now.steps)
  {
    
    lambda.sd = (lambda.max - lambda.min) * 0.4  
    log.prob = log(pnorm(new.lambda - old.lambda,mean = 0, sd = lambda.sd))
    return(log.prob)
  }

log.post.marginal.prob.delta <- function(firing.time,amp.time,spike.label,delta,lambda,
                                          delta.min,delta.max,lambda.min,lambda.max,
                                          num.neurons,num.spikes)
  {
    ## posterior marginal robability about delta
    ## P(delta , lambda | firing.time,amp.time,lambda)

    Integral = 0
    ##partitions = 10

     ## four values to four vectors
    #delta.min <- 
    #delta.max <- c(delta.max,delta.max,delta.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    

    ##integrate on spike label
    ## spike label

    ## spike label is a probability variable ,but it is not explicity contained in posterior probability

    
    ##integrate on lambda
    ##for(i in 1:partitions){
    ## d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

  
    ## log.post.value = log.post.prob(firing.time,amp.time,spike.label, delta,d.lambda
    ##   ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
    ## Integral = log.post.value*(lambda.max-lambda.min)/partitions + Integral
    ##}

    ##Integrate on lambda
    Integral  = log.post.prob(firing.time,amp.time,spike.label,delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max,
      num.neurons,num.spikes,case=2)

     
    log.post.value = log.post.prob(firing.time,amp.time,spike.label,delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max,
      num.neurons,num.spikes)

  
    
    post.marginal.prob.value = log.post.value - Integral

    
    
    return(post.marginal.prob.value)
  }

log.post.marginal.prob.lambda <- function(firing.time,amp.time,spike.label,delta,lambda,
                                          delta.min,delta.max,lambda.min,lambda.max,
                                      num.neurons,num.spikes)
  {
    ## posterior marginal probability about delta
    ## P(delta , lambda | firing.time,amp.time,lambda)

    ## posterior marginal probability about delta can not be written by closed form
    ## Pouzat calculates(approximates) it useing piecewise linear function(p2923)
    
    Integral = 0
    partitions = 10

     ## four values to four vectors
    #delta.min <- c(delta.min,delta.min,delta.min)
    #delta.max <- c(delta.max,delta.max,delta.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    


    ##integrate on spike label
    ## spike label

    ## spike label is a probability variable ,but it is not explicity contained in posterior probability

    
    ##integrate on delta
    ##for(i in 1:partitions){
    ## d.delta = delta.min + (delta.max - delta.min)*(i-1)/partitions

    ## post.value = log.post.prob(firing.time,amp.time,spike.label,d.delta,lambda
    ##  ,delta.min,delta.max,lambda.min,lambda.max,
    ## num.neurons,num.spikes)


    ##Integral = post.value/partitions + Integral
    ##}

    ##integrate on delta
    
    Integral = log.post.prob(firing.time,amp.time,spike.label,delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max,
      num.neurons,num.spikes,case=4)

    
    post.value = log.post.prob(firing.time,amp.time,spike.label,delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max,
      num.neurons,num.spikes)

    post.marginal.prob.value = post.value/Integral

    return(post.marginal.prob.value)
  }


post.marginal.prob.spike.label <- function(firing.time,amp.time,spike.label,
                                           delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)
  {
    ##posterior marginal probabilyty about the spike.label
    ##P(spike.label | firing.time,amp.time,delta,lambda)

   
    
    Integral = 0

     ## four values to four vectors
    delta.min <- c(delta.min,delta.min,delta.min)
    delta.max <- c(delta.max,delta.max,delta.max)
    lambda.min <- c(lambda.min,lambda.min,lambda.min)
    lambda.max <- c(lambda.max,lambda.max,lambda.max)
    
    
    partitions = 10 # how many splits the integral region
    ##integrate on delta
    for(i in 1:partitions){
      d.delta = delta.min + (delta.max - delta.min)*(i-1)/partitions
      
      post.value = post.prob(firing.time,amp.time,spike.label, d.delta,lambda
        ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)


      Integral = post.value/partitions + Integral
    }
    ##integrate on lambda
    for(i in 1:partitions){
      d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

      
      post.value = post.prob(firing.time,amp.time,spike.label, delta,d.lambda
        ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      Integral = post.value/partitions + Integral
    }

    post.value = post.prob(firing.time,amp.time,spike.label, delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)
    
    post.marginal.prob.value = post.value/Integral

    return(post.marginal.prob.value)
  }

log.post.marginal.prob.spike.label <- function(firing.time,amp.time,spike.label,
                                               delta,lambda,
                                           delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)
  {
    ##posterior marginal probabilyty about the spike.label
    ##P(spike.label | firing.time,amp.time,delta,lambda)

    ##delta,lambda a vector
    ##return a value
    
    Integral = 0
    ## four values to four vectors
    #delta.min <- c(delta.min,delta.min,delta.min)
    #delta.max <- c(delta.max,delta.max,delta.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    
    #partitions = 10 # how many splits the integral region
    ##integrate on delta
   
    #for(i in 1:partitions){
    #  d.delta = delta.min + (delta.max - delta.min)*(i-1)/partitions
      
    #  log.post.value = log.post.prob(firing.time,amp.time,spike.label, d.delta,lambda
    #    ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
    #  Integral = log.post.value/partitions + Integral
     
    #}
 
    ##integrate on lambda
    #for(i in 1:partitions){
     # d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

      #log.post.value = log.post.prob(firing.time,amp.time,spike.label, delta,d.lambda
       # ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
     # Integral = log.post.value/partitions + Integral
    
   # }

    Integral =log.post.prob(firing.time,amp.time,spike.label, delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes,case=4)


    log.post.value = log.post.prob(firing.time,amp.time,spike.label, delta,lambda
      ,delta.min,delta.max,lambda.min,lambda.max ,num.neurons,num.spikes)

    
    post.marginal.prob.value = log.post.value - Integral

    
    return(post.marginal.prob.value)
  }
    


 Energy <- function(firing.time,amp.time,spike.label)
   {
     likelihood = calculate.Likelihood(firing.time,amp.time,spike.label)

     energy = likelihood * prior.prob(params)
     Energy = log(energy)

     return(Energy)
   }
        

#proposed transision matrix
log.proposed.prob.spike.label.i<- function(label.after,label.before,label.number,firing.time,amp.time,spike.label,
                                       delta,lambda)
  {
    ##compites the proposed probability when (label.number)-th spike label are proposal label(label.after)
    ##see eq(24)
    ##This function returns a NOT-normalized probability value. When we compute MCMC ,we don't need a normalized value. 
    

    ##post.prob() returs too small value to compute,because we use log.post.prob()
    new.spike.label <- spike.label
    new.spike.label[label.number] = label.after

    
    denumerator=0
    for(i in 1:num.neurons)
      {
        if(i != label.number){
          new.spike.label[label.number] = i
          ##delta.value = delta[i]
          ##lambda.value = lambda[i]
          delta.value <- delta
          lambda.value <- lambda
          denumerator = denumerator +
            log.post.prob(firing.time,amp.time,
                          new.spike.label,delta.value,lambda.value
                          ,delta.min,delta.max,lambda.min,lambda.max,
                          num.neurons,num.spikes)
        }
      }
    #denumerator = denumerator[label.number]
    new.spike.label[label.number] = label.after
    ##delta.value = delta[label.after]
    ##lambda.value = lambda[label.after]

    
    numerator =
      log.post.prob(firing.time,amp.time,
                    new.spike.label,delta.value,lambda.value
                    ,delta.min,delta.max,lambda.min,lambda.max,
                    num.neurons,num.spikes)
    #numerator = numerator[label.number]
    
    result = numerator - denumerator

  
    return(result[1])
  }

generate.random.walk.spike.label <- function(spike.label,label.number,num.neurons)
  {
    ##generate (label.number)-th proposal spike label from 
    
    old.label = spike.label[label.number]
    new.label = old.label
    while(old.label == new.label){
      new.label = round(runif(1)*(num.neurons-1)+1)
    }
     
    return(new.label)
  }

log.proposed.random.walk.spike.label <- function(num.neurons)
  {
    return(-log(num.neurons-1))
  }

isi.histgram.trans.kernel <- function(delta.min,delta.max)
  {
    ##This function dosen't use now
    
    ##we generate proposal amplitude params(delta,lambda foe example)
    
    ##See page 2925,Pouzat writes amplitude parameter spike kernel  as inverse-Gamma.
    ##To satisfy with the equillibrium condition, amplitude parameters must be generated from inverse-gamma   
    
    ##inverse-gamma
    mean.log.i <- mean.log.isi(firing.time,spike.label,num.neurons)

    zeta=lambda.min -1 #smaller than lambda.min
    amplitude.gamma = (n.q/2) - 1 
    scale.gamma = (mean.log.i - log.s.q)*(mean.log.i -log.s.q)*n.1/2

    while(delta.min < zeta && delta.max<zeta){
      
      u = pgamma(1,amplitude=gamma,scale=scale.gamma)
      zeta = sqrt(1/u)
    }
      
    result(zeta)
  }


mean.log.isi <- function(firing.time,spike.label,num.neurons)
  {
    ##this function dosen't use now
    isi <- interspike.interval(firing.time,spike.label,num.neurons)
    res <- matrix(nrow=num.neurons,ncol=1)
    for(i in 1:num.neurons)
      {
        v <- isi[[i]]
        v <- log(v)
        v <- sum(v)
        len.v = length(v)

        res[i,] = v / len.v
      }
    return(res)
        
  }
Accept.prob.spike.label <- function(proposed.spike.label,before.spike.label,label.number,
                                    firing.time,amp.time,
                                    delta,lambda,
                                    delta.min,delta.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
  {
    ##see eq(A3)

    ##post.prob() return too small value to compute, we use log.post.prob().
    ##Therefore accept condition are different with an orthodox method.

    log.post.marginal.prob.numerator = log.post.marginal.prob.spike.label(firing.time,amp.time,
      proposed.spike.label,
      delta,lambda,
      delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)

    
    log.post.marginal.prob.denumerator = log.post.marginal.prob.spike.label(firing.time,amp.time,
      before.spike.label,
      delta,lambda,
      delta.min,delta.max,lambda.min,lambda.max,num.neurons,num.spikes)


    log.proposed.prob.numerator =
      log.proposed.random.walk.spike.label(num.neurons)

    log.proposed.prob.denumerator =
      log.proposed.random.walk.spike.label(num.neurons)

    
    

    
    numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
    denumerator= log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

    accept=numerator - denumerator

    
    return(min(0,accept))  
  }


Accept.prob.delta <- function(label.number,
                               proposed.delta,before.delta,
                               firing.time,amp.time,
                               spike.label,
                               delta,delta.min,delta.max,
                               lambda,lambda.min,lambda.max,
                               num.neurons,num.spikes,
                               now.steps)
  {
    
   
    
   log.post.marginal.prob.numerator =
     log.post.marginal.prob.delta(firing.time,amp.time,spike.label,
                                   proposed.delta, lambda,
                                   delta.min,delta.max,lambda.min,lambda.max,
                                   num.neurons,num.spikes)
                        
   log.post.marginal.prob.denumerator =
      log.post.marginal.prob.delta(firing.time,amp.time,spike.label,
                                    before.delta, lambda,
                                    delta.min,delta.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
  
   
   log.proposed.prob.numerator =
      log.proposed.prob.random.walk.delta(label.number,spike.label,
                                           before.delta,proposed.delta,
                                           lambda,delta.min,delta.max,now.steps)    

    
   log.proposed.prob.denumerator =
      log.proposed.prob.random.walk.delta(label.number,spike.label,
                                           proposed.delta,before.delta,
                                           lambda,delta.min,delta.max,now.steps)    


   numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
   denumerator = log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

   accept= numerator- denumerator

   
   return(min(0,accept))  
    
  }


Accept.prob.lambda <- function(label.number,
                               proposed.lambda,before.lambda,
                               firing.time,amp.time,
                               spike.label,
                               delta,delta.min,delta.max,
                               lambda,lambda.min,lambda.max,
                               num.neurons,num.spikes,
                               now.steps)
  {
    
 
    log.post.marginal.prob.numerator =
      log.post.marginal.prob.lambda(firing.time,amp.time,spike.label,delta,
                                proposed.lambda,
                                delta.min,delta.max,lambda.min,lambda.max,
                                num.neurons,num.spikes)

   
    
    log.post.marginal.prob.denumerator =
      log.post.marginal.prob.lambda(firing.time,amp.time,spike.label,delta,
                                    before.lambda,
                                    delta.min,delta.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)

     
    log.proposed.prob.numerator =
      log.proposed.prob.randaom.walk.lambda(label.number,spike.label,
                                            delta,
                                            before.lambda,proposed.lambda,
                                            lambda.min,lambda.max,now.steps)

    
    log.proposed.prob.denumerator =
      log.proposed.prob.randaom.walk.lambda(label.number,spike.label,
                                            delta,
                                            proposed.lambda,proposed.lambda,
                                            lambda.min,lambda.max,now.steps)

    
  

    numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
    denumerator= log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

    accept=numerator - denumerator
      
    return(min(0,accept))  
    
  }
