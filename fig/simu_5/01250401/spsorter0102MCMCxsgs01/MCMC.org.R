library(adapt)
 

spikesort.MCMC <- function()
  {
    ##set variables
    init.MCMC <- init
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
    lambda.history <- matrix(nrow=MCMC.steps,ncol=num.neurons)
    deltha = init.MCMC$deltha;deltha.min = init.MCMC$deltha.min;deltha.max = init.MCMC$deltha.max   
    ##deltha.history for each neuron 
    deltha.history <- matrix(nrow=MCMC.steps,ncol=num.neurons)

    

    ##Energy
    energy = 0
    energy.history <- matrix(nrow=MCMC.steps,nco=1)
    
    ##for chooose proper parameter(sd) in generate.random.walk
    ##a proper sd generates an accepted deltha about 40%. 
    count.accept.deltha = 0;count.reject.deltha =0
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
             
            cat(now.label.number,"\n")
            before.spike.label <- spike.label
            proposed.spike.label <- spike.label
            proposed.spike.label[now.label.number] <-
              generate.random.walk.spike.label(spike.label, now.label.number , num.neurons)
            
            ##accept or reject the transision
            u = log(runif(1)) #to use log lilklihood,we need a log value by generated from unif distribution
           
            
            Accept.proposed.label = Accept.prob.spike.label(proposed.spike.label,before.spike.label,now.label.number,
                                    firing.time,amp.time,
                                    deltha,lambda,
                                    deltha.min,deltha.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
            #cat(Accept.proposed.label,"\n")
            
            if( Accept.proposed.label > u) # accept the prposed spike label
              {
                spike.label <- proposed.spike.label
                spike.label.history[now.step,] <- spike.label 
              }
            else
              {
                spike.label.history[now.step,] <- spike.label 
              }
            ##When all spike's label are inquired (accept or reject),add spikelabel to history
            spike.label.history[now.step,] <- spike.label
          }

        cat("deltha transition \n")
        ##deltha for each neuron
        for(label in 1:num.neurons){
          u.deltha = log(runif(1)) #use a log likelihood  
          before.deltha <- deltha
          before.deltha[label] <- deltha[label]
          proposed.deltha <- deltha
          proposed.deltha[label] =
            generate.random.walk.deltha(label,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)
          ##generate.piece.wise.deltha(label,firing.time,amp.time,
          ##spike.label,deltha,lambda,deltha.min,deltha.max,now.step)

          Accept = Accept.prob.deltha(label,proposed.deltha,before.deltha,
            firing.time,amp.time,spike.label,
            deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max,
            num.neurons,num.spikes,
            now.steps)

          
          if( Accept > u.deltha)
            {
              deltha <- proposed.deltha
              deltha.history[now.step,label] <- deltha[label]

              count.accept.deltha = count.accept.deltha + 1
            }
          else
            {
              deltha.history[now.step,label] <- deltha[label]
              count.reject.deltha = count.reject.deltha + 1
            }
        }

        cat("lambda transition \n")
        
        ##amplitude parameter lambda for each neuron transision
        for(label in 1:num.neurons){
          ##accept or reject the proposed lambda
          u.lambda = log(runif(1))
          before.lambda <- lambda
          before.lambda[label] <- lambda[label]
          proposed.lambda <- lambda
          ##proposed.lambda =generate.piece.wise.lambda(label,firing.time,amp.time,spike.label,deltha,lambda,lambda.min,lambda.max,now.step)
          proposed.lambda[label] =
            generate.random.walk.lambda(label,spike.label,deltha,lambda,lambda.min,lambda.max,now.steps)
          
        
           Accept = Accept.prob.lambda(label,proposed.lambda,before.lambda,
             firing.time,amp.time,spike.label,
             deltha,deltha.min,deltha.max,
             lambda,lambda.min,lambda.max,
             num.neurons,num.spikes,
             now.steps)
          
          
          if( Accept > u.lambda)
            {
              lambda = proposed.lambda
              lambda.history[now.step,label] <- lambda[label]

              count.accept.lambda = count.accept.lambda + 1
            }
          else
            {
              lambda.history[now.step,label] <- lambda[label]
              count.reject.lambda = count.reject.lambda + 1
            }
        }
      

    ##calculate Energy
    energy = calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,
                             deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)
 
    ##energy = Energy(firing.time,amp.time,spike.label)
    energy.history[now.step] = energy
      
      }
    browser()
    write(spike.label.history,file="./result/spike_label_history.txt")
    write(deltha.history,file="./result/deltha_history.txt")
    write(lambda.history,file="./result/lambda_history.txt")
    write(energy.history,file="./result/energy_history.txt")
    
    #browser()
    return(list("spike.label.history"=spike.label.history,
           "deltha.history"=deltha.history,
           "lambda.history"=lambda.history,
           "energy.history"=energy.history))
  }

#
#Y: observed datas(firing.time, amp.time)
#C: spike label
#params:  parameters(deltha,lambda for amplitude attenuation)
#


#
#Y: observed datas(firing.time, amp.time)
#C: spike label
#params:  parameters(deltha,lambda for amplitude attenuation)
#

init <- function()
  {
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
   spike.label <- scan("./simu_data1205/answer.dat")[1:num.spikes]

   ##Read datas
   ##firing time
   firing.time <- scan("./simu_data1205/all_firing_time.dat")[1:num.spikes]
   amp.time <- scan("./simu_data1205/max_amp.dat")[1:num.spikes] * 0.01
   #firing.time <- scan("./data/all_firing_time.dat")
   #amp.time <- scan("./data/max_amp.dat")

   ##amplitude parameters
   ##amplitude or datas which are attenuated with interspike-interval time
   lambda.min=1
   lambda.max=200#s^-1
   deltha.min=0.1
   deltha.max=0.9
   ## first lambda settings for each neuron
   first.lambda = (lambda.min + lambda.max)/2
   lambda <- c(first.lambda,first.lambda,first.lambda)
   ## first lambda settings for each neuron
   first.deltha = (deltha.min + deltha.max)/2
   deltha <- c(first.deltha,first.deltha,first.deltha)
   

   MCMC.steps = 100

   
   return(list("neurons"=num.neurons,"spikes"=num.spikes,
               "spike.label"=spike.label,
               "firing.time"=firing.time,
               "amp.time"=amp.time,
               "lambda"=lambda,
               "lambda.min"=lambda.min,
               "lambda.max"=lambda.max,
               "deltha"=deltha,
               "deltha.min"=deltha.min,
               "deltha.max"=deltha.max,
               "MCMC.steps"=MCMC.steps
               ))
   
  }

model.amplitude <- function(lapsed.time,amp.value,deltha.value,lambda.value)
{
  return(amp.value - amp.value*(1-deltha.value*exp(-lambda.value*lapsed.time)))
}

prob.amplitude <- function(lapsed.time,amp.value,deltha,lambda)
  {
    ##prior probability density for amplitude
    ##P(amp | isi, deltha,lambda)
    ##
    ##spike amplitude attenuates with lapsed time from which former spike fired
    ##
    mean = amp.value*(1-deltha*exp(-lambda*lapsed.time))
    var = 1
    result = dnorm(1,mean,var)
    return(result)              
  }
prob.marginal.amplitude.deltha <- function(lapsed.time,amp.value,deltha,
                                           lambda.min,lambda.max)
  {
    ##prob.amplitude integrated on lambda
    lapse <<- lapsed.time
    amp <<- amp.value
    deltha.value <-  deltha
    prob <- function(v)
      {    return(prob.amplitude(lapse,amp, deltha.value, v)) }
    
    result = integrate(prob,lambda.min,lambda.max)
    return(result$value)
  }

prob.marginal.amplitude.lambda <- function(lapsed.time,amp.value,lambda,
                                           deltha.min,deltha.max)
  {
    ##prob.amplitude integrated on deltha 
    lapse <<- lapsed.time
    amp <<- amp.value
    lambda.value <-  lambda
    prob <- function(v)
      {    return(prob.amplitude(lapse,amp, v ,lambda.value)) }
    
    result = integrate(prob,deltha.min,deltha.max)
    return(result$value)
  }

prob.marginal.amplitude.spike.label <- function(lapsed.time,amp.value,
                                                deltha.min,deltha.max,lambda.min,lambda.max)
  {
    
    lapse <<-lapsed.time
    amp <<- amp.value
  
    prob <- function(v)
      {    return(prob.amplitude(lapse,amp, v[1], v[2])) }

    v.lower <-c(deltha.min,lambda.min)
    v.upper <- c(deltha.max,lambda.max)
    
    result = adapt(2,v.lower,v.upper,functn=prob)
   
    return(result$value)
  }


prior.prob.deltha <- function(deltha,deltha.min,deltha.max)
  {
    ##amplitude parameter
    ##P(deltha)
    ##we suppose uniform distribution
  
    return(dunif(deltha,min=deltha.min,max=deltha.max))
    
  }
           
prior.prob.lambda <- function(lambda,lambda.min,lambda.max)
  {
    ##amplitude parameter
    ##P(lambda)
    ##we suppose uniform ditribution

    return(dunif(lambda,min=lambda.min,max=lambda.max))
    
  }
           
prior.prob <- function(deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max)
  {
    return(prior.prob.deltha(deltha,deltha.min,deltha.max)*prior.prob.lambda(lambda,lambda.min,lambda.max))
  }

         
prob.model <- function(label.number,firing.time,amp.value,spike.label,
                       deltha,lambda,deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,
                       case=1)
  {
    ##This function is used when calculate likelihood 
    
    ##P(isi,amp | deltha,lambda)
    ##my model: PriorProb(isi,amp|deltha,lambda,mean,var) == P(amp | deltha,lambda)
    ##
    ##my model have no assumption on a interspike-interval histogram shape.
    ##
    
    isi <- interspike.interval(firing.time,spike.label,num.neurons,num.spikes)
    lapsed.time = isi[label.number] 
    
    ##prior.prob = prior.prob.attenuation(isi,amp,deltha,lambda) * prob.amplitude(deltha,lambda)
    prob.value =
      switch(case,
             prob.amplitude(lapsed.time,amp.value,deltha,lambda),
             ##case 2:integrated on lambda
             prob.marginal.amplitude.deltha(lapsed.time,amp.value,deltha,lambda.min,lambda.max),
             ##case 3:integrated on deltha
             prob.marginal.amplitude.lambda(lapsed.time,amp.value,lambda,deltha.min,deltha.max),
             ##case 4:integrated on deltha,lambda
             prob.marginal.amplitude.spike.label(lapsed.time,amp.value,deltha.min,deltha.max,lambda.min,lambda.max),
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
calculate.Likelihood <- function(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
  {
    ##deltha,lambda are given as vector
    
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
            amp.value = amp.time[j]
            likelihood[i] = likelihood[i] *
              prob.model(i,firing.time,amp.time[i],spike.label,deltha[spike.label[j]],
                         lambda[spike.label[j]],num.neurons,num.spikes)
            before.time = firing.time[j]
            }
       }

        Likelihood = Likelihood * likelihood[i]
      }

    return(Likelihood)
  }
calculate.log.Likelihood <- function(firing.time,amp.time,spike.label,deltha,lambda,
                                     deltha.min=0.1,deltha.max=0.9,lambda.min=10,lambda.max=200,
                                     num.neurons,num.spikes,case=1)
  {
    ##deltha,lambda are given as vector
    ##return log.likelihood value
    
    ##log.Likelihood for all neuron's spikes
    log.Likelihood = 0
    ##log.Likelihood for each neuron's spikes
    log.likelihood <- rep(1,num.neurons) 

    ##pre-spike amp value,I suppose first spike's pre-spike 's amp is 1
    before.amp.value <- rep(1,num.neurons)
    
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
            
                
                log.likelihood[i] = log.likelihood[i] +
                  switch(case,
                         log(prob.model(i,firing.time,before.amp.value[i],spike.label,
                                        deltha[i],lambda[i],
                                        deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)),
                         ##integrate on lambda
                         log(prob.model(i,firing.time,before.amp.value[i],spike.label,
                                        deltha[i],lambda[i],
                                        deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=2)),
                         ##integrate on deltha
                         log(prob.model(i,firing.time,before.amp.value[i],spike.label,
                                        deltha[i],lambda[i],
                                        deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=3)),
                         ##integrate on deltha,lambda
                         log(prob.model(i,firing.time,before.amp.value[i],spike.label,
                                        deltha[i],lambda[i],
                                        deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=4))
                         )
                
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
calculate.Likelihood.amp <- function(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
  {
    ##deltha,lambda are given as vector
    
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
              prob.amplitude(isi.time,amp.time[i],deltha[spike.label[j]],lambda[spike.label[j]]) 
            before.time = firing.time[j]
            }
       }

        Likelihood = Likelihood * log.likelihood[i]
      }
    
    return(Likelihood)
  }

#L(Y,C|params)             
log.Likelihood.amp <- function(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
  {
    ##Likelihoo is too small,so we compute log likelihood
    
    ##deltha,lambda are given as vector
    
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
              log(prob.amplitude(isi.time,amp.time[i],deltha[spike.label[j]],lambda[spike.label[j]]))
            before.time = firing.time[j]
            }
       }

        log.Likelihood = log.Likelihood + log.likelihood[i]
      }
    return(log.Likelihood)
  }

post.prob <- function(firing.time,amp.time,spike.label,deltha,lambda,
                      deltha.min,deltha.max,lambda.min,lambda.max,
                      num.neurons,num.spikes
                      )
  {
    ##deltha,lambda are given as vector
    
    ##P(spike.label,params | firing.time,amp.time) not normolized
    ##posterior probability are only used when we calc the Accept probability
    ##we do not need to calc Z(normolization constant)
     
    ##posterior prob density is proportial to Likelihood * prior.prob    
    Likelihood = calculate.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
    
    return(Likelihood * prior.prob(deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max))
  }


log.post.prob <- function(firing.time,amp.time,spike.label,deltha,lambda,
                          deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=1)
{
  
  log.likelihood = switch(case,
    calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,
                             deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes),
    calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,
                             deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=2),
    calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,
                             deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=3),
     calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,
                             deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes,case=4)
    )

  log.likelihood = log.likelihood + log( prior.prob(deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max) )
  
  return(log.likelihood[1])
}

likelihood.amplitude <- function(label.number,firing.time,amp.time,spike.label,deltha,lambda,
                                 num.neurons,num.spikes)
  {
    ## log-likelihood function from a single neuron
    ## See the 'Technique for spikesorting' p11,(11)
    ##
    ##L(y,C| deltha.lambda)
    likelihood = 0
    isi <- interspike.interval(firing.time,spike.label,num.neurons,num.spikes)    
    likelihood=0
    for(i in 1:num.spikes){
      if(spike.label[i] == label.number){
        amp.value = amp.time[i]
        lapsed.time = isi[i]
        deltha.value = deltha[label.number]
        lambda.value = lambda[label.number]
        likelihood = likelihood +
              model.amplitude(lapsed.time,amp.value,deltha.value,lambda.value)*model.amplitude(lapsed.time,amp.value,deltha.value,lambda.value) 
    #    likelihood = model.amplitude(lapsed.time,amp.value,deltha.value,lambda.value) 
      }
    }
   
    return(exp(-0.5*likelihood))
  }

generate.piece.wise.lambda <-function(label.number,firing.time,amp.time,spike.label,
                                      deltha,lambda,lambda.min,lambda.max,MCMC.steps)
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
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,deltha,d.lambda[i],num.neurons,num.spikes)
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
    S <- rep(0,partitions+1) # cmf value in deltha[i]
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
                                  spike.label,deltha,lambda,lambda.min,lambda.max,MCMC.steps)
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
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,deltha,d.lambda[i],num.neurons,num.spikes)
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
    S <- rep(0,partitions+1) # cmf value in deltha[i]
    for(i in 1:partitions)
      {
        N = N + (likelihood.table[i+1] + likelihood.table[i])*(d.lambda[i+1] -d.lambda[i])*0.5
        S[i] = N
      }
    gen = gen/N
    
    
    return(gen)
  }
  
generate.piece.wise.deltha <-
  function(label.number,firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)
  {
    ##proposed a new deltha from piece-wice approximate Likelihood
    ##
    ##deltha:vector
    ##
    partitions = 100
    
    if(now.steps > 100)
      partitions = 15

    deltha.number = spike.label[label.number]
    d.deltha <- seq(deltha.min,deltha.max,by=(deltha.max-deltha.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    deltha.value = deltha[deltha.number]
    
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] = 
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,d.deltha[i],lambda,num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.deltha[i] <= deltha.value && d.deltha[i+1] >= deltha.value)
          nearest.deltha = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
    
    
    i = nearest.deltha
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.deltha[i+1] -d.deltha[i])*(deltha.value-d.deltha[i])
       
    N=0 # Normoarized factor
    S <- rep(0,partitions+1) # cmf value in deltha[i]
    for(i in 1:partitions)
      {
        N = N + (likelihood.table[i+1] + likelihood.table[i])*(d.deltha[i+1] -d.deltha[i])*0.5
        S[i] = N
      }
    gen = gen/N
    S <- S/N
    u <- runif(1)
   
    for(i in 1:partitions)
      {
        if((S[i]<u) && (u<S[i+1]))
           new.deltha = (d.deltha[i+1] + d.deltha[i])*0.5
      }
    
    return(new.deltha)
  }

proposed.prob.deltha <-
  function(label.number,firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)
  {
    
    ##proposed a new deltha from piece-wice approximate Likelihood
    ##
    ##deltha:vector
    ##
    partitions = 100
    
    if(now.steps > 100)
      partitions = 15

    deltha.number = spike.label[label.number]
    d.deltha <- seq(deltha.min,deltha.max,by=(deltha.max-deltha.min)/(partitions))
    likelihood.table <- rep(1,partitions+1)
    deltha.value = deltha[deltha.number]
   
    ##set tables
    for(i in 1:partitions+1)
      {
        likelihood.table[i] = 
          likelihood.amplitude(label.number,firing.time,amp.time,spike.label,d.deltha[i],lambda,num.neurons,num.spikes)
      }
     for(i in 1:partitions)
      {
        if(d.deltha[i] <= deltha.value && d.deltha[i+1] >= deltha.value)
          nearest.deltha = i
      }
    ##piece-wice approximate
    ##gen is may not necessary
     gen = 0
  
    i=nearest.deltha
    gen <-  likelihood.table[i] +
      (likelihood.table[i+1] - likelihood.table[i])/(d.deltha[i+1] -d.deltha[i])*(deltha.value-d.deltha[i])
  
    
    N=0 # Normoarized factor
    Integral=0
    S <- rep(0,partitions+1) # cmf value in deltha[i]
    for(i in 1:partitions)
      {
        Integral = Integral + (likelihood.table[i+1] + likelihood.table[i])*(d.deltha[i+1] -d.deltha[i])*0.5
        S[i] = Integral
      }
    gen = gen/Integral
   
    
    
   
    return(gen)
  }

generate.random.walk.deltha <-
  function(label.number,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)
  {
    ##generate new deltha from random walk process
    old.deltha = deltha[spike.label[label.number]]
    ##you choose (search for) a proper value
    deltha.sd = (deltha.max - deltha.min) * 0.4
    
    new.deltha = old.deltha + rnorm(1,mean=0,sd=deltha.sd)
    while(new.deltha < deltha.min || new.deltha > deltha.max)
      {
      new.deltha = old.deltha + rnorm(1,mean=0,sd=deltha.sd)
    }

    return(new.deltha)
  }

generate.random.walk.lambda <-
    function(label.number,spike.label,deltha,lambda,lambda.min,lambda.max,now.steps)
  {
    ##generate new lambda from random walk process
    old.lambda = lambda[spike.label[label.number]]
    ##you choose (search for) a proper value
    lambda.sd = (lambda.max - lambda.min) * 0.4

    
    new.lambda = old.lambda + rnorm(1,mean=0,sd=lambda.sd)
    while((new.lambda < lambda.min)||(new.lambda > lambda.max))
      {
        new.lambda = old.lambda + rnorm(1,mean=0,sd=lambda.sd)
      }
    
    return(new.lambda)
  }

log.proposed.prob.random.walk.deltha <-function(label.number,spike.label,
                                                new.deltha,old.deltha,lambda,deltha.min,deltha.max,now.steps)
  {
    
    deltha.sd = (deltha.max - deltha.min) * 0.4  
    log.prob = log(pnorm(new.deltha - old.deltha,mean = 0, sd = deltha.sd))
    return(log.prob)
  }

log.proposed.prob.randaom.walk.lambda <- function(label.number,spike.label,deltha,
                                                  new.lambda,old.lambda,lambda.min,lambda.max,now.steps)
  {
    
    lambda.sd = (lambda.max - lambda.min) * 0.4  
    log.prob = log(pnorm(new.lambda - old.lambda,mean = 0, sd = lambda.sd))
    return(log.prob)
  }

log.post.marginal.prob.deltha <- function(firing.time,amp.time,spike.label,deltha,lambda,
                                          deltha.min,deltha.max,lambda.min,lambda.max,
                                          num.neurons,num.spikes)
  {
    ## posterior marginal robability about deltha
    ## P(deltha , lambda | firing.time,amp.time,lambda)

    Integral = 0
    ##partitions = 10

     ## four values to four vectors
    #deltha.min <- 
    #deltha.max <- c(deltha.max,deltha.max,deltha.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    

    ##integrate on spike label
    ## spike label

    ## spike label is a probability variable ,but it is not explicity contained in posterior probability

    
    ##integrate on lambda
    ##for(i in 1:partitions){
    ## d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

  
    ## log.post.value = log.post.prob(firing.time,amp.time,spike.label, deltha,d.lambda
    ##   ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
    ## Integral = log.post.value*(lambda.max-lambda.min)/partitions + Integral
    ##}

    ##Integrate on lambda
    Integral  = log.post.prob(firing.time,amp.time,spike.label,deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max,
      num.neurons,num.spikes,case=2)

     
    log.post.value = log.post.prob(firing.time,amp.time,spike.label,deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max,
      num.neurons,num.spikes)

  
    
    post.marginal.prob.value = log.post.value - Integral

    
    #browser()
    return(post.marginal.prob.value)
  }

log.post.marginal.prob.lambda <- function(firing.time,amp.time,spike.label,deltha,lambda,
                                          deltha.min,deltha.max,lambda.min,lambda.max,
                                      num.neurons,num.spikes)
  {
    ## posterior marginal probability about deltha
    ## P(deltha , lambda | firing.time,amp.time,lambda)

    ## posterior marginal probability about deltha can not be written by closed form
    ## Pouzat calculates(approximates) it useing piecewise linear function(p2923)
    
    Integral = 0
    partitions = 10

     ## four values to four vectors
    #deltha.min <- c(deltha.min,deltha.min,deltha.min)
    #deltha.max <- c(deltha.max,deltha.max,deltha.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    


    ##integrate on spike label
    ## spike label

    ## spike label is a probability variable ,but it is not explicity contained in posterior probability

    
    ##integrate on deltha
    ##for(i in 1:partitions){
    ## d.deltha = deltha.min + (deltha.max - deltha.min)*(i-1)/partitions

    ## post.value = log.post.prob(firing.time,amp.time,spike.label,d.deltha,lambda
    ##  ,deltha.min,deltha.max,lambda.min,lambda.max,
    ## num.neurons,num.spikes)


    ##Integral = post.value/partitions + Integral
    ##}

    ##integrate on deltha
    
    Integral = log.post.prob(firing.time,amp.time,spike.label,deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max,
      num.neurons,num.spikes,case=4)

    
    post.value = log.post.prob(firing.time,amp.time,spike.label,deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max,
      num.neurons,num.spikes)

    post.marginal.prob.value = post.value/Integral

    return(post.marginal.prob.value)
  }


post.marginal.prob.spike.label <- function(firing.time,amp.time,spike.label,
                                           deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)
  {
    ##posterior marginal probabilyty about the spike.label
    ##P(spike.label | firing.time,amp.time,deltha,lambda)

   
    
    Integral = 0

     ## four values to four vectors
    deltha.min <- c(deltha.min,deltha.min,deltha.min)
    deltha.max <- c(deltha.max,deltha.max,deltha.max)
    lambda.min <- c(lambda.min,lambda.min,lambda.min)
    lambda.max <- c(lambda.max,lambda.max,lambda.max)
    
    
    partitions = 10 # how many splits the integral region
    ##integrate on deltha
    for(i in 1:partitions){
      d.deltha = deltha.min + (deltha.max - deltha.min)*(i-1)/partitions
      
      post.value = post.prob(firing.time,amp.time,spike.label, d.deltha,lambda
        ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)


      Integral = post.value/partitions + Integral
    }
    ##integrate on lambda
    for(i in 1:partitions){
      d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

      
      post.value = post.prob(firing.time,amp.time,spike.label, deltha,d.lambda
        ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      Integral = post.value/partitions + Integral
    }

    post.value = post.prob(firing.time,amp.time,spike.label, deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)
    
    post.marginal.prob.value = post.value/Integral

    return(post.marginal.prob.value)
  }

log.post.marginal.prob.spike.label <- function(firing.time,amp.time,spike.label,
                                               deltha,lambda,
                                           deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)
  {
    ##posterior marginal probabilyty about the spike.label
    ##P(spike.label | firing.time,amp.time,deltha,lambda)

    ##deltha,lambda a vector
    ##return a value
    
    Integral = 0
    ## four values to four vectors
    #deltha.min <- c(deltha.min,deltha.min,deltha.min)
    #deltha.max <- c(deltha.max,deltha.max,deltha.max)
    #lambda.min <- c(lambda.min,lambda.min,lambda.min)
    #lambda.max <- c(lambda.max,lambda.max,lambda.max)
    
    #partitions = 10 # how many splits the integral region
    ##integrate on deltha
   
    #for(i in 1:partitions){
    #  d.deltha = deltha.min + (deltha.max - deltha.min)*(i-1)/partitions
      
    #  log.post.value = log.post.prob(firing.time,amp.time,spike.label, d.deltha,lambda
    #    ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
    #  Integral = log.post.value/partitions + Integral
     
    #}
 
    ##integrate on lambda
    #for(i in 1:partitions){
     # d.lambda = lambda.min + (lambda.max - lambda.min)*(i-1)/partitions

      #log.post.value = log.post.prob(firing.time,amp.time,spike.label, deltha,d.lambda
       # ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)

      
     # Integral = log.post.value/partitions + Integral
    
   # }

    Integral =log.post.prob(firing.time,amp.time,spike.label, deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes,case=4)


    log.post.value = log.post.prob(firing.time,amp.time,spike.label, deltha,lambda
      ,deltha.min,deltha.max,lambda.min,lambda.max ,num.neurons,num.spikes)

    
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
                                       deltha,lambda)
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
          ##deltha.value = deltha[i]
          ##lambda.value = lambda[i]
          deltha.value <- deltha
          lambda.value <- lambda
          denumerator = denumerator +
            log.post.prob(firing.time,amp.time,
                          new.spike.label,deltha.value,lambda.value
                          ,deltha.min,deltha.max,lambda.min,lambda.max,
                          num.neurons,num.spikes)
        }
      }
    #denumerator = denumerator[label.number]
    new.spike.label[label.number] = label.after
    ##deltha.value = deltha[label.after]
    ##lambda.value = lambda[label.after]

    
    numerator =
      log.post.prob(firing.time,amp.time,
                    new.spike.label,deltha.value,lambda.value
                    ,deltha.min,deltha.max,lambda.min,lambda.max,
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

isi.histgram.trans.kernel <- function(deltha.min,deltha.max)
  {
    ##This function dosen't use now
    
    ##we generate proposal amplitude params(deltha,lambda foe example)
    
    ##See page 2925,Pouzat writes amplitude parameter spike kernel  as inverse-Gamma.
    ##To satisfy with the equillibrium condition, amplitude parameters must be generated from inverse-gamma   
    
    ##inverse-gamma
    mean.log.i <- mean.log.isi(firing.time,spike.label,num.neurons)

    zeta=lambda.min -1 #smaller than lambda.min
    amplitude.gamma = (n.q/2) - 1 
    scale.gamma = (mean.log.i - log.s.q)*(mean.log.i -log.s.q)*n.1/2

    while(deltha.min < zeta && deltha.max<zeta){
      
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
                                    deltha,lambda,
                                    deltha.min,deltha.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
  {
    ##see eq(A3)

    ##post.prob() return too small value to compute, we use log.post.prob().
    ##Therefore accept condition are different with an orthodox method.

    log.post.marginal.prob.numerator = log.post.marginal.prob.spike.label(firing.time,amp.time,
      proposed.spike.label,
      deltha,lambda,
      deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)

    
    log.post.marginal.prob.denumerator = log.post.marginal.prob.spike.label(firing.time,amp.time,
      before.spike.label,
      deltha,lambda,
      deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)


    log.proposed.prob.numerator =
      log.proposed.random.walk.spike.label(num.neurons)

    log.proposed.prob.denumerator =
      log.proposed.random.walk.spike.label(num.neurons)

    
    

    
    numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
    denumerator= log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

    accept=numerator - denumerator
    
    return(min(0,accept))  
  }


Accept.prob.deltha <- function(label.number,
                               proposed.deltha,before.deltha,
                               firing.time,amp.time,
                               spike.label,
                               deltha,deltha.min,deltha.max,
                               lambda,lambda.min,lambda.max,
                               num.neurons,num.spikes,
                               now.steps)
  {
    
   
    
   log.post.marginal.prob.numerator =
     log.post.marginal.prob.deltha(firing.time,amp.time,spike.label,
                                   proposed.deltha, lambda,
                                   deltha.min,deltha.max,lambda.min,lambda.max,
                                   num.neurons,num.spikes)
                        
   log.post.marginal.prob.denumerator =
      log.post.marginal.prob.deltha(firing.time,amp.time,spike.label,
                                    before.deltha, lambda,
                                    deltha.min,deltha.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)
  
   
   log.proposed.prob.numerator =
      log.proposed.prob.random.walk.deltha(label.number,spike.label,
                                           before.deltha,proposed.deltha,
                                           lambda,deltha.min,deltha.max,now.steps)    

    
   log.proposed.prob.denumerator =
      log.proposed.prob.random.walk.deltha(label.number,spike.label,
                                           proposed.deltha,before.deltha,
                                           lambda,deltha.min,deltha.max,now.steps)    


   numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
   denumerator = log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

   accept= numerator- denumerator

   
   return(min(0,accept))  
    
  }


Accept.prob.lambda <- function(label.number,
                               proposed.lambda,before.lambda,
                               firing.time,amp.time,
                               spike.label,
                               deltha,deltha.min,deltha.max,
                               lambda,lambda.min,lambda.max,
                               num.neurons,num.spikes,
                               now.steps)
  {
    
 
    log.post.marginal.prob.numerator =
      log.post.marginal.prob.lambda(firing.time,amp.time,spike.label,deltha,
                                proposed.lambda,
                                deltha.min,deltha.max,lambda.min,lambda.max,
                                num.neurons,num.spikes)

   
    
    log.post.marginal.prob.denumerator =
      log.post.marginal.prob.lambda(firing.time,amp.time,spike.label,deltha,
                                    before.lambda,
                                    deltha.min,deltha.max,lambda.min,lambda.max,
                                    num.neurons,num.spikes)

     
    log.proposed.prob.numerator =
      log.proposed.prob.randaom.walk.lambda(label.number,spike.label,
                                            deltha,
                                            before.lambda,proposed.lambda,
                                            lambda.min,lambda.max,now.steps)

    
    log.proposed.prob.denumerator =
      log.proposed.prob.randaom.walk.lambda(label.number,spike.label,
                                            deltha,
                                            proposed.lambda,proposed.lambda,
                                            lambda.min,lambda.max,now.steps)

    
  

    numerator=  log.post.marginal.prob.numerator  + log.proposed.prob.numerator
    denumerator= log.post.marginal.prob.denumerator + log.proposed.prob.denumerator

    accept=numerator - denumerator
      
    return(min(0,accept))  
    
  }
