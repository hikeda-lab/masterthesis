source("MCMC.R")

##set variables
num.neurons = init()$neurons
num.spikes = init()$spikes
MCMC.steps = init()$MCMC.steps

##observed spike datas
firing.time <- init()$firing.time
amp.time <- init()$amp.time
spike.label <- init()$spike.label
spike.label.history <- matrix(nrow=MCMC.steps,ncol=num.spikes)

##Attenuation parameters
lambda = init()$lambda;lambda.min = init()$lambda.min;lambda.max = init()$lambda.max
##lambda history for each neuron
lambda.history <- matrix(nrow=MCMC.steps,ncol=num.neurons)
deltha = init()$deltha;deltha.min = init()$deltha.min;deltha.max = init()$deltha.max   
##deltha.history for each neuron 
deltha.history <- matrix(nrow=MCMC.steps,ncol=num.neurons)
    

##Energy
energy = 0
energy.history <- matrix(nrow=MCMC.steps,nco=1)
    


##Test functions
lapsed.time=firing.time[10] - firing.time[9]
amp.value=amp.time[10]
deltha.value = deltha[1]
lambda.value=lambda[1]
#model.amp=model.amplitude(lapsed.time,amp.value,deltha.value,lambda.value)
##likelihood <- calculate.Likelihood.amp(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
##likelihood <- calculate.Likelihood.amp(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)



##log.likelihood <- calculate.log.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
#likelihood <- calculate.Likelihood(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)

##likelihood causes errors
#likelihood <- log.Likelihood.amp(firing.time,amp.time,spike.label,deltha,110,num.neurons,num.spikes)

#likelihood.amp = likelihood.amplitude(1,firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)

##generate.piece.wise.deltha causes 6 warnins.
#new.lambda = generate.piece.wise.lambda(1,firing.time,amp.time,spike.label,deltha,lambda,lambda.min,lambda.max,120)
#new.deltha = generate.piece.wise.deltha(1,firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,120)
now.steps=10
label.number = 1
#proposed.prob = proposed.prob.deltha(label.number,firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)

##isi <- interspike.interval(firing.time,spike.label,num.neurons,num.spikes)
#prior.prob.value=prior.prob(deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max)


proposed.lambda.value =
  generate.random.walk.lambda(label.number,spike.label,deltha,lambda,lambda.min,lambda.max,now.steps)
         
 
#marginal.prob.numerator = post.marginal.prob.spike.label(
#  firing.time,amp.timespike.label,deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)

#post.prob.value = post.prob(firing.time,amp.time,spike.label,deltha,lambda,num.neurons,num.spikes)
#log.post.prob.value = log.post.prob(firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,
#  lambda.min,lambda.max,num.neurons,num.spikes)
#post.marginal.prob.value = post.marginal.prob.deltha(firing.time,amp.time,spike.label,deltha,lambda,
#  lambda.min,lambda.max,deltha.min,deltha.max
#  ,num.neurons,num.spikes)

proposed.deltha <- deltha
proposed.deltha[1] = 0.7
log.post.marginal.prob.value =
  log.post.marginal.prob.deltha(firing.time,amp.time,spike.label,
                                proposed.deltha, lambda,
                                deltha.min,deltha.max,lambda.min,lambda.max,
                                num.neurons,num.spikes)




#label.number=14
#accept = Accept.prob.deltha(label.number,
#                   proposed.deltha,before.deltha,
#                   firing.time,amp.time,
#                   spike.label,lambda,lambda.min,lambda.max
#                   ,num.neurons,num.spikes)


## this function returns NULL
#log.post.marginal.prob.value = log.post.marginal.prob.spike.label(firing.time,amp.time,spike.label,
#  deltha.min,deltha.max,lambda.min,lambda.max,num.neurons,num.spikes)

#spike.label[10] = 1
#label.after = 2
#label.before = 1
#label.number = 10
##this function return 3 values vector, I need one value
#trans.value = log.proposed.prob.spike.label.i(label.after,label.before,label.number,firing.time,amp.time,spike.label,
#  deltha,lambda)


#prob.model.value = prob.model(label.number,firing.time,amp.value,deltha,lambda,num.neurons,num.spikes)


#enegy.history <- spikesort.MCMC()

#label=1
#for(i in 1:1){
  
#  before.deltha <- deltha
#  before.deltha[label] <- deltha[label]
#  proposed.deltha <- deltha
#  proposed.deltha[label] =
#    generate.random.walk.deltha(label,spike.label,deltha,lambda,deltha.min,deltha.max,now.steps)
  ##generate.piece.wise.deltha(label,firing.time,amp.time,spike.label,deltha,lambda,deltha.min,deltha.max,now.step)

#  Accept = Accept.prob.deltha(label,proposed.deltha,before.deltha,
#    firing.time,amp.time,spike.label,
#    deltha,deltha.min,deltha.max,lambda,lambda.min,lambda.max,
#    num.neurons,num.spikes,
#    now.steps)
#  cat(Accept,"\n")
#}

#for(i in 1:10){
#  u.lambda = log(runif(1))
#  before.lambda <- lambda
#  before.lambda[label] <- lambda[label]
  ##proposed.lambda =generate.piece.wise.lambda(label,firing.time,amp.time,spike.label,deltha,lambda,lambda.min,lambda.max,now.step)
#  proposed.lambda <- lambda
#  proposed.lambda[label] =
#    generate.random.walk.lambda(label,spike.label,deltha,lambda,lambda.min,lambda.max,now.steps)
          
          
#  Accept = Accept.prob.lambda(label,proposed.lambda,before.lambda,
#    firing.time,amp.time,spike.label,
#    deltha,deltha.min,deltha.max,
#    lambda,lambda.min,lambda.max,
#    num.neurons,num.spikes,
#    now.steps)

#  cat(Accept,"\n")
#}
