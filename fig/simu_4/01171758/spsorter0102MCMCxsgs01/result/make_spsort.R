burn.in <- 90:100
num.neurons=3

spikesort <-vector(length=num.spikes)
vote <- matrix(nrow=num.spikes,ncol=num.neurons)
for(i in 1:num.spikes)
  for(j in 1:num.neurons)
  vote[i,j] = 0

for(i in burn.in)
  {
    for(j in 1:num.spikes)
      {
        if(spike.label.history[i,j] == 1)
          {vote[j,1] = vote[j,1] + 1}

        if(spike.label.history[i,j] == 2)
          {vote[j,2] = vote[j,2] + 1}
              
        if(spike.label.history[i,j] == 3)
          {vote[j,3] = vote[j,3] + 1}
      }
  }

for( i in 1:num.spikes)  
  spikesort[i] = which(vote[j,] == max(vote[j,]))
