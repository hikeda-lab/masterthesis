source("MCMC.R")
times <- system.time(result <- spikesort.MCMC())

write(times,file="times.txt")
