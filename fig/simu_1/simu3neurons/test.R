source("template.R")
source("test_init.R")
source("spikegen.R")

res.init <- init()
res <- spikegen(res.init)
