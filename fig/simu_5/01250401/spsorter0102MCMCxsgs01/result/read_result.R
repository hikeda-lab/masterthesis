library(MASS)
source("check_answer.R")

num.neurons=3
num.spikes=100

delta.history <- matrix(scan(file="delta_history.txt"),ncol=num.neurons)
lambda.history <- matrix(scan(file="lambda_history.txt"),ncol=num.neurons)
spike.label.history <- matrix(scan(file="spike_label_history.txt"),ncol=num.spikes)
init.spike.label <- scan(file="SNR1top10Dkmeans_spikelabel.dat")
##this is posterior prob value
energy.history <- scan(file="energy_history.txt")

answer=scan(file="SNR1answer.dat")

x <- 1:10000
x1 <- 1:10000

cex.val = 1.4
pdf("energy_history.pdf")
par(cex=cex.val)
plot(x,energy.history,type="l",
     xlab="MC step",
     ylab="log-posterior prob")
dev.off()

pdf("delta1_history.pdf")
par(cex=cex.val)
plot(x1,delta.history[,1],type="l",
     xlab="MC step",
     title="neuron1 :delta path")
dev.off()

pdf("delta2_history.pdf")
par(cex=cex.val)
plot(x1,delta.history[,2],type="l",
     xlab="MC step",
     title="neuron2 :delta path")
dev.off()

pdf("delta3_history.pdf")
par(cex=cex.val)
plot(x1,delta.history[,3],type="l",
     xlab="MC step",
     title="neuron3 :delta path")
dev.off()

pdf("lambda1_history.pdf")
par(cex=cex.val)
plot(x1,lambda.history[,1],type="l",
     xlab="MC step",
     title="neuron1 :lambda path")
dev.off()

pdf("lambda2_history.pdf")
par(cex=cex.val)
plot(x1,lambda.history[,2],type="l",
     xlab="MC step",
     title="neuron2 :lambda path")
dev.off()

pdf("lambda3_history.pdf")
par(cex=cex.val)
plot(x1,lambda.history[,3],type="l",
     xlab="MC step",
     title="neuron3 :lambda path")
dev.off()


##probably densities
y = 9000:10000

pdf("hist_delta1.pdf")
par(cex=cex.val)
truehist(delta.history[y,1],
         xlab="delta")
dev.off()

pdf("hist_delta2.pdf")
par(cex=cex.val)
truehist(delta.history[y,2],
         xlab="delta")
dev.off()


pdf("hist_delta3.pdf")
par(cex=cex.val)
truehist(delta.history[y,3],
         xlab="delta")
dev.off()

pdf("hist_lambda1.pdf")
par(cex=cex.val)
truehist(lambda.history[y,1],
         xlab="lambda")
dev.off()

pdf("hist_lambda2.pdf")
par(cex=cex.val)
truehist(lambda.history[y,2],
         xlab="lambda")
dev.off()

pdf("hist_lambda3.pdf")
par(cex=cex.val)
truehist(lambda.history[y,3],
         xlab="lambda")
dev.off()
