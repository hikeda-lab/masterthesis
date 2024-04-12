
#04-9A_ep013.005.ev_spwaveform_0.dat_template.dat
#03-35Bep008.004.ev_spwaveform_1.dat_template.dat
#0310d_ep031.003.ev_spwaveform_0.dat_template.dat


sp1 <- res.init$template[1,]
sp2 <- res.init$template[2,]
sp3 <- res.init$template[3,]

cex.val=1.4

pdf("template.pdf")
par(cex=cex.val)
plot(1:40,sp1,col=1,type="l",lwd=5,ylim=c(-1400,1000),xlab="Data Samples",ylab="Amplitude")
points(1:40,sp2,col=2,type="l",lwd=5)
points(1:40,sp3,col=3,type="l",lwd=5)
dev.off()

jpeg("template.jpeg")
par(cex=cex.val)
plot(1:40,sp1,col=1,type="l",lwd=5,ylim=c(-1400,1000),xlab="Data Samples",ylab="Amplitude")
points(1:40,sp2,col=2,type="l",lwd=5)
points(1:40,sp3,col=3,type="l",lwd=5)
dev.off()
