#rkcy00ep050.001.ev_spwaveform_0.dat_template.dat
#0318a_ep005.011.ev_spwaveform_0.dat_template.dat
#05-9C_ep031.011.ev_spwaveform_3.dat_template.dat


sp1 <- res.init$template[1,]
sp2 <- res.init$template[2,]
sp3 <- res.init$template[3,]

pdf("template_01151414.pdf")
plot(1:40,sp1,col=1,type="l",lwd=5,ylim=c(-1400,1400),xlab="samples",ylab="Amplitude")
points(1:40,sp2,col=2,type="l",lwd=5)
points(1:40,sp3,col=3,type="l",lwd=5)
dev.off()

jpeg("template_01151414.jpeg")
plot(1:40,sp1,col=1,type="l",lwd=5,ylim=c(-1400,1400),xlab="samples",ylab="Amplitude")
points(1:40,sp2,col=2,type="l",lwd=5)
points(1:40,sp3,col=3,type="l",lwd=5)
dev.off()
