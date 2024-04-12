files <- list.files()
for(i in 1:90){
  v <- scan(file = files[i])
  pdf(sprintf("%s%s",files[i],".pdf"))
  plot(1:40,v,type="l",lwd=10,xlab="points",ylab="amplitude")
  dev.off()
}
