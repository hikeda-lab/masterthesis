plot.enlarge <- function(filename)
{
  sp <- scan(file=filename)

  pdf(sprintf("%s%s%s","enlarge_",filename,".pdf"))
  plot(1:40,sp,type="l",lwd=10,xlab="Samples",ylab="Amplitude")
  dev.off()

}

      

 
