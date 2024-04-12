spwave.template <- function()
  {
    SAMPLES = 40
    num.read = 100 #I read 100 spikes

    
    chosen.files.number <- scan(file="chosen_files_number.txt")

    

    ##shuffle template-files  in the template folder
    ##use 3 templates 
    sp1.num = chosen.files.number[7]
    sp2.num = chosen.files.number[8]
    sp3.num = chosen.files.number[9]

    
    setwd("./select_templates")
    files <- list.files()
    
    sp1 <- matrix(scan(file=files[sp1.num]),ncol=SAMPLES)
    sp2 <- matrix(scan(file=files[sp2.num]),ncol=SAMPLES)
    sp3 <- matrix(scan(file=files[sp3.num]),ncol=SAMPLES)

    setwd("../")

        
    return(list("files"=files,
                "sp1"=sp1,"sp2"=sp2,"sp3"=sp3))
  }


