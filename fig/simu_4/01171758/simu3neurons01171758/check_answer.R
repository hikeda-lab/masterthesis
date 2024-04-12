 check.answer <- function(answer.vec,cluster.vec,num.vec,num.label=2)
  {
    ##count increases count[j,k] +1 when answer.vec[i] == j and cluster.vec[i] == k,
    ##which means i-th spike is labeled k.
    count <- matrix(rep(0,num.label*num.label),ncol=num.label,nrow=num.label)
    

    for(i in 1:num.vec){

      for(j in 1:num.label){
        if(j == answer.vec[i]){
          
          for(k in 1:num.label)
            if(k == cluster.vec[i]){               
              count[j,k] = count[j,k] + 1
            }
        }
      }
    }
              
    return(list("count"=count))
  }
