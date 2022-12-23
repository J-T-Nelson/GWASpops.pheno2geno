# testing new benchmarking function for basic behavior. going to use it to do some optimization
#




dogmode <- function(){
  tic("sleep 1")
  Sys.sleep(20)
  toc()

  tic('sleep2')
  Sys.sleep(40)
  toc()
}

library(tictoc)

dogmode()


