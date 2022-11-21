# Practicing with trycatch()


vec1 <- c("5", "5", "5", "fice", "6")
vec2 <- c(0, 0, 0, 0, 0)

for(i in 1:length(vec1)){
  vec2[i] <- bad.numeric(vec1[i])
}

bad.numeric <- function(value){
  v <- as.numeric(value)
  if(is.na(v)){
    stop("Ya funked up brotha, that shit can't be made into a number!")
  }
  return(v) #.. having no return results in a loop stopping error, which is what I need. For some reason having bad.numeric throw an error doesn't stop the execution though?
}



# working on tryCatch now -------------------------------------------------

# for any value that throws an error with bad.numeric() we should see 20 in that position and the loop should complete.
tryCatchTest <- function(vector1, vector2){

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){ vector2[i] <- 20; message('the ERROR condition executed!') },
      warning = function(W){ vector2[i] <- 20; message('the warning condition executed!')}
            )
  }
}

v <- tryCatchTest(vec1, vec2)

vec3 <- c("5", "5", "5", "7", "6")
v2 <- tryCatchTest(vec1, vec3) #verified that the tryCatch block isn't breaking basic funcitonality.

debug(tryCatchTest)
v <- tryCatchTest(vec1, vec2)

vec4 <- as.numeric(vec1) # its vectorized... but why cant I do a single value in my function? after one call of as.numeric on a single value in my vector, I am getting a fully processed output.. not value by value

val1 <- as.numeric(vec1[1L]) #this is working for outputting a single value ... not sure why the funny business is happening in the for loop...

tryCatchTest2 <- function(vector1){

  vector2 <- list()

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){ vector2[i] <- 20; message('the ERROR condition executed!') },
      warning = function(W){ vector2[i] <- 20; message('the warning condition executed!')},
      finally = { vector2[i] <- 20; cat("the finally block executed!!\n") }
    )
  }
  return(vector2)
}


vec5 <- tryCatchTest2(vec1) # assignment doesn't appear to be possible within the error / warning blocks. It is possible within the finally block though... Even in the debugger I never see the list updated for the value that should be added when an error or warning occurs. I am a bit confused about how to properly use the tryCatch() if I can dynamically respond to the error with more code.

debug(tryCatchTest2)
undebug(tryCatchTest2)
?as.vector

tryCatchTest3 <- function(vector1){

  vector2 <- list()

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){  message('the ERROR condition executed!'); return(20) },# hoping the return will assign 20 to the list
      warning = function(W){ message('the warning condition executed!'); return(20) }
    )
  }
  return(vector2)
} # getting closer I think to understanding how to use this tool... I obviously need more practice, as I cannot get it to do what I want... which is simple assignment of value to a for loop specified position in a list of some default value upon catching errors or warnings.

vec6 <- tryCatchTest3(vec1)

