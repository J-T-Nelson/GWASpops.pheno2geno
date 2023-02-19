
# testing the use of supplied options to an R script.... Just want to workout the basics


print(arg1)
print(arg2)

if(is.numeric(arg1) && is.numeric(arg2)){
  print(arg1 + arg2)
}


# ok the call: `source('testArgsScript.R', arg1 = 4, arg2 = 6)` doesnt work.. we get the error:
# Error in source("tScript.R", arg1 = 4, arg2 = 6) :
#  unused arguments (arg1 = 4, arg2 = 6)


