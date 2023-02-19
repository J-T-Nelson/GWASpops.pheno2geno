
args <- commandArgs()

arg1 <- as.numeric(args[1])
arg2 <- as.numeric(args[2])

print(arg1)
print(arg2)

if(is.numeric(arg1) && is.numeric(arg2)){
  print(arg1 + arg2)
}

