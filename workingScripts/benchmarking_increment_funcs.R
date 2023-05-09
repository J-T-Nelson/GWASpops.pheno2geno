inc <- function(x) {eval.parent(substitute(x <- x + 1))}

inc_but_good <- function(x, inc_val = 1){
  x <- x + inc_val
  }

inc_but_better <- function(x){
  x <- x + 1
  }

# Test One (function)
TestOne <- function() {
  x = 1
  for (i in 1:1e4) {
    inc(x)
  }
}


# Test Two (reference)
TestTwo <- function() {
  x = 1
  for (i in 1:1e4) {
    x <- x + 1
  }
}

# Test Two (reference)
TestThree <- function() {
  x = 1
  for (i in 1:1e4) {
    inc_but_good(x)
  }
}

# Test Two (reference)
TestFour <- function() {
  x = 1
  for (i in 1:1e4) {
    inc_but_better(x)
  }
}

microbenchmark::microbenchmark(
  a = TestOne(),
  b = TestTwo(),
  c = TestThree(),
  d = TestFour()
)

# Unit: microseconds
# expr     min       lq      mean   median       uq     max neval
# a 15216.8 15522.70 19901.039 16013.05 16435.90 69434.0   100
# b    89.9    90.95   107.047    92.10    93.05  1170.4   100
# c  2168.5  2214.95  2327.065  2261.85  2315.60  3769.1   100
# d  1807.3  1849.55  1944.652  1877.10  1930.75  3969.4   100
