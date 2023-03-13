# chat gpt help


#I would like the x axis to be a range from 1 to 220000, and the y axis to be 0 or 1, where it will have a point for a given x axis position at 1 on the y axis when the input vector has that value contained within, how can I do that in R?



# Create a vector of indices
indices <- 1:220000

# Create a vector of values
#values <- sample(1:1000, 100, replace = TRUE)

# Create a binary vector indicating presence or absence of values
presence <- ifelse(indices %in% dupeIndex, 1, 0)

# Plot the presence of values at each index
plot(indices, presence, pch = 19, col = "blue", xlab = "Index", ylab = "Presence")


?ifelse


dupeList <- split(dupeIndex,ceiling( seq_along(dupeIndex)/1000 ))
dupeList2 <- split(dupeIndex, 14) # doesn't work
ceiling( seq_along(dupeIndex)/1000 ) # creates a vector of the same length of dupeIndex, which contains numbers that correspond to the bucket each element of dupeIndex will be placed into.
seq_along(dupeIndex)/1000
seq_along(dupeIndex) # creates an index vector, 1 index value for each element in dupeIndex ... i.e, 1:length(dupeIndex)








# learning to inspect objects in memeory to compare for sameness.  --------


.Internal(inspect(i)) # @0x000001ce37a47cc0 13 INTSXP g1c1 [MARK,REF(8)] (len=1, tl=0) 23
a <- i
.Internal(inspect(a))# @0x000001ce37a47cc0 13 INTSXP g1c1 [MARK,REF(8)] (len=1, tl=0) 23
inspect(i) # not a func on its own.. hm ... what does .Internal do?





















































































































