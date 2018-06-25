

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------- INTERESTING FUNCTIONS ---------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# gets all combinations
expand.grid(bases = DNA_BASES, numbers = 1:2, letter = letters[1:3])


# wrapper for C-style string formatting commands
# %s for string, %i for integer, %f for float.....
verb = c("ho","hai","hanno","sono","vengono")
number = c(2,4,6,1,7)
thing = c("criceti", "dita", "bici", "idiota", "fratelli")
sprintf("%s %i %s", verb, number, thing)
