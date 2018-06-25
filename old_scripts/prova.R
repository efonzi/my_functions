#-------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------- 170612 - redirect STDOUT ---------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

m = system('ls /Users', intern = T)

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------- 170430 - hclust() test ---------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
x = c(1,1.5,2,2.5,2.5,4,4,5,6,6,7)
y = c(1,2,1,3,4,1,4,5,4,6,1)
m = matrix(c(x,y), ncol=2)
rownames(m) = letters[1:11]

plot(m)

d = dist(m, method="euclidean")

fit <- hclust(d, method="ward.D")
plot(fit, main="ward.D")

fit <- hclust(d, method="ward.D2")
plot(fit, main="ward.D2")

fit <- hclust(d, method="single")
plot(fit, main="single")

fit <- hclust(d, method="complete")
plot(fit, main="complete")

fit <- hclust(d, method="average") # UPGMA
plot(fit, main="average")

fit <- hclust(d, method="mcquitty") # WPGMA
plot(fit, main="mcquitty")

fit <- hclust(d, method="median") # WPGMC
plot(fit, main="median")

fit <- hclust(d, method="centroid") # UPGMC
plot(fit, main="centroid")






#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------- dgv sources stats ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# load and edit file
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/dgv")
dgv_papers = read.csv("dgv_papers.csv", header=T, sep=",", dec=".", quote = "'", stringsAsFactors = F, check.names = F)
colnames(dgv_papers)[2] = "var"
# transform counts to proportions and barplot
dgv_papers$proportions = dgv_papers$var/sum(dgv_papers$var)
prop = dgv_papers$proportions
jpeg(filename="barplot_DGVvariant_study.jpeg", width=9, height=6, units="in", res=175, quality=100)
par(mar=c(5,10,0,1)) # set margins
barplot(prop, names.arg = dgv_papers$study, horiz=T, las=1, cex.names = 0.5, xlab="Proportion of DGV variants per study")
dev.off()
# proportion of first 7 sources together? ---> 0.987
sum(prop[1:7])/sum(prop)




#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------- PARALLEL PROGRAMMING ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
library(parallel)
cores = detectCores() -1 # detect number of total cores - 1
cluster = makeCluster(cores) # build cluster on "cores" number of cores
# syntax to call a function in parallel is:
# parLapply(cluster, levels, function())
stopCluster(cluster) # close cluster

prova <- function(exponent){
  base = 3
  base^exponent}


parLapply(cluster, 2:4, prova)

stopCluster(cl)




#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------ OTHER -----------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# in case of error such as: "Error in plot.new() : figure margins too large"
# to check the current margins, usually "[1] 5.1 4.1 4.1 2.1"
par("mar")
# to change the margins as desired
# the number indicate margins in clockwise order, startin from lower margin
par(mar=c(3,3,2,1))

# binomial test
# null: probability of success is equal to "p"
# binom.test(number of successes, total number of trials, p)
binom.test(24,1000,p=0.015)


## calculate probability in poisson distribution
## p(x;λ) = (λ^x * e^-λ)/x!
# λ is the expected value (mean) and x is a specific observation
# ex: we expect 6 customers (this is λ) every 30 minutes (fixed amount of time); what is the probability (p) of 
# having 3 customers (this is x) in 30 minutes?
# in R terms ---> ((λ^x) * (exp(1)^-λ)) / factorial(x)



# redirect STDIN to R from command line
basename(commandArgs()[6])
dirname(commandArgs()[6])


# df and matrices for trial
dat1 <- data.frame(var1 = rnorm(10), var2 = (1:10), var3 = gl(2, 5, labels = c("red", "blue")))
dat1[,1] = as.character(dat1[,1])
class(dat1[,1])

a = matrix(c(1, 2, 3, 4, 5, 6), nrow=3, ncol=2)
b = matrix(c(7, 8, 9, 10, 11, 12), nrow=3, ncol=2)
c = matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0), nrow=5, ncol=5)

m = matrix(runif(100),10,10)
image(c, breaks=1,2,3)

# counts how many zeroes and ones(or higher) are in each row of a matrix and summarize it in a data.frame
zeroesOnesGainAll = data.frame(matrix(rep(0, dim(gainAll)[1]*2), ncol=2))
colnames(zeroesOnesGainAll) = c("zeroes", "ones")
rownames(zeroesOnesGainAll) = rownames(gainAll)
for(i in 1:dim(gainAll)[1]){
  trueFalse = table(as.numeric(gainAll[i,]) > 0)
  zeroesOnesGainAll[i,1] = as.numeric(trueFalse)[1]
  zeroesOnesGainAll[i,2] = as.numeric(trueFalse)[2]
}

# useful for trials when filtering out specific genes from gene lists
head(sort(colnames(countTable)[grep("[-]", colnames(countTable))]), n=50)
tail(sort(colnames(countTable)[grep("[.]IT", colnames(countTable))]), n=50)
length(sort(colnames(countTable)[grep("[-]", colnames(countTable))]))
table(grepl("AS7", sort(colnames(countTable)[grep("[.]AS", colnames(countTable))])))
false = which(grepl("AS", sort(colnames(countTable)[grep("[.]AS", colnames(countTable))])) == FALSE)
sort(colnames(countTable)[grep("AS", colnames(countTable))])[false]



# ?? anna tp53?
min(1-cumsum(dhyper(0:(3-1),40,19960,300)))

