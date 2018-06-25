## TO CALCULATE HYPERGEOMETRIC P-VALUE (python2)
## (check 'http://blog.nextgenetics.net/?e=94')

# RUN AS: "python script.py [total number of genes in the list] [total number of genes in the list with condition A] [total number of genes in the list with condition B] [number of genes with both condition A and B]"

# which, for us, means: "python script.py [total number of genes] [total number of genes belonging to the pathway] [total number of altered genes] [number of altered genes belonging to the pathway]"

# or, following the standard nomenclature: "python script.py [M] [n] [N] [x]"

# the function "hypergeom.cdf()" computes the cumulative probability of having a number of altered genes in pathway <= of the number we actually found (i.e. if we found 4 genes, it computes P(x=0) + P(x=1) + P(x=2) + P(x=3) + P(x=4))

# the function "hypergeom.sf()" computes the probability of having x > than the number we found (i.e. if  we found 4 genes, it computes P(x=5) + P(x=6) + ... + P(x=n))

# Therefore ---> 1 = hypergeom.cdf(x, M, n, N) + hypergeom.sf(x, M, n, N)


# IF SO, I D0N'T UNDERSTAND WHY HERE IS "hypergeom.cdf(x+1, M, n, N)", instead of "hypergeom.cdf(x, M, n, N)"

import sys
import scipy.stats as stats

print
print 'total number in population: ' + sys.argv[1]
print 'total number with condition in population: ' + sys.argv[2]
print 'number in subset: ' + sys.argv[3]
print 'number with condition in subset: ' + sys.argv[4]
print
print 'p-value <= ' + sys.argv[4] + ': ' + str(stats.hypergeom.cdf(int(sys.argv[4]) + 1,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])))
print 
print 'p-value >= ' + sys.argv[4] + ': ' + str(stats.hypergeom.sf(int(sys.argv[4]) - 1,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])))
print