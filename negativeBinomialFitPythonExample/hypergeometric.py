from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import numpy as np


M=13
n=8
N=7
x=3

table=np.array([[x,n-x],[N-x,M-n-N+x]])
#table = np.array([[6, 2], [1, 4]])
#table = np.array([[15, 35], [135, 215]])

M = table.sum()
n = table[0].sum()
N = table[:, 0].sum()
start, end = hypergeom.support(M, n, N)
#pb1=hypergeom.pmf(np.arange(start, end+1), M, n, N)

xr=np.arange(min(n,N)+1)
pb1=hypergeom.pmf(xr,M,n,N)

print(M,n,N,xr,np.sum(pb1),sum(pb1[0:16]),sum(pb1[16:]))
print(pb1)

oddsr, p = fisher_exact(table, alternative='two-sided')
print('two sideded',p)

oddsr, p = fisher_exact(table, alternative='greater')
print('greater',p) # if x=3 then greater means 3 or greater than 3

oddsr, p = fisher_exact(table, alternative='less')
print('less',p) # if x=3 then less means 3 or less than 3 

print(table)
