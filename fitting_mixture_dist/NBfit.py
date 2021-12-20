import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma, factorial
from statsmodels.discrete.discrete_model import NegativeBinomial as NB

def NBpmf(x,r,p):
    #x must be integer
    disA=[]
    for k in range(len(x)):
        const=np.math.gamma(r+k)/(np.math.factorial(k)*np.math.gamma(r))
        value=const* ((1-p)**k)*(p**r)
        disA.append(value)
    disA=np.array(disA)
    return disA
  
  
data2=np.random.negative_binomial(2, 0.15, size=(10000,1)) # parameter is n and p
 
hist=plt.hist(data2,bins=np.arange(50),density=True,facecolor='g',alpha=0.75)
#x = [hist[0], 0.5*(hist[1][1:]+hist[1][:-1])]; xdata=x[1];ydata=x[0];

xdata=hist[1][0:-1];ydata=hist[0];
print(xdata)

popt, pcov = curve_fit(NBpmf, xdata, ydata)
print('point1',popt)
plt.plot(xdata, NBpmf(xdata, *popt), 'k-', label='fit: r=%5.3f, p=%5.3f' % tuple(popt))
#parameters, covariance = curve_fit(negativeBinomial, xdata, ydata)
#print(parameters)
#yfit=negativeBinomial(xdata,1)

popt, pcov = curve_fit(NBpmf, xdata, ydata, bounds=(0, [4,0.5]))
print('point2',popt)
plt.plot(xdata, NBpmf(xdata, *popt), 'g--', label='fit: r=%5.3f, p=%5.3f' % tuple(popt))

#plt.plot(xdata, yfit, '-', label='fit')

plt.legend(fontsize=10, loc='best')

plt.show()

  
  
