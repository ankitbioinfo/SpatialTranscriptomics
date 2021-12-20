import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def bimodal(x,mu1,mu2,s1,s2,p1):
    disA=np.exp(-0.5*(x-mu1)**2 / s1**2.) / (s1*np.sqrt(2.0 * np.pi))
    disB=np.exp(-0.5*(x-mu2)**2 / s2**2.) / (s2*np.sqrt(2.0 * np.pi))
    return p1*disA+(1-p1)*disB


#xdata = np.linspace(0, 4, 50)
#y = func(xdata, 2.5, 1.3, 0.5)
#rng = np.random.default_rng()
#y_noise = 0.2 * rng.normal(size=xdata.size)
#ydata = y + y_noise
#plt.plot(xdata, ydata, 'b-', label='data')

mu=6
sigma=2
data1=mu+sigma*np.random.randn(10000,1)

mu=1
sigma=1
data2=mu+sigma*np.random.randn(10000,1)
#data2=np.random.exponential(5, size=(10000,1))
#data2=np.random.poisson(5, (10000,1))
#data2=np.random.beta(1,10, (10000,1))

data=np.vstack((data1,data2))
print(data2.shape,data.shape)

hist=plt.hist(data,bins=30,density=True,facecolor='g',alpha=0.75)
x = [hist[0], 0.5*(hist[1][1:]+hist[1][:-1])]
xdata=x[1]
ydata=x[0]

#plt.plot(xdata,ydata,'k-')

popt, pcov = curve_fit(bimodal, xdata, ydata)
print('point',popt)
plt.plot(xdata, bimodal(xdata, *popt), 'k-', label='fit: mu1=%5.3f, mu2=%5.3f, s1=%5.3f,s2=%5.3f,p=%5.3f' % tuple(popt))


popt, pcov = curve_fit(bimodal, xdata, ydata, bounds=(0, [6., 1., 2,1,0.1]))
print(popt)
#plt.plot(xdata, bimodal(xdata, *popt), 'g--', label='fit: mu1=%5.3f, mu2=%5.3f, s1=%5.3f,s2=%5.3f,p=%5.3f' % tuple(popt))


plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
