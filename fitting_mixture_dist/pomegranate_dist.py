

import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import math
#from scipy.special import gamma, factorial
#from scipy.optimize import curve_fit
import xlsxwriter
import scipy.stats as stats
from types import SimpleNamespace
from scipy.stats import nbinom
import random
from fitter import Fitter
from pomegranate import *



mu=10
sigma=2
data1=mu+sigma*np.random.randn(10000,1)

mu=1
sigma=1
#data2=mu+sigma*np.random.randn(10000,1)
data2=np.random.exponential(5, size=(10000,1))
#data2=np.random.poisson(5, (10000,1))
#data2=np.random.beta(1,10, (10000,1))

data=np.vstack((data1,data2))
print(data2.shape,data.shape)



print('d',data.shape)

#x = numpy.arange(0, 10, .01)

#model = GeneralMixtureModel.from_samples(NormalDistribution, 2, data)
#model = GeneralMixtureModel.from_samples([ExponentialDistribution, NormalDistribution], n_components=2, X=data)
model = GeneralMixtureModel.from_samples([ExponentialDistribution, NormalDistribution], n_components=2, X=data)
#model = GeneralMixtureModel.from_samples([PoissonDistribution, NormalDistribution], n_components=2, X=data)
#model = GeneralMixtureModel.from_samples([BetaDistribution, NormalDistribution], n_components=2, X=data)


plt.figure(figsize=(8, 4))
height,xrange,patches=plt.hist(data,bins=50,density=True,facecolor='g',alpha=0.75)
plt.plot(xrange, model.distributions[0].probability(xrange), label="Distribution 1")
plt.plot(xrange, model.distributions[1].probability(xrange), label="Distribution 2")
plt.plot(xrange, model.probability(xrange), label="Mixture")
plt.legend(fontsize=14, loc='best')
plt.show()
