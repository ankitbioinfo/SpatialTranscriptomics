

import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.special import gamma, factorial
from scipy.optimize import curve_fit
import xlsxwriter


def negativeBinomial(k,r,mu):
    # NB (r,p)
    # p is a function of r and mu
    const=gamma(r+k)/(factorial(k)*gamma(r))
    ft=(r/(r+mu))**r  # p = r/(r+mu)
    st=(mu/(r+mu))**k   # 1-p = mu/(r+mu)
    return const*ft*st

k=np.arange(51)
mu=15
r=[1,10,20,40]


fig,ax=plt.subplots(2,2,figsize=(7,6))


count=0
for i in range(2):
    for j in range(2):
        y=negativeBinomial(k,r[count],mu)
        p=r[count]/(mu+r[count])
        ax[i][j].plot(k,y,'b.-')
        ax[i][j].set_title('r = '+str(r[count])+',p = '+str('%0.2f'%p) +', mu = '+str(mu) )
        count+=1

plt.savefig('negativeBinomial.png')
