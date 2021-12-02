import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.distributions.discrete as distr
from statsmodels.discrete.discrete_model import NegativeBinomialP, NegativeBinomial, Poisson, GeneralizedPoisson
from statsmodels.discrete.count_model import (ZeroInflatedNegativeBinomialP, ZeroInflatedPoisson,
                                              ZeroInflatedGeneralizedPoisson)
import statsmodels.discrete._diagnostics_count as dia

import statsmodels.api as sm


df=pd.read_csv('sampleNBdata.dat')
#print(df.keys())
data=pd.concat((df,pd.get_dummies(df['prog'],drop_first=False)),axis=1)
endog=data['daysabs']
data['intercept'] = 1
exog=data.drop(['prog','daysabs','id','gender','Unnamed: 0','General'],axis=1)
#exog=pd.concat(df['math'],df['prog'])
#exog=np.concatenate(df['math'],df['prog'])


#print(data.keys())
#print(exog)

print("endog",endog.shape)
print("exog",exog.shape)


#model_nb = NegativeBinomial(endog, exog, loglike_method='nb2')
model_nb = NegativeBinomialP(endog, exog, p=2)
res_nb = model_nb.fit(method='bfgs', maxiter=5000, maxfun=5000)
print(res_nb.summary())
