import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.distributions.discrete as distr
from statsmodels.discrete.discrete_model import NegativeBinomialP, NegativeBinomial, Poisson, GeneralizedPoisson
from statsmodels.discrete.count_model import (ZeroInflatedNegativeBinomialP, ZeroInflatedPoisson,
                                              ZeroInflatedGeneralizedPoisson)
import statsmodels.discrete._diagnostics_count as dia

import statsmodels.api as sm

data = sm.datasets.scotland.load()
data.exog = sm.add_constant(data.exog)

print(data.endog.shape)
print(data.exog.shape)

gamma_model = sm.GLM(data.endog, data.exog, family=sm.families.Gamma())
gamma_results = gamma_model.fit()
print(gamma_results.summary())

expected_params = [1, 1, 0.5]
np.random.seed(987123)
nobs = 10
exog = np.ones((nobs, 2))
exog[:nobs//2, 1] = 0
# offset is used to create misspecification of the model
# for predicted probabilities conditional moment test
#offset = 0.5 * np.random.randn(nobs)
range_mix = 0.5
offset = -range_mix / 2 + range_mix * np.random.rand(nobs)
offset = 0
mu_true = np.exp(exog.dot(expected_params[:-1]) + offset)
prob_infl = 0.15
endog = distr.zinegbin.rvs(mu_true, expected_params[-1],
                           2, prob_infl, size=mu_true.shape)

#model_nb = NegativeBinomial(endog, exog, loglike_method='nb2')
model_nb = NegativeBinomialP(endog, exog, p=2)

res_nb = model_nb.fit(method='bfgs', maxiter=5000, maxfun=5000)

print(endog)
print(exog)

#print(res_nb)
#probs_nb = res_nb.predict(which='prob')
#probsm_nb = probs_nb.mean(0)

#print(probs_nb)
#print(probsm_nb)

print("endog",endog.shape)
print("exog",exog.shape)
print(res_nb.summary())
