
#from bayes_opt import BayesianOptimization
from skopt import Optimizer, gp_minimize 
import numpy as np
from pomegranate import *
import matplotlib.pyplot as plt



# example data of beta/gaussian distribution
data = np.hstack((np.random.beta(1, 10, size=20000),
                  np.random.randn(10000) * 0.2 + 0.6))
data = data[np.logical_and(data >= 0.0, data <= 1.0)]

print(data.shape)

def loss_bimodal(a, b, mu, sigma, w1):
    beta_loss = BetaDistribution(a, b).probability(data)
    norm_loss = NormalDistribution(mu, sigma).probability(data)
    return np.log(w1 * beta_loss + (1 - w1) * norm_loss).sum()

def pdf_bimodal(a, b, mu, sigma, w1, x=np.arange(0, 1, 0.01)):
    return w1 * BetaDistribution(a, b).probability(x) + \
        (1 - w1) * NormalDistribution(mu, sigma).probability(x)

optimizer = Optimizer(
    base_estimator=loss_bimodal,
    dimensions={'mu': (0., 2.),
             'sigma': (0., 2.),
             'a': (0, 5),
             'b': (1, 25),
             'w1': (0., 1.)},
    random_state=1
)
optimizer.maximize(
    init_points=50,
    n_iter=1000
)

a, b, mu, sigma, w1 = [v for k, v in optimizer.max['params'].items()]
x = np.arange(0, 1, 0.01)
plt.plot(x, pdf(a, b, mu, sigma, w1, x))
plt.hist(data, density=True, bins=100)
plt.show()
