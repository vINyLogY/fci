#!/usr/bin/env python
from scipy.sparse import random
from scipy import stats
class CustomRandomState(object):
    def randint(self, k):
        i = np.random.randint(k)
        return i - i % 2


rs = CustomRandomState()
rvs = stats.poisson(25, loc=10).rvs
S = random(3, 4, density=0.25, random_state=rs, data_rvs=rvs)

print S

