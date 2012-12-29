from numpy import log, sqrt
from fcython import fcython_quad

f = lambda x: log(x) / sqrt(x)

print fcython_quad(f, 0, 1) 

class fc_integrand(object):
    def __init__(self, alpha):
        self.alpha = alpha
        
    def __call__(self, x):
        return log(self.alpha * x) / sqrt(x)
        
fc_igrd = fc_integrand(1)
print fcython_quad(fc_igrd, 0, 1) 

print fcython_quad(f, 0, 1)

from scipy.integrate import quad

print quad(lambda x: log(x) / sqrt(x), 0, 1)[0] # -4.0

class integrand(object):
    def __init__(self, alpha):
        self.alpha = alpha
        
    def __call__(self, x):
        return log(self.alpha * x) / sqrt(x)
        
igrd = integrand(1)
print quad(igrd, 0, 1)[0] # -4.0