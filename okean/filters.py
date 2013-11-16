'''
plfilt

apply pl33 33 hour filter to timeseries.
'''

import numpy as np

class plfilt(object):
    """docstring for plfilt"""
    
    _pl33 = np.array([-0.00027, -0.00114, -0.00211, -0.00317, -0.00427, -0.00537,
                      -0.00641, -0.00735, -0.00811, -0.00864, -0.00887, -0.00872,
                      -0.00816, -0.00714, -0.0056 , -0.00355, -0.00097,  0.00213,
                       0.00574,  0.0098 ,  0.01425,  0.01902,  0.024  ,  0.02911,
                       0.03423,  0.03923,  0.04399,  0.04842,  0.05237,  0.05576,
                       0.0585 ,  0.06051,  0.06174,  0.06215,  0.06174,  0.06051,
                       0.0585 ,  0.05576,  0.05237,  0.04842,  0.04399,  0.03923,
                       0.03423,  0.02911,  0.024  ,  0.01902,  0.01425,  0.0098 ,
                       0.00574,  0.00213, -0.00097, -0.00355, -0.0056 , -0.00714,
                      -0.00816, -0.00872, -0.00887, -0.00864, -0.00811, -0.00735,
                      -0.00641, -0.00537, -0.00427, -0.00317, -0.00211, -0.00114,
                      -0.00027], dtype='d')
    
    _dt = np.linspace(-33, 33, 67)
    
    def __init__(self, dt=1.0):
        
        if np.isscalar(dt):
            self.dt = float(dt)
        else:
            self.dt = np.diff(dt).mean()
        
        filter_time = np.arange(0.0, 33.0, 1.0/self.dt, dtype='d')
        self.Nt = len(filter_time)
        self.filter_time = np.hstack((-filter_time[-1:0:-1], filter_time))
        
        self.pl33 = np.interp(self.filter_time, self._dt, self._pl33)
        self.pl33 /= self.pl33.sum()
    
    def __call__(self, u, t=None, mode='valid'):
        uf = np.convolve(u, self.pl33, mode=mode)
        if t is None:
            return uf
        else:
            tf = t[self.Nt-1:-self.Nt+1]
            return uf, tf



def pl33(u,dt=1,t=None):
  a=plfilt(dt)
  return a(u,t)

#if __name__ == '__main__':
#    u = np.zeros((1000,),'d')
#    u[300] = 1.0
#    t = np.linspace(0, 100, 1000)
#    pl33 = plfilt(t)
#    uf, tf = pl33(u, t)
