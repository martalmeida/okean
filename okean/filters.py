'''
plfilt

apply pl33 33 hour filter to timeseries.
'''

import numpy as np
from scipy.signal import detrend as dtrend

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
        
        #filter_time = np.arange(0.0, 33.0, 1.0/self.dt, dtype='d')
        filter_time = np.arange(0.0, 33.0, self.dt, dtype='d')
        self.Nt = len(filter_time)
        self.filter_time = np.hstack((-filter_time[-1:0:-1], filter_time))
        
        self.pl33 = np.interp(self.filter_time, self._dt, self._pl33)
        self.pl33 /= self.pl33.sum()
    
    def __call__(self, u, t=None, mode='valid', detrend = False):
        if detrend:
            u_dt = dtrend(u, type = 'linear')
            du = u - u_dt
            u = u_dt
        else:
            du = 0
        uf = np.convolve(u, self.pl33, mode=mode)+du
        if t is None:
            return uf
        else:
            tf = {'valid' : t[self.Nt-1:-self.Nt+1], 'same' : t}
            return uf, tf[mode]



def pl33(u,dt=1,t=None, extend_mode = 'same', detrend = False):
    """Applies a tide-killer type convolution filter. 
    
    The pl33 filter is described on p. 21, Rosenfeld (1983), WHOI
    Technical Report 85-35. Filter half amplitude period = 33 h., half
    power period = 38 h.
    
     Parameters
    ----------
    u : (M,) array_like
        Time series signal to filter.
    dt : float (optional)
        Sampling rate in multiples of the unit used in `t`. Default
        is 1. This usually is represented in multiples of hrs. If the
        series is sampled every hour and is set to 1 h, the filter will apply
        a 33 h half amplitude period (`T`). `dt` can be set different from the
        actual sampling rate to change the cutoff period, e.g. if sampling
        rate is equal to 1h and dt is set to 0.5, `T` will be equal to 66 h.
    t : (M,) array_like (optional)
        Time stamp series with the same length of `u`.
    extend_mode : {`valid`, `same`}, (optional)
        Keyword for the type of padding made by the convolution function.
        For a window W = 1/dt*33, N = M-2*W.
        `same`:
            Mode `same` returns output of length max(M, N). Boundary effects
            are still visible.
        `valid`:
            Mode `valid` returns output of length max(M, N) - min(M, N) + 1.
            The convolution product is only given for points where the signals
            overlap completely. Values outside the signal boundary have no
            effect.
    detrend : bool
        Apply linear detrend prior to applying convolution.

    Returns
    -------
    uf : numpy.ndarray
        Filtered signal.
    tf : numpy.ndarray
        Time stamps cut to the same length as the signal.
    """
    
    a=plfilt(dt)
    return a(u,t, mode = extend_mode, detrend = detrend)

#if __name__ == '__main__':
#    u = np.zeros((1000,),'d')
#    u[300] = 1.0
#    t = np.linspace(0, 100, 1000)
#    pl33 = plfilt(t)
#    uf, tf = pl33(u, t)
