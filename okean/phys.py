import numpy as np
from okean.roms import roms_tools as rt

def speed(*uv):
 if len(uv)==1: return np.abs(uv[0]) # assume it is imaginary
 elif len(uv)==2: return np.sqrt(uv[0]**2+uv[1]**2)

def ke(*uv):
 if len(uv)==1: return 0.5*(uv[0].real**2+uv[0].imag**2) # assume it is imaginary
 elif len(uv)==2: return 0.5*(uv[0]**2+uv[1]**2)


def okubo(u,v,pm,pn):
  '''
  Compute the Okubo-Weiss parameter
  '''

  Mp,Lp=pm.shape
  L=Lp-1
  M=Mp-1
  Lm=L-1
  Mm=M-1

  uom = 2.*u/(pm[:,:L]+pm[:,1:Lp])
  uon = 2.*u/(pn[:,:L]+pn[:,1:Lp])
  von = 2.*v/(pn[0:M,:]+pn[1:Mp,:])
  vom = 2.*v/(pm[0:M,:]+pm[1:Mp,:])

  mn=pm*pn
  mn_p=(mn[:M,:L]+mn[:M,1:Lp]+mn[1:Mp,1:Lp]+mn[1:Mp,:L])/4.

  # relative vorticity:
  xi=mn*rt.psi2rho(von[:,1:Lp]-von[:,:L]-uom[1:Mp,:]+uom[:M,:])

  # Sigma_T:
  ST=mn*rt.psi2rho(von[:,1:Lp]-von[:,:L]+uom[1:Mp,:]-uom[:M,:])

  # Sigma_N:
  SN=np.ma.zeros((Mp,Lp),u.dtype)
  SN[1:-1,1:-1]=mn[1:-1,1:-1]*(uon[1:-1,1:]
                              -uon[1:-1,:-1]
                              -vom[1:,1:-1]
                              +vom[:-1,1:-1])
  return SN**2+ST**2-xi**2


def vorticity(u,v,pm,pn):
  '''
  Compute the relative vorticity
  '''

  Mp,Lp=pm.shape
  L=Lp-1
  M=Mp-1
  Lm=L-1
  Mm=M-1

  uom = 2.*u/(pm[:,:L]+pm[:,1:Lp])
  uon = 2.*u/(pn[:,:L]+pn[:,1:Lp])
  von = 2.*v/(pn[0:M,:]+pn[1:Mp,:])
  vom = 2.*v/(pm[0:M,:]+pm[1:Mp,:])

  mn=pm*pn
  mn_p=(mn[:M,:L]+mn[:M,1:Lp]+mn[1:Mp,1:Lp]+mn[1:Mp,:L])/4.

  # relative vorticity:
  xi=mn*rt.psi2rho(von[:,1:Lp]-von[:,:L]-uom[1:Mp,:]+uom[:M,:])

  return xi


def hor_div(u,v,pm,pn):
  '''
  Compute the horizontal divergence
  '''
  Mp,Lp=pm.shape
  L  = Lp-1
  M  = Mp-1
  Lm = L-1
  Mm = M-1
  Dtype=u.dtype

  uom = 2.*u/(pm[:,:L]+pm[:,1:Lp])
  uon = 2.*u/(pn[:,:L]+pn[:,1:Lp])
  von = 2.*v/(pn[0:M,:]+pn[1:Mp,:])
  vom = 2.*v/(pm[0:M,:]+pm[1:Mp,:])

  mn  = pm*pn
  # Horizontal divergence:
  hdiv = mn[1:-1,1:-1]* ((vom[1:,1:L]-vom[:-1,1:L])+(uon[1:M,1:]-uon[1:M,:-1]))

  return hdiv
