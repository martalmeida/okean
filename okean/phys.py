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

  dx=1/pm
  dy=1/pn

  dx=(dx[:,1:]+dx[:,:-1])/2 # at u
  dx=(dx[1:]+dx[:-1])/2 # at psi

  dy=(dy[:,1:]+dy[:,:-1])/2 # at u
  dy=(dy[1:]+dy[:-1])/2 # at psi

  dvdx=(v[:,1:]-v[:,:-1])/dx
  dudy=(u[1:]-u[:-1])/dy

  return dvdx-dudy

  ## or, using roms_tools:

  #Mp,Lp=pm.shape
  #L=Lp-1
  #M=Mp-1
  #Lm=L-1
  #Mm=M-1
  #
  #uom = 2.*u/(pm[:,:L]+pm[:,1:Lp])
  #uon = 2.*u/(pn[:,:L]+pn[:,1:Lp])
  #von = 2.*v/(pn[0:M,:]+pn[1:Mp,:])
  #vom = 2.*v/(pm[0:M,:]+pm[1:Mp,:])
  #
  #mn=pm*pn
  #mn_p=(mn[:M,:L]+mn[:M,1:Lp]+mn[1:Mp,1:Lp]+mn[1:Mp,:L])/4.
  #
  #xi=mn_p*(von[:,1:Lp]-von[:,:L]-uom[1:Mp,:]+uom[:M,:])
  ## at rho:
  #xi=rt.psi2rho(xi)
  #
  #return xi


def hor_div(u,v,pm,pn):
  '''
  Compute the horizontal divergence
  '''
  dx=1/pm
  dy=1/pn

  dudx=(u[:,1:]-u[:,:-1])/dx[:,1:-1] # at rho, excluding 1st and last cols
  dvdy=(v[1:]-v[:-1])/dy[1:-1] # at rho, excluding 1st and last rows

  dx=(dx[:,1:]+dx[:,:-1])/2 # at u
  dx=(dx[1:]+dx[:-1])/2 # at psi

  dy=(dy[:,1:]+dy[:,:-1])/2 # at u
  dy=(dy[1:]+dy[:-1])/2 # at psi

  dvdx=(v[:,1:]-v[:,:-1])/dx
  dudy=(u[1:]-u[:-1])/dy

  hdiv=dudx[1:-1]+dvdy[:,1:-1]
  return hdiv

  ## or, using roms_tools:

  #Mp,Lp=pm.shape
  #L  = Lp-1
  #M  = Mp-1
  #Lm = L-1
  #Mm = M-1
  #
  #uom = 2.*u/(pm[:,:L]+pm[:,1:Lp])
  #uon = 2.*u/(pn[:,:L]+pn[:,1:Lp])
  #von = 2.*v/(pn[0:M,:]+pn[1:Mp,:])
  #vom = 2.*v/(pm[0:M,:]+pm[1:Mp,:])
  #
  #mn  = pm*pn
  ## Horizontal divergence:
  #hdiv = mn[1:-1,1:-1]* ((vom[1:,1:L]-vom[:-1,1:L])+(uon[1:M,1:]-uon[1:M,:-1]))
  #return hdiv
