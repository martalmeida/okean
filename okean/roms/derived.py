import numpy as np
import roms_tools as rt


def speed(u,v): return np.sqrt(u**2+v**2)
def ke(u,v): return 0.5*(u**2+v**2)

def okubo(u,v,pm,pn):
  '''
  Compute the Okubo-Weiss parameter
  '''

  Mp,Lp=pm.shape
  L=Lp-1
  M=Mp-1
  Lm=L-1
  Mm=M-1
  Dtype=u.dtype

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
  SN=np.zeros((Mp,Lp),Dtype)
  SN[1:-1,1:-1]=mn[1:-1,1:-1]*(uon[1:-1,1:]
                              -uon[1:-1,:-1]
                              -vom[1:,1:-1]
                              +vom[:-1,1:-1])
  return SN**2+ST**2-xi**2


class Derived:
  def slice_derived(self,var,ind,time,**opts):
    if var in ('speed','ke','okubo'):
      if not isinstance(ind,dict):
        x,y,z,u,v=self.sliceuv(ind,time,**opts)
        if   var=='speed': res=speed(u,v)
        elif var=='ke':    res=ke(u,v)
        elif var=='okubo':
          res=okubo(rt.psi2uvr(u,'u'),rt.psi2uvr(v,'v'),self.grid.pm,self.grid.pn)
          x,y=self.grid.lon,self.grid.lat
        return x,y,z,res
      else: # slice ll:
        X=ind['x']
        Y=ind['y']
        dist,z,u=self.slicell('u',X,Y,time,dist=True,**opts)
        dist,z,v=self.slicell('v',X,Y,time,dist=True,**opts)

        if   var=='speed': res=speed(u,v)
        elif var=='ke':    res=ke(u,v)
        elif var=='okubo':
          res=okubo(rt.psi2uvr(u,'u'),rt.psi2uvr(v,'v'),self.grid.pm,self.grid.pn)
          x,y=self.grid.lon,self.grid.lat

        return dist,z,res
