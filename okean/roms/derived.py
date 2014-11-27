import numpy as np
import roms_tools as rt
from okean import phys



class Derived:
  def slice_derived(self,var,ind,time,**opts):
      if not isinstance(ind,dict): # horizontal slices

        if var in ('speed','ke'):
          x,y,z,u,v=self.sliceuv(ind,time,**opts)
        elif var in ('okubo','hdiv','vort'):
          if ind>=0: slc=self.slicek
          else: slc=self.slicez
          xu,yu,zu,u=slc(u,ind,time,**opts)
          xv,yv,zv,v=slc(v,ind,time,**opts)

        if   var=='speed': res=phys.speed(u,v)
        elif var=='ke':    res=phys.ke(u,v)
        elif var=='okubo':
          res=phys.okubo(u,v,self.grid.pm,self.grid.pn)
          x,y=self.grid.lon,self.grid.lat
        elif var=='hdiv':
          res=phys.hor_div(u,v,self.grid.pm,self.grid.pn)
          x,y=self.grid.lon,self.grid.lat
        elif var=='rvort':
          res=phys.vorticity(u,v,self.grid.pm,self.grid.pn)
          x,y=self.grid.lon,self.grid.lat

        return x,y,z,res

      else: # slice ll:
        X=ind['x']
        Y=ind['y']
        if var=='okubo':
          x,y,z,u,v=self.sliceuv(ind,time,**opts)
          data=okubo(rt.psi2uvr(u,'u'),rt.psi2uvr(v,'v'),self.grid.pm,self.grid.pn)
          res=slicell(self,'r',X,Y,time,dist=True,data=data,**opts):

        elif var in ['speed','ek']:
          dist,z,u=self.slicell('u',X,Y,time,dist=True,**opts)
          dist,z,v=self.slicell('v',X,Y,time,dist=True,**opts)

          if   var=='speed': res=speed(u,v)
          elif var=='ke':    res=ke(u,v)

        return dist,z,res
