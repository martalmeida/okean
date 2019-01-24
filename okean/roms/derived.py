import numpy as np
from . import roms_tools as rt
from .. import phys
from ..vis import Data


class Derived:
  def slice_derived(self,var,ind,time,**opts):
    plot= opts.get('plot',False)
    ax  = opts.get('ax',None)
    if not isinstance(ind,dict): # slicez,k
      coords=opts.get('coords','x,y,t').split(',')
    else: # slicell
      coords=opts.get('coords','d,z,t').split(',')

    out=Data()
    out.msg=self.check_slice('u',t=time) # any time dependent var could be
                                         # used here!
    if out.msg: return out


    if var in ('div','divergence','diverg'): var='hdiv'
    if var in ('vort','vorticity','rel vort'): var='rvort'

    if not isinstance(ind,dict): # horizontal slices, else, slicell

      if var in ('speed','ke'):
        outuv=self.sliceuv(ind,time,**opts)
        if outuv.msg: return outuv

        if   var=='speed': res=phys.speed(outuv.v)
        elif var=='ke':    res=phys.ke(outuv.v)

        out=outuv
        out.v=res
        out.info['v']['name']=var
        if var=='speed': out.info['v']['units']='metre second-1'
        elif var=='ke':  out.info['v']['units']='metre2 second-2'


        return out

      elif var in ('okubo','hdiv','rvort'):
        if ind=='bar':
          slc,uname,vname,ind = self.slicek, 'ubar','vbar', 9999
        elif ind in ('s','surf','surface'):
          slc,uname,vname,ind = self.slicek,  'u','v', -1
        elif ind>=0:
          slc,uname,vname,ind = self.slicek,  'u','v', ind
        elif ind <0:
          slc,uname,vname,ind = self.slicez, 'u','v', ind

##        if ind>=0: slc=self.slicek
##        else: slc=self.slicez
##
##
##        if ind=='bar': uname,vname,ind='ubar','vbar',0
##        else: uname,vname='u','v'

        outu=slc(uname,ind,time,**opts)
        outv=slc(vname,ind,time,**opts)

        if   outu.msg: return outu
        elif outv.msg: return outv

        # coords at rho needed for okubo','hdiv','rvort, so:
        if any([i in coords for i in 'xy']):
          x,y,h,m=self.grid.vars(ruvp='r')

        if 'x' in coords:
          if self.grid.spherical:
            out.x=x
            out.info['x']=dict(name='Longitude',units=r'$\^o$E')
          else:
            out.x=x/1000.
            out.info['x']=dict(name='Distance',units='km')

        if 'y' in coords:
          if self.grid.spherical:
            out.y=y
            out.info['y']=dict(name='Latitude',units=r'$\^o$N')
          else:
            out.y=y/1000.
            out.info['y']=dict(name='Distance',units='km')

        if 't' in coords and 't' in self.vaxes(uname): out.t=self.time[time]

        if 'z' in coords and self.hasz(uname):
          if ind>0: out.z=self.s_levels(time,k=ind,loc='rr')
          else: out.z=ind+np.zeros(out.v.shape)

          out.info['z']=dict(name='Depth',units='m')


        if var=='okubo':
          res=phys.okubo(outu.v,outv.v,self.grid.pm,self.grid.pn)

          out.v=res
          out.info['v']['name']=var
          out.info['v']['units']='second-2'

        elif var=='hdiv':
          res=phys.hor_div(outu.v,outv.v,self.grid.pm,self.grid.pn)
          out.v=res
          out.info['v']['name']=var
          out.info['v']['units']='second-1'

          if not out.x is None: out.x=out.x[1:-1,1:-1]
          if not out.y is None: out.y=out.y[1:-1,1:-1]
          if not out.z is None: out.z=out.z[1:-1,1:-1]


        elif var=='rvort':
          res=phys.vorticity(outu.v,outv.v,self.grid.pm,self.grid.pn)
          out.v=res
          out.info['v']['name']=var
          out.info['v']['units']='second-1'


        out.coordsReq=','.join(sorted(coords))

        out.set_projection(self.grid.proj_info['basemap_opts'])

        return out



    else: # slice ll:
      print('TODO !!!!')
      return out
#        X=ind['x']
#        Y=ind['y']
#        if var=='okubo':
#          x,y,z,u,v=self.sliceuv(ind,time,**opts)
#          data=okubo(rt.psi2uvr(u,'u'),rt.psi2uvr(v,'v'),self.grid.pm,self.grid.pn)
#          res=slicell(self,'r',X,Y,time,dist=True,data=data,**opts)
#
#        elif var in ['speed','ek']:
#          dist,z,u=self.slicell('u',X,Y,time,dist=True,**opts)
#          dist,z,v=self.slicell('v',X,Y,time,dist=True,**opts)
#
#          if   var=='speed': res=speed(u,v)
#          elif var=='ke':    res=ke(u,v)
#
#        return dist,z,res
