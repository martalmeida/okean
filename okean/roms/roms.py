import numpy as np
#import pylab as pl
import os

from okean import calc, netcdf, ticks
import roms_tools as rt
from derived import Derived
from okean.vis import Data

try:
  from mpl_toolkits.basemap import Basemap
except:
  Basemap=False

__all__='His','Grid'

#class Aux:
#  def __init__(self):
#    self.x=None
#    self.y=None
#    self.d=None
#    self.z=None
#    self.t=None
#    self.tnum=None # used with time_series (2d times)
#    self.msg=''
#
#  def info(self):
#    for i in ['x','y','d','z','t','tnum']:
#      a=getattr(self,i)
#      try: print '%5s shape=%12s size=%d'%(i,str(a.shape),a.size)
#      except: print '%5s %s'%(i,str(a))
#
#    if not self.msg: msg='EMPTY'
#    else: msg=self.msg
#    print ' == msg == %s'%msg
    

class Common:
  def load_vars(self,vars):
    #self.quiet=0
    var=netcdf.var(self.nc)
    varnames=var.keys()
    for i in vars:
      if not self.quiet: print ':: loading var ',  i

      if isinstance(i,dict):
        iname, iopts = i.items()[0]
      elif isinstance(i,basestring):
        iname=i
        if iname.startswith('@'): iopts=iname[1:]
        else: iopts=iname

      if not isinstance(iopts,tuple):
        iopts=(iopts,)

      found=0
      if not hasattr(self,'var_as'): self.var_as={}
      for j in iopts:
        if j in varnames:
          if iname.startswith('@'):
            iname=iname[1:]
            setattr(Common, iname, property(fget=lambda self: self.use(j)))
            # not working as expected! self.use is executed in the end, for last variable only!
            # check later...
          else:
            setattr(self,iname,var[j][:])

          found=1
          if not self.quiet: print '  - found %s as %s'  % (iname,j)

          # store original name and other info, like units
          try: units=var[j].atts['units'].value
          except: units=''
          self.var_as[iname]={'name':j,'units':units}
          break
      if not found and not self.quiet: print ':: %s not found in %s' % (str(iopts),self.name)


  def load_dims(self):
    dms=netcdf.fdim(self.nc)
    for k in dms.keys():
      if k in ('ocean_time','time','scrum_time'):
        setattr(self,'TIME',dms[k])
      else:
        setattr(self,k.upper(),dms[k])

  def load_atts(self):
    atts=netcdf.fatt(self.nc)
    self.atts={}
    for k in atts.keys(): self.atts[k]=atts[k].value

  def load_grid(self,grd=False):
    ats=netcdf.fatt(self.nc)
    if not grd:
      if 'grd_file' in ats.keys(): grd=ats['grd_file'].value
      if grd and not os.path.isfile(grd): grd=False

    if grd: self.grid=Grid(grd,self.quiet)
    else:
      try:
        self.grid=Grid(self.name)
      except:
        self.grid=False

  def show(self): netcdf.show(self.nc)
  def showvar(self,vname): netcdf.showvar(self.nc,vname)
  def use(self,varname,**kargs): return netcdf.use(self.nc,varname,**kargs)

  def hasz(self,v):
    '''
    True if variable has depth dimension (s_rho or s_w)
    '''

    dnames=netcdf.vdim(self.nc,v)
    return 's_rho' in dnames or 's_w' in dnames

  def hast(self,v):
    '''
    True if variable has *time* dimension
    '''

    dnames=netcdf.vdim(self.nc,v)
    for d in dnames:
      if d.find('time')>=0: return True

    return False

  def var_at(self,v):##,_3d=True):
    '''
    returns location of the variable in the 2d/3d grid
    Possible values are u,v,r(rho), p(psi)
    If _3d, also checks if at w points
    '''
    dims=netcdf.vdim(self.nc,v).keys()

    hLoc='r' # r,u,v
    vLoc='r' # r,w

    if 's_w' in dims: vLoc='w'
    if 'xi_u'     in dims or   'eta_u'   in dims: hLoc='u' # or instead of and cos of agrif ...
    elif 'xi_v'   in dims or   'eta_v'   in dims: hLoc='v'
    elif 'xi_psi' in dims and  'eta_psi' in dims: hLoc='p'
    return hLoc,vLoc

    #if _3d and 's_w' in dims: return 'w'
    #elif 'xi_u'   in dims or   'eta_u'   in dims: return 'u'
    #elif 'xi_v'   in dims or   'eta_v'   in dims: return 'v'
    #elif 'xi_psi' in dims and  'eta_psi' in dims: return 'p'
    #elif 'xi_rho' in dims and  'eta_rho' in dims: return 'r'
    #
    #return False


class Grid(Common):
  def __init__(self,grd,quiet=True):
    if grd.startswith('http'):
      self.name=grd
      self.isremote=True
    else:
      self.name=os.path.realpath(grd)
      self.isremote=False

    self.type='grid'
    self.quiet=quiet

    self.nc=netcdf.ncopen(self.name)
    self.load()
    self.load_dims()
    self.load_uvp()
    self.load_atts()


  def load(self):
    self.isremote=1
    if self.isremote: # load as few vars as possible! @ means loaded only when needed!
      vars={'lon':('lon_rho','x_rho')},{'lat':('lat_rho','y_rho')},\
           {'mask':'mask_rho'},'h','angle','pm','pn'
#           {'mask':'mask_rho'},'h','@angle','@pm','@pn'

    else:
      vars={'lon':('lon_rho','x_rho')}, {'lonu':('lon_u','x_u')},\
           {'lonv':('lon_v','x_v')},    {'lonp':('lon_psi','x_psi')},\
           {'lat':('lat_rho','y_rho')}, {'latu':('lat_u','y_u')},\
           {'latv':('lat_v','y_v')},    {'latp':('lat_psi','y_psi')},\
           {'mask':'mask_rho'},{'masku':'mask_u'},{'maskv':'mask_v'},{'maskp':'mask_psi'},\
           'h','angle', 'pm','pn'
#           'h','@angle', '@pm','@pn'

    self.load_vars(vars)

    # some vars may not be present... set defaults:
    if not hasattr(self,'mask'): self.mask=np.ones(self.h.shape)
    if not hasattr(self,'angle'): self.angle=np.zeros(self.h.shape)

    # add spherical attr:
    sph=self.use('spherical')
    # may be 1, T, an array with 'T','','',... etc
    try:
      sph=''.join(sph).strip()
      if len(sph): sph=sph[0]
    except: pass
    trueVals=[1,'T']
    # let us assume that empty sph means True!! (may not be a very
    # smart idea!
    trueVals+=['']
    if sph in trueVals: self.is_spherical=True
    else: self.is_spherical=False


  def load_grid(): self.load()

  def load_uvp(self):
    if not hasattr(self,'lonu'): self.lonu=rt.rho2uvp(self.lon,'u')
    if not hasattr(self,'lonv'): self.lonv=rt.rho2uvp(self.lon,'v')
    if not hasattr(self,'lonp'): self.lonp=rt.rho2uvp(self.lon,'p')

    if not hasattr(self,'latu'): self.latu=rt.rho2uvp(self.lat,'u')
    if not hasattr(self,'latv'): self.latv=rt.rho2uvp(self.lat,'v')
    if not hasattr(self,'latp'): self.latp=rt.rho2uvp(self.lat,'p')

    mu,mv,mp=rt.uvp_mask(self.mask)
    if not hasattr(self,'masku'): self.masku=mu
    if not hasattr(self,'maskv'): self.maskv=mv
    if not hasattr(self,'maskp'): self.maskp=mp

    if not hasattr(self,'XI_U'): self.XI_U=self.XI_RHO-1
    if not hasattr(self,'XI_V'): self.XI_V=self.XI_RHO
    if not hasattr(self,'ETA_U'): self.ETA_U=self.ETA_RHO
    if not hasattr(self,'ETA_V'): self.ETA_V=self.ETA_RHO-1


  def border_multi_corner(self):
    mc=self.use('maskc')
    l,n=mc.shape
    print 
    m=np.zeros((l+2,n+2))
    x=np.zeros((l+2,n+2))
    y=np.zeros((l+2,n+2))

    m[1:-1,1:-1]=mc
    x[1:-1,1:-1]=self.lon
    y[1:-1,1:-1]=self.lat

    x[:,0]  = x[:,1]
    x[:,-1] = x[:,-2]
    x[0]    = x[1]
    x[-1]   = x[-2]

    y[:,0]  = y[:,1]
    y[:,-1] = y[:,-2]
    y[0]    = y[1]
    y[-1]   = y[-2]

    import pylab as pl
    f=pl.figure()
    f.set_visible(0)
    ax=pl.axes()
    cs=ax.contour(x,y,m,[.99999])
    p = cs.collections[0].get_paths()[0]
    v = p.vertices
    pl.close(f)

    return v[:,0],v[:,1]


  def border(self,type='r',margin=0,di=1,dj=1):
    if   type == 'r': lon,lat=self.lon, self.lat
    elif type == 'u': lon,lat=self.lonu,self.latu
    elif type == 'v': lon,lat=self.lonv,self.latv
    elif type == 'p': lon,lat=self.lonp,self.latp

    if margin>0:
      m,n=lon.shape
      lon=lon[margin:m-margin,margin:n-margin]
      lat=lat[margin:m-margin,margin:n-margin]
    elif margin<0: # expand grid in multiples of last cell
      n=-margin
      x=lon.copy()
      y=lat.copy()

      x[0,:]=(n+1)*x[0,:]-n*x[1,:]
      y[0,:]=(n+1)*y[0,:]-n*y[1,:]

      x[-1,:]=(n+1)*x[-1,:]-n*x[-2,:]
      y[-1,:]=(n+1)*y[-1,:]-n*y[-2,:]

      x[:,0]=(n+1)*x[:,0]-n*x[:,1]
      y[:,0]=(n+1)*y[:,0]-n*y[:,1]

      x[:,-1]=(n+1)*x[:,-1]-n*x[:,-2]
      y[:,-1]=(n+1)*y[:,-1]-n*y[:,-2]

      lon=x
      lat=y

    xb=calc.var_border(lon,di,dj)
    yb=calc.var_border(lat,di,dj)
    return xb,yb

  def vars(self,ruvp='r',i=False,j=False):
    ruvp=ruvp[0]
    isU=ruvp=='u'
    isV=ruvp=='v'
    isP=ruvp=='p'


    if isU:
      x = self.lonu
      y = self.latu
      h = rt.rho2uvp(self.h,'u')
      m = self.masku
    elif isV:
      x = self.lonv
      y = self.latv
      h = rt.rho2uvp(self.h,'v')
      m = self.maskv
    elif isP:
      x = self.lonp
      y = self.latp
      h = rt.rho2uvp(self.h,'p')
      m = self.maskp
    else:
      x = self.lon
      y = self.lat
      h = self.h
      m = self.mask

    if j is False: j=slice(None)
    if i is False: i=slice(None)

    return x[j,i],y[j,i],h[j,i],m[j,i]

#####  def s_levels(self,sparams,zeta=0,h=False,ruvpw='r',i=False,j=False,k=False):
  def s_levels(self,sparams,zeta=0,h=False,loc='rr',i=False,j=False,k=False):
    hLoc,vLoc=loc
####    ruvpw=ruvpw[0]
####    isW=ruvpw=='w'

    if h is False:
      h=self.h

    try:
      zeta.shape==h.shape
    except:
      zeta=np.tile(zeta,h.shape).astype(h.dtype)

    if h.ndim==2:
      h=rt.rho2uvp(h,hLoc)
      zeta=rt.rho2uvp(zeta,hLoc)

    z=rt.s_levels(h,zeta,sparams,rw=vLoc)
###    zr,zw=rt.s_levels(h,zeta,sparams)
###    if isW:z=zw
###    if vLoc=='w': z=zw
###    else: z=zr

    if k is False: k=slice(None)
    if j is False: j=slice(None)
    if i is False: i=slice(None)

    return np.squeeze(z[k,j,i])


#  def inside(self,x,y,data,ij=2):
#    d1=(x-self.lon.min())**2 + (y-self.lat.min())**2
#    d2=(x-self.lon.max())**2 + (y-self.lat.min())**2
#    d3=(x-self.lon.max())**2 + (y-self.lat.max())**2
#    d4=(x-self.lon.min())**2 + (y-self.lat.max())**2
#
#    i1,j1=np.where(d1==d1.min())
#    i2,j2=np.where(d2==d2.min())
#    i3,j3=np.where(d3==d3.min())
#    i4,j4=np.where(d4==d4.min())
#
#    i1,j1=i1[0],j1[0]
#    i2,j2=i2[0],j2[0]
#    i3,j3=i3[0],j3[0]
#    i4,j4=i4[0],j4[0]
#
#    i1=min((i1,i2,i3,i4))
#    i2=max((i1,i2,i3,i4))
#    j1=min((j1,j2,j3,j4))
#    j2=max((j1,j2,j3,j4))
#
#    i1=max((0,i1-ij))
#    j1=max((0,j1-ij))
#    i2=min((x.shape[0],i2+ij))
#    j2=min((x.shape[1],j2+ij))
#
#    return x[i1:i2,j1:j2],y[i1:i2,j1:j2],data[i1:i2,j1:j2]


  def ingrid(self,x,y):#,retInds=False):
    '''in polygon of region border'''
    xb,yb=self.border()
    return calc.inpolygon(x,y,xb,yb)
#    i=calc.inpolygon(x,y,xb,yb)
#    if retInds:
#      return x[i!=0], y[i!=0], np.where(i)[0]
#    else: return x[i!=0], y[i!=0]


  def resolution(self):
    dx=1./self.pm
    dy=1./self.pn
    dx=np.ma.masked_where(self.mask==0,dx)
    dy=np.ma.masked_where(self.mask==0,dy)
    return self.lon,self.lat,dx,dy


  def plot3d(self,type='surface',**kargs):
    nmax=kargs.get('nmax',200) # max wireframe lines

    if type=='surface':

      # mayavi and pylab version<3 needs the environment variable
      # QT_API=pyqt, otherwise there will be the error "ValueError: API
      # 'QString' has already been set to version 1
      # This is because ipython uses pyqt4 using the version 1 of the api
      # see: http://mail.scipy.org/pipermail/ipython-user/2012-February/009478.html

      import sys
      pver=eval(sys.version[:3]) 
      qtapi=os.environ.get('QT_API','')

      if pver<3 and qtapi!='pyqt':
        print 'warning: environmen variable QT_API may need to be set as pyqt '\
              'in python <3 so that mayavi can be used after pylab' \
 
        a=raw_input('Wanna continue ([n],y) ?')
        if a!='y': return 
      
      from mayavi import mlab

      z=-self.h.copy()
      z[self.mask==0]=np.nan

      rx=self.lon.max()-self.lon.min()
      ry=self.lat.max()-self.lat.min()
      r=2*(self.h.max()-self.h.min())/(0.5*(rx+ry))

      res=mlab.mesh(self.lon,self.lat,z/r)
      return res

    else: # wireframe
      from mpl_toolkits.mplot3d import axes3d
      import pylab as pl
      fig=pl.figure()

      ax=fig.add_subplot(111, projection='3d')

      ny,nx=self.h.shape
      dx=max(int(np.round(nx/float(nmax))),1)
      dy=max(int(np.round(ny/float(nmax))),1)

      x=self.lon[::dy,::dx]
      y=self.lat[::dy,::dx]
      h=-self.h[::dy,::dx]

      h[self.mask[::dy,::dx]==0]=np.nan

      ax.plot_wireframe(x,y,h,lw=0.5,alpha=.5)
      ax.contour3D(self.lon,self.lat,self.mask,[.5])
      return ax


  def plot(self,**kargs):
    from okean import vis
    from matplotlib.cm import gist_earth_r
    h=np.ma.masked_where(self.mask==0,self.h)
    a=vis.Data(x=self.lon,y=self.lat,v=h)
    ht=ticks.loose_label_n(self.h.min(),self.h.max(),7)
    a.set_param(field__plot='contourf',field__cvals=ht,field__cmap=gist_earth_r)
    a.set_param(**kargs)

    # also show domain boundary:
    xb,yb=self.border()
    b=vis.Data(x=xb,v=yb)
    b.set_param(d1_line__options=dict(lw=0.5,color='k',ls='-'))
    a.extra=[b]

    a.plot(labels=0)

    return a


  def plot_old(self,**kargs):
    import pylab as pl
    bathy      = kargs.get('bathy','contourf')
    hvals      = kargs.get('hvals','auto')
    resolution = kargs.get('res','i')
    rivers     = kargs.get('rivers',False)
    parallels  = kargs.get('parallels','auto')
    meridians  = kargs.get('meridians','auto')
    scale      = kargs.get('scale',False)
    states     = kargs.get('states',False)
    xlim       = kargs.get('xlim',False)
    ylim       = kargs.get('ylim',False)
    title      = kargs.get('title','auto')
    cmap       = kargs.get('cmap',pl.cm.gist_earth_r)
    proj       = kargs.get('proj','merc')
    ax         = kargs.get('ax',False)

    if not ax:
      fig=pl.figure()
      ax=pl.axes()
    else: fig=ax.figure

    h=np.ma.masked_where(self.mask==0,self.h)
    xb,yb=self.border()

    xt=ticks.loose_label(self.lon.min(),self.lon.max())
    yt=ticks.loose_label(self.lat.min(),self.lat.max())
    ht=ticks.loose_label_n(self.h.min(),self.h.max(),7)

    if xlim: Lonlims=xlim
    else:  Lonlims=xt[0],xt[-1]

    if ylim: Latlims=ylim
    else: Latlims=yt[0],yt[-1]

    if hvals=='auto': hvals=ht

    if parallels=='auto': parallels=yt
    if meridians=='auto': meridians=xt

    if Basemap and proj:

      key=Lonlims,Latlims,proj,resolution
      if not hasattr(self,'_proj'):
        self._proj={}
      
      if self._proj.has_key(key):
        m=self._proj[key]
      else:
        m = Basemap(projection=proj,lat_ts=0.0,
                  #lon_0=self.lon[0].mean(),
                  resolution=resolution,
                  urcrnrlon=Lonlims[1], urcrnrlat=Latlims[1],
                  llcrnrlon=Lonlims[0], llcrnrlat=Latlims[0])

        self._proj[key]=m

      m.drawcoastlines(ax=ax)
      m.fillcontinents(ax=ax,color=(0.7604,    1.0000,    0.7459))
      if rivers: m.drawrivers(color='b')

      m.drawcountries(ax=ax)
      # m.drawlsmask()
      # m.drawmapboundary()
      if scale:
        dx=Lonlims[1]-Lonlims[0]
        dy=Latlims[1]-Latlims[0]
        lon=Lonlims[1]-dx/10
        lat=Latlims[0]+dy/10
        lon0=lon
        lat0=lat
        length=100

        m.drawmapscale(lon,lat,lon0,lat0,length,ax=ax)

      if states: m.drawstates(ax=ax)

      m.drawparallels(parallels, labels=[1,0,0,0],ax=ax)
      m.drawmeridians(meridians, labels=[0,0,0,1],ax=ax)

      x, y = m(self.lon, self.lat)
      pch=False
      if bathy in ['pcolor','pcolormesh']:
        pch = ax.pcolormesh(x,y,h,cmap=cmap)
      elif bathy=='contour':
        pch = ax.contour(x,y,h,hvals,cmap=cmap)
      elif bathy=='contourf':
        pch = ax.contourf(x,y,h,hvals,cmap=cmap)

      if pch:
        if m.xmax-m.xmin>m.ymax-m.ymin: orientation='horizontal'
        else: orientation='vertical'

        cbh = pl.colorbar(pch,orientation=orientation,ax=ax)
        cbh.set_label('Depth',fontsize=10)

        if title=='auto': title=self.name
        if len(title)>80: font=8
        elif len(title)>60: font=9
        else: font=10
        ax.set_title(title,fontsize=font)

      xb, yb = m(xb,yb)
      ax.plot(xb,yb)
      ax.axis([m.xmin, m.xmax, m.ymin, m.ymax])
      return m

    else:
      ax.pcolormesh(self.lon,self.lat,h,cmap=cmap)
      ax.contour(self.lon,self.lat,self.h,hvals,colors='w')
      ax.plot(xb,yb)


class His(Common,Derived):
  def __init__(self,his,grd=False,quiet=True):
    #try: self.name=os.path.realpath(his)
    #except: self.name=his

    self.name=his
    if isinstance(his,basestring) and not his.startswith('http:'):
      self.name=os.path.realpath(his) # local file

    self.type='his'
    self.quiet=quiet

    self.nc=netcdf.ncopen(self.name)
    self.load_grid(grd)
    self.load()
    self.load_dims()
    self.load_atts()

  def load_grid(self,grd=False):
    Common.load_grid(self,grd)

  def load(self):
    # time:
    vars={'time':('time','ocean_time','scrum_time','clim_time')},
    self.load_vars(vars)

    # time:
    try:
      # units must have the format <units> since <date>
      self.time=netcdf.num2date(self.time,self.var_as['time']['units'])
    except: pass

    # s params:
    self.s_params=rt.s_params(self.nc)


  def path_s_levels(self,time,x,y,rw='r',**kargs):
    inds=kargs.get('inds',None)

    xr,yr,h,m=self.grid.vars()
    zeta=self.use('zeta',SEARCHtime=time)

    if inds:
      i0,i1=inds['xi']
      j0,j1=inds['eta']

      xr=xr[j0:j1,i0:i1]
      yr=yr[j0:j1,i0:i1]
      h=h[j0:j1,i0:i1]
      m=m[j0:j1,i0:i1]
      zeta= zeta[j0:j1,i0:i1]

    if np.ma.isMA(zeta): h=np.ma.masked_where(zeta.mask,h)
    else:# probably a clm/ini file. Mask maskc point (pygridgen complex grid):
      if 'maskc' in netcdf.varnames(self.grid.nc):
        mc=self.grid.use('maskc')
        if inds: mc=mc[j0:j1,i0:i1]
        h=np.ma.masked_where(mc==0,h)
        zeta=np.ma.masked_where(mc==0,zeta)

    h    = calc.griddata(xr,yr,h,x,y,extrap=False)
    zeta = calc.griddata(xr,yr,zeta,x,y,extrap=False)

    z=rt.s_levels(h,zeta,self.s_params,rw=rw)
    z=np.squeeze(z)
    return np.ma.masked_where(np.isnan(z),z)


###  def s_levels(self,time,ruvpw='r',i=False,j=False,k=False,extrapZeta=False):
  def s_levels(self,time,loc='rr',i=False,j=False,k=False,extrapZeta=False):
##    ruvpw=ruvpw[0]
    hLoc,vLoc=loc

    h=self.grid.h
    zeta=self.use('zeta',SEARCHtime=time)

    if extrapZeta:
      if not calc.ismarray(zeta): zeta=np.ma.masked_where(self.grid.mask==0,zeta)
      zeta=calc.mask_extrap(self.grid.lon,self.grid.lat,zeta)

    h=rt.rho2uvp(h,hLoc) ##########ruvpw)
    zeta=rt.rho2uvp(zeta,hLoc) ###########ruvpw)

    z=rt.s_levels(h,zeta,self.s_params,rw=vLoc) ##########ruvpw)

    if k is False: k=slice(None)
    if j is False: j=slice(None)
    if i is False: i=slice(None)

    return z[k,j,i]

# --------

  def check_slice(self,varname,**dims):
    msg=''

    # check varname:
    if varname not in netcdf.varnames(self.nc):
      return ':: variable %s not found' % varname

    isU=varname in ('u','ubar')
    isV=varname in ('v','vbar')
    isW=varname=='w'

    for dim,ind in dims.items():
      # check time:
      if dim=='t':
        if self.hast(varname) and ind>=self.TIME:
          msg='t = %d exceeds TIME dimension (%d)' % (ind,self.TIME)

      # check dim k:
      elif dim=='k':
        if self.hasz(varname):
          if isW: s='S_W'
          else: s='S_RHO'
          indMax=getattr(self,s)
          if ind >= indMax: msg='k = %d exceeds %s dimension (%d)' % (ind,s,indMax)

      # check dim i:
      elif dim=='i':
        xi='XI_RHO'
        if isU: xi='XI_U'
        elif isV and hasattr(self,'XI_V'): xi='XI_V' # needed for roms-agrif
        indMax=getattr(self,xi)
        if ind >= indMax:  msg='i = %d exceeds %s dimension (%d)' % (ind,xi,indMax)

      # check dim j:
      elif dim=='j':
        eta='ETA_RHO'
        if isU and hasattr(self,'ETA_U'): eta='ETA_U' # needed for roms-agrif
        elif isV: eta='ETA_V'
        indMax=getattr(self,eta)
        if ind >= indMax: msg='j = %d exceeds %s dimension (%d)' % (ind,eta,indMax)

    return msg

  @staticmethod
  def _default_coords(slc):
    if   slc=='slicei':      res='d,z,t'
    elif slc=='slicej':      res='d,z,t'
    elif slc=='slicek':      res='x,y,t'
    elif slc=='slicez':      res='x,y,t'
    elif slc=='slicell':     res='d,z,t'
    elif slc=='sliceuv':     res='x,y,t'
    elif slc=='sliceiso':    res='x,y,t'
    elif slc=='time_series': res='x,y,z,t'
##    return sorted(res)
    return ','.join(sorted(res.split(',')))

  def slicei(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicei')).split(',')

    out=Data() 
    out.msg=self.check_slice(varname,t=time,i=ind)
    if out.msg: return out
    
    v=self.use(varname,SEARCHtime=time,xi_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v):
      m=self.grid.vars(ruvp=self.var_at(varname)[0],i=ind)[-1]
      if v.ndim==2: m=np.tile(m,(v.shape[0],1))
      v=np.ma.masked_where(m==0,v)

    out.v=v
    out.info['v']['name']=varname
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
    if 'z' in coords and v.ndim==2:
#######   out.z=self.s_levels(time=time,ruvpw=self.var_at(varname),i=ind)
      out.z=self.s_levels(time=time,loc=self.var_at(varname),i=ind)
      out.info['z']=dict(name='Depth',units='m')

    if any([i in coords for i in 'xyd']):
      x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0],i=ind)

    if 'd' in coords:
      d=calc.distance(x,y)
      if d[-1]-d[0]>1e4:
        d=d/1000.
        dunits='km'
      else: dunits='m'

      if v.ndim==2: d=np.tile(d,(v.shape[0],1))
      out.d=d
      out.info['d']=dict(name='Distance',units=dunits)

    if 'x' in coords:
      if v.ndim==2: x=np.tile(x,(v.shape[0],1))
      out.x=x
      out.info['x']=dict(name='Longitude',units=r'$\^o$E')

    if 'y' in coords:
      if v.ndim==2: y=np.tile(y,(v.shape[0],1))
      out.y=y
      out.info['y']=dict(name='Latitude',units=r'$\^o$N')

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[Data()]
      if 'd' in coords: out.extra[0].x=out.d[0]
      if 'y' in coords: out.extra[0].x=out.y[0]
      if 'x' in coords: out.extra[0].y=out.x[0]
      out.extra[0].v=-h
      out.extra[0].config['d1.plot']='fill_between'
      out.extra[0].config['d1.y0']=-h.max()-(h.max()-h.min())/20.
      out.extra[0].alias='bottom'

    out.coordsReq=','.join(sorted(coords))
    return out 


  def slicej(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicej')).split(',')

    out=Data()
    out.msg=self.check_slice(varname,t=time,j=ind)
    if out.msg: return out

    v=self.use(varname,SEARCHtime=time,eta_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v):
      m=self.grid.vars(ruvp=self.var_at(varname)[0],j=ind)[-1]
      if v.ndim==2: m=np.tile(m,(v.shape[0],1))
      v=np.ma.masked_where(m==0,v)

    out.v=v
    out.info['v']['name']=varname
    out.info['v']['slice']='j=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.name,varname,'units')
    except: pass

    # coords:
    if 'z' in coords and v.ndim==2:
######      out.z=self.s_levels(time=time,ruvpw=self.var_at(varname),j=ind)
      out.z=self.s_levels(time=time,loc=self.var_at(varname),j=ind)
      out.info['z']=dict(name='Depth',units='m')

    if any([i in coords for i in 'xyd']):
      x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0],j=ind)

    if 'd' in coords:
      d=calc.distance(x,y)
      if d[-1]-d[0]>1e4:
        d=d/1000.
        dunits='km'
      else: dunits='m'

      if v.ndim==2: d=np.tile(d,(v.shape[0],1))
      out.d=d
      out.info['d']=dict(name='Distance',units=dunits)

    if 'x' in coords:
      if v.ndim==2: x=np.tile(x,(v.shape[0],1))
      out.x=x
      out.info['x']=dict(name='Longitude',units=r'$\^o$E')

    if 'y' in coords:
      if v.ndim==2: y=np.tile(y,(v.shape[0],1))
      out.y=y
      out.info['y']=dict(name='Latitude',units=r'$\^o$N')

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[Data()]
      if 'd' in coords: out.extra[0].x=out.d[0]
      if 'x' in coords: out.extra[0].y=out.x[0]
      if 'y' in coords: out.extra[0].x=out.y[0]
      out.extra[0].v=-h
      out.extra[0].config['d1.plot']='fill_between'
      out.extra[0].config['d1.y0']=-h.max()-(h.max()-h.min())/20.
      out.extra[0].alias='bottom'

    out.coordsReq=','.join(sorted(coords))
    return out


  def slicek(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicek')).split(',')

    out=Data()
    out.msg=self.check_slice(varname,t=time,k=ind)
    if out.msg: out

    v=self.use(varname,SEARCHtime=time,s_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v): 
      m=self.grid.vars(ruvp=self.var_at(varname)[0])[-1]
      v=np.ma.masked_where(m==0,v)

    out.v=v
    out.info['v']['name']=varname
    if self.hasz(varname): out.info['v']['slice']='k=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
    if 'z' in coords and self.hasz(varname):
#####      out.z=self.s_levels(time,k=ind,ruvpw=self.var_at(varname))
      out.z=self.s_levels(time,k=ind,loc=self.var_at(varname))
      out.info['z']=dict(name='Depth',units='m')


    if any([i in coords for i in 'xy']):
      x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0])

    if 'x' in coords:
       if self.grid.is_spherical:
         out.x=x
         out.info['x']=dict(name='Longitude',units=r'$\^o$E')
       else:
         out.x=x/1000.
         out.info['x']=dict(name='Distance',units='km')
         
    if 'y' in coords:
       if self.grid.is_spherical:
         out.y=y
         out.info['y']=dict(name='Latitude',units=r'$\^o$N')
       else:
         out.y=y/1000.
         out.info['y']=dict(name='Distance',units='km')

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.linecolors']='k'
      out.extra[0].alias='bathy'


    out.coordsReq=','.join(sorted(coords))
    return out


  def slicez(self,varname,ind,time=0,**opts):
    if not self.hasz(varname):
      return self.slicek(varname,ind,time,**opts)

    surf_nans=opts.get('surf_nans',True)
    spline=opts.get('spline',True)
    coords=opts.get('coords',self._default_coords('slicez')).split(',')

    out=Data()
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out##None,au

    v=self.use(varname,SEARCHtime=time)
    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0])
    zeta=self.use('zeta',SEARCHtime=time)
    zeta=rt.rho2uvp(zeta,varname)

    out.v=rt.slicez(v,m,h,zeta,self.s_params,ind,surf_nans,spline)

    out.info['v']['name']=varname
    if calc.isarray(ind):
      out.info['v']['slice']='z= array %.2f to %.2f'%(ind.min(),ind.max())
    else:
      out.info['v']['slice']='z=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
    if 'x' in coords:
       if self.grid.is_spherical:
         out.x=x
         out.info['x']=dict(name='Longitude',units=r'$\^o$E')
       else:
         out.x=x/1000.
         out.info['x']=dict(name='Distance',units='km')

    if 'y' in coords:
       if self.grid.is_spherical:
         out.y=y
         out.info['y']=dict(name='Latitude',units=r'$\^o$N')
       else:
         out.y=y/1000.
         out.info['y']=dict(name='Distance',units='km')

    if 'z' in coords:
      out.z=ind+np.zeros(out.v.shape)
      out.info['z']=dict(name='Depth',units='m')

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.linecolors']='k'
      out.extra[0].alias='bathy'


    out.coordsReq=','.join(sorted(coords))
    return out


  def slicell(self,varname,X,Y,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicell')).split(',')

    data      = opts.get('data',False)
    extrap    = opts.get('extrap',False)
    maskLimit = opts.get('lmask',0.5) # points where interpolated mask are above
                                      # this value are considered as mask!
                                      # Most strict value is 0

    out=Data()
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out#None,aux

    X=np.asarray(X)
    Y=np.asarray(Y)
    if X.ndim>1: X=np.squeeze(X)
    if Y.ndim>1: Y=np.squeeze(X)

    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0])
    if True: # extrat only portion of data needed:
      i0,i1,j0,j1=calc.ij_limits(x, y, (X.min(),X.max()),(Y.min(),Y.max()), margin=1)
      xi='%d:%d'%(i0,i1)
      eta='%d:%d'%(j0,j1)

      if data is False: V=self.use(varname,SEARCHtime=time,xi_SEARCH=xi,eta_SEARCH=eta)
      else: v=data[...,j0:j1,i0:i1]

      x=x[j0:j1,i0:i1]
      y=y[j0:j1,i0:i1]
      #h=h[j0:j1,i0:i1] # not used
      m=m[j0:j1,i0:i1]

    else:
      if data is False: V=self.use(varname,SEARCHtime=time)
      else: v=data

    if V.ndim==3:
      v=calc.griddata(x,y,V,X,Y,extrap=extrap,mask2d=m==0, keepMaskVal=maskLimit)
    elif V.ndim==2:
      v=calc.griddata(x,y,np.ma.masked_where(m==0,V),X,Y,extrap=extrap, keepMaskVal=maskLimit)

    out.v=v
    out.info['v']['name']=varname
    out.info['v']['slice']='path npts=%d'%X.size
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
    if 'z' in coords and V.ndim==3:
      inds=dict(xi=(i0,i1),eta=(j0,j1))
#########      out.z=self.path_s_levels(time,X,Y,rw=varname[0],inds=inds)
      out.z=self.path_s_levels(time,X,Y,rw=self.var_at(varname)[1],inds=inds)

    if 'd' in coords:
      d=calc.distance(X,Y)
      if v.ndim==2: d=np.tile(d,(v.shape[0],1))
      out.d=d

    if 'x' in coords:
      if v.ndim==2: X=np.tile(X,(v.shape[0],1))
      out.x=X

    if 'y' in coords:
      if v.ndim==2: Y=np.tile(Y,(v.shape[0],1))
      out.y=Y

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    out.coordsReq=','.join(sorted(coords))
    return out


  def sliceuv(self,ind,time=0,**opts):#plot=False,**opts):
###    plot= opts.pop('plot',False)
###    opts['plot']=False
###    ax  = opts.get('ax',None)
    coords=opts.get('coords',self._default_coords('sliceuv')).split(',')

####    out=Data()
####
#    savename=opts.pop('savename',False)
#    retMsg=opts.pop('msg',False)

    if ind=='bar':
      slc,uname,vname,ind = self.slicek, 'ubar','vbar', 9999
    elif ind in ('s','surf','surface'):
      slc,uname,vname,ind = self.slicek,  'u','v', -1
    elif ind>=0:
      slc,uname,vname,ind = self.slicek,  'u','v', ind
    elif ind <0:
      isK=False
      slc,uname,vname,ind = self.slicez, 'u','v', ind


##    if isK:
#      xu,yu,zu,u,msg1=self.slicek(uname,ind,time,msg=True,**opts)
#      xv,yv,zv,v,msg2=self.slicek(vname,ind,time,msg=True,**opts)
##    u,auxu=slc(uname,ind,time,**opts)
##    v,auxv=slc(vname,ind,time,**opts)
    outu=slc(uname,ind,time,**opts)
    outv=slc(vname,ind,time,**opts)

    if   outu.msg: return outu##None,None,auxu
    elif outv.msg: return outv##None,None,auxv

    # at psi:
    u,v=outu.v,outv.v
    u=( u[1:,:]+ u[:-1,:])/2.
    v=( v[:,1:]+ v[:,:-1])/2.

    # rotate uv:
    ang=rt.rho2uvp(self.grid.angle,'p')
    u,v=calc.rot2d(u,v,-ang)


##    else:
##      xu,yu,zu,u,msg1=self.slicez(uname,ind,time,msg=True,**opts)
##      xv,yv,zv,v,msg2=self.slicez(vname,ind,time,msg=True,**opts)
##
##    if msg1: msg=msg1+' and '+msg2
##    else: msg=msg2

    out=outu
    out.v=u,v
    out.info['v']['name']=uname+vname

    if 'z' in coords: out.z=(out.z[1:,:]+out.z[:-1,:])/2.
    if 'x' in coords: out.x=(out.x[1:,:]+out.x[:-1,:])/2.
    if 'y' in coords: out.y=(out.y[1:,:]+out.y[:-1,:])/2.

    out.coordsReq=','.join(sorted(coords))
    return out###u,v,aux


##    # coords:
##    if 'z' in coords:
##      out.z=(outu.z[1:,:]+outu.z[:-1,:])/2.
##
##    if 'x' in coords: out.x=self.grid.lonp
##    if 'y' in coords: out.y=self.grid.latp
##    if 't' in coords: out.t=self.time[time]
##
##
##    if plot:
##      if not ax:
##        ax=pl.gca()
##        opts['ax']=ax
##
##      prj=self.grid.plot(bathy=None,**opts)
##      xm, ym = prj(out.x,out.y)
##
##      mm=(0*xm).astype('bool')
##      nMax=150
##      nj,ni=xm.shape
##      dj=np.ceil(nj/float(nMax))
##      di=np.ceil(ni/float(nMax))
##      mm[::dj,::di]=1
######      mm=mm*self.grid.maskp==1
##
##      s=np.sqrt(u**2+v**2)
##      q=ax.quiver(xm[mm],ym[mm],u[mm],v[mm],s[mm])
##      pl.colorbar(q,shrink=.7)
##      ax.axis([prj.xmin, prj.xmax, prj.ymin, prj.ymax])
##      qk = ax.quiverkey(q, .05, .01, 0.2, '0.2 m/s',labelpos='E',
##                                              coordinates='axes')
####      if savename:
####        pylab.savefig(savename,dpi=pylab.gcf().dpi)
####        pylab.close(pylab.gcf())
##
####    if retMsg: return x,y,z,u,v, msg
####    else: return x,y,z,u,v
##
##    out.v=u,v
##    out.coordsReq=','.join(sorted(coords))
##    return out###u,v,aux


  def sliceiso(self,varname,iso,time,**opts):
    '''
    Depths where variable, increasing/dec with depth, has some value.

    Output is masked where surface is higher/lower than value or where all
    water column is lower/higher than value.
    '''


    coords=opts.get('coords',self._default_coords('sliceiso')).split(',')

    out=Data()
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out

    if not self.hasz(varname):
      out.msg='a depth dependent variable is needed for sliceiso!'
      return out

    v=self.use(varname,SEARCHtime=time)
####    z=self.s_levels(time=time,ruvpw=self.var_at(varname))
    z=self.s_levels(time=time,loc=self.var_at(varname))
    v=rt.depthof(v,z,iso)

    out.v=v
    #out.info['v']['name']=varname
    out.info['v']['name']='depth'

    if isinstance(iso,np.ndarray) and iso.ndim==2: siso='2d'
    else: siso=str(iso)

    try:
      out.info['v']['slice']='%s (%s) iso=%s'%(varname,netcdf.vatt(self.nc,varname,'units'),siso)
    except:
      out.info['v']['slice']='%s iso=%s'%(varname,siso)
  
    out.info['v']['units']='metre'


    # coords:
    if any([i in coords for i in 'xy']):
      x,y,h=self.grid.vars(ruvp=self.var_at(varname))[:3]

    if 'x' in coords:
       if self.grid.is_spherical:
         out.x=x
         out.info['x']=dict(name='Longitude',units=r'$\^o$E')
       else:
         out.x=x/1000.
         out.info['x']=dict(name='Distance',units='km')

    if 'y' in coords:
       if self.grid.is_spherical:
         out.y=y
         out.info['y']=dict(name='Latitude',units=r'$\^o$N')
       else:
         out.y=y/1000.
         out.info['y']=dict(name='Distance',units='km')

    if 'z' in coords:
      # makes no sense... v is the depth!
      coords.remove('z')

    if 't' in coords and self.hast(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.linecolors']='k'
      out.extra[0].alias='bathy'

    out.coordsReq=','.join(sorted(coords))
    return out


  def time_series(self,varname,x,y,times=None,depth=None,**opts):
    coords=opts.get('coords',self._default_coords('time_series')).split(',')

    if times is None: times=range(0,self.time.size)

    # depth or s_level: check is is float or if is negative!
    isDepth=False
    if not depth is None:
       if calc.isiterable(depth): depth=np.asarray(depth)
       if calc.isarray(depth):
         isDepth=np.any(depth<0) or depth.kind!='i' 
       else: isDepth=depth<0 or np.asarray(depth).dtype.kind!='i'

    out=Data()
    if not depth is None and not isDepth:
      out.msg=self.check_slice(varname,t=np.max(times),k=depth) 
    else:
      out.msg=self.check_slice(varname,t=np.max(times)) 

    if out.msg: return out

    # find nearest point:
    lon,lat,hr,mr=self.grid.vars(ruvp=self.var_at(varname))
    dist=(lon-x)**2+(lat-y)**2
    i,j=np.where(dist==dist.min())
    i,j=i[0],j[0]

    if not depth is None and not isDepth: arg={'s_SEARCH':depth}
    else: arg={}
    v=self.use(varname,xiSEARCH=j,etaSEARCH=i,SEARCHtime=times,**arg).T

    # calculate depths:
    if self.hasz(varname):
      h=self.grid.h[i,j]
      zeta=self.use('zeta',xiSEARCH=j,etaSEARCH=i,SEARCHtime=times)
      h=h+0*zeta
####      z=rt.s_levels(h,zeta,self.s_params,rw=varname)
      z=rt.s_levels(h,zeta,self.s_params,rw=self.var_at(varname)[1])
      z=np.squeeze(z)

    # depth slice:
    if isDepth and self.hasz(varname):
      if v.ndim==2:
        # could use calc.griddata, but better use slicez cos data at
        # different times may be independent!
        if 0:
          from matplotlib.dates import date2num
          t=np.tile(date2num(self.time[times]),(v.shape[0],1))
          v=calc.griddata(t,z,v,t[0],depth+np.zeros(t[0].shape),
            extrap=opts.get('extrap',False),norm_xy=opts.get('norm_xy',False))
            # norm_xy True may be needed!
            # extrap also may be needed cos the 1st and last value may be masked!

        else:
          nt=len(times)
          land_mask=np.ones((nt,1),dtype=v.dtype) # needed for slicez... not used here!        

          v,vm=rt.slicez(v[...,np.newaxis],land_mask,
               self.grid.h[i,j]*np.ones((nt,1),dtype=v.dtype), # bottom depth
               zeta[:,np.newaxis],self.s_params,depth,spline=opts.get('spline',True))

          v=np.ma.masked_where(vm,v)
          v=v[...,0]

      else: # one time only
        v=np.interp(depth,z,v,left=np.nan,right=np.nan)
        v=np.ma.masked_where(np.isnan(v),v) 


    out.v=v
    out.info['v']['name']=varname
    out.info['v']['slice']='time series'
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass
     
    # coords
    if 't' in coords and self.hast(varname):
      if v.ndim==2:
        out.t=np.tile(self.time[times],(v.shape[0],1))
        from matplotlib.dates import date2num
        out.tnum=np.tile(date2num(self.time[times]),(v.shape[0],1))
      else: out.t=self.time[times]
      out.info['t']['name']='Time'
      out.info['tnum']=dict(name='Time',units=self.var_as['time']['units'])

    if 'z' in coords and self.hasz(varname):
      if not depth is None:
        if not isDepth: out.z=z[depth,...]
        else: out.z=depth+0*v
      else: out.z=z
      out.info['z']=dict(name='Depth',units='m')

    if 'x' in coords:
      out.x=lon[i,j]
      if self.grid.is_spherical:
         out.info['x']=dict(name='Longitude',units=r'$\^o$E')
      else:
        out.x=x/1000.
        out.info['x']=dict(name='X-position',units='km')

    if 'y' in coords:
      out.y=lat[i,j]
      if self.grid.is_spherical:
        out.info['y']=dict(name='Latitude',units=r'$\^o$N')
      else:
        out.y=y/1000.
        out.info['y']=dict(name='Y-position',units='km')


    out.coordsReq=','.join(sorted(coords))
    return out


  def extrap_at_mask(self,vname,time,**opts): 
    method=opts.get('method','delaunay') # also available 'easy'
    quiet=opts.get('quiet',True)

    if method=='delaunay': extrap=calc.mask_extrap
    elif method=='easy': extrap=calc.easy_extrap

    if not quiet: print 'Extrap at mask:\n  loading'
    v=self.use(vname,SEARCHtime=time)
    x,y,h,m=self.grid.vars(vname)

    if self.hasz(vname): # 3d
      N=v.shape[0]
      for n in range(N):
        if not quiet: print '  extrap level %d of %d' % (n,N)
        v[n,...]=extrap(x,y,np.ma.masked_where(m==0,v[n,...]))
    else: # 2d
        v=extrap(x,y,np.ma.masked_where(m==0,v))

    if np.ma.isMA(v) and v.count()==v.size: v=v.data # no need for masked array
    return v

# ------------------------------------------------------------------------ OLD CODE FROM HERE
# may not be working anymore....

class Flt(Common):
  def __init__(self,flt,grd=False,quiet=True):
    self.name=os.path.realpath(flt)
    self.type='flt'
    self.quiet=quiet

    self.load_grid(grd)
    self.load()
    self.load_dims()
    self.load_atts()


  def load(self):
    # time:
    vars={'time':('time','scrum_time','ocean_time')},'lon','lat','depth'

    self.load_vars(vars)

    if len(self.time)>1:
      self.dt=self.time[1]-self.time[0]
    else: self.dt=0

    self.tdays=self.time/86400.


  def order_dist(self,quiet=True):
    '''
    order floats by distance travelled
    returns floats indices and distances, ie,
    returns dist.argsort(), dist
    '''
    d=zeros(self.DRIFTER)
    D=zeros(self.DRIFTER,'i')
    for i in range(self.DRIFTER):
      tmp=np.where(self.lon[:,i]>1000)[0]
      if len(tmp):ie=tmp[0]-1
      else: ie=-1

      d[i]=distance(self.lon[[0,ie],i],self.lat[[0,ie],i])[1]

    D=d.argsort()

    if not quiet:
      for i in range(len(D)):
        print '%5d %8.2f' % (D[i],d[D[i]])

    return D[::-1],d[::-1]


class Sta(Common):
  def __init__(self,sta,grd=False,quiet=True):
    self.name=os.path.realpath(sta)
    self.type='sta'
    self.quiet=quiet

    self.load_grid(grd)
    self.load()
    self.load_dims()
    self.load_atts()

  def load(self):
    # time:
    vars={'time':('time','scrum_time','ocean_time')},'lon','lat','h'

    self.load_vars(vars)

    if len(self.time)>1:
      self.dt=self.time[1]-self.time[0]
    else: self.dt=0

    self.tdays=self.time/86400.

    # s params:
    self.s_params=s_params(self.name)


  def s_levels(self,itime,isw=False):
    tts,ttb,hc,N=self.s_params
    zeta=self.use('zeta',ftime=itime)

    zr,zw=rt.s_levels(self.h,zeta,hc,tts,ttb,N)
    if isw:z=zw
    else: z=zr

    return squeeze(z)

  def get(self,what,**kargs):
    time=False
    day=False
    quiet=True

    keys=kargs.keys()
    if 'time'  in keys: time  = kargs['time']
    if 'day'   in keys: day   = kargs['day']
    if 'quiet' in keys: quiet = kargs['quiet']


    if what=='data':
      out={}
      args={}
      if not time is False:
        args['ftime']=time

      if not quiet:
        print 'loading Sta data time=',time
        print '* * * * * * * * * * * *.'

      if not quiet: print '*',
      out['lon']  = self.use('lon',  **args)
      if not quiet: print '*',
      out['lat']  = self.use('lat',  **args)
      if not quiet: print '*',
      out['h']  = self.use('h',  **args)
      if not quiet: print '*',
      out['zeta'] = self.use('zeta', **args)
      if not quiet: print '*',
      out['ubar'] = self.use('ubar', **args)
      if not quiet: print '*',
      out['vbar'] = self.use('vbar', **args)

    statime=self.use('scrum_time',**args)
    if 't0' in self.atts.keys():
      t0=self.atts['t0']
    else: t0=0

    return out,statime,t0



class Blk(Common):
  def __init__(self,flt,grd=False,quiet=True):
    self.name=os.path.realpath(flt)
    self.type='blk'
    self.quiet=quiet

    self.load_grid(grd)
    self.load()
    self.load_dims()
    self.load_atts()

  def load(self):
    # time:
    vars={'time':('bulk_time','time','scrum_time','ocean_time')},

    self.load_vars(vars)

    # netcdf time:
    try:
      # units must have the format <units> since <date>
      self.time=netcdf.num2date(self.time,self.var_as['time']['units'])
    except: pass
#    try:
#      # units must have the format <units> since <date>
#      self.datetime=netcdf.num2date(self.time,self.var_as['time']['units'])
#    except: self.datetime=False
#
#    # time is usually seconds, but may be days!!, so:
#    if self.var_as['time']['units'].strip().startswith('days'): self.time=self.time*86400.
#
#    if len(self.time)>1:
#      self.dt=self.time[1]-self.time[0]
#    else: self.dt=0
#
#    self.tdays=self.time/86400.


  def get(self,what,**kargs):
    time=False
    day=False
    quiet=True

    keys=kargs.keys()
    if 'time'  in keys: time  = kargs['time']
    if 'day'   in keys: day   = kargs['day']
    if 'quiet' in keys: quiet = kargs['quiet']


    if what=='wind_ts':
      if 'xi' in kargs.keys() and 'eta' in kargs.keys():
        xi=kargs['xi']
        eta=kargs['eta']
      elif 'lon' in kargs.keys() and 'lat' in kargs.keys():
        d=(self.grid.lon-kargs['lon'])**2+(self.grid.lat-kargs['lat'])**2
        i,j=np.where(d==d.min())
        eta=i[0]
        xi=j[0]

      u=self.use('uwnd',xi_rho=xi,eta_rho=eta)
      v=self.use('vwnd',xi_rho=xi,eta_rho=eta)
      ang=self.grid.angle[eta,xi]*np.ones(u.shape)
      u,v=calc.rot2d(u,v,-ang)
      tdays=self.tdays

      if not day is False:
        ndays=self.tdays[-1]-self.tdays[0]
        n=int(86400/self.dt)
        u=u[day*n:(day+1)*n]
        v=v[day*n:(day+1)*n]
        tdays=tdays[day*n:(day+1)*n]

      return tdays,u,v


    elif what=='wind':
      if 'uwnd' in netcdf.varnames(self.name):
        u=self.use('uwnd')
        v=self.use('vwnd')
      else:
        u=self.use('Uwind')
        v=self.use('Vwind')

      if time=='mean' and u.ndim==3:
        u=u.mean(0)
        v=v.mean(0)
      elif time=='dailyMean' and u.ndim==3:
        ndays=self.tdays[-1]-self.tdays[0]
        n=int(86400/self.dt)

        if not day is False:
          if day==-1: day=ndays-1
          # this is wrong if there is data missing!
          u=u[day*n:(day+1)*n,...].mean(0)
          v=v[day*n:(day+1)*n,...].mean(0)
        else:

          U=range(ndays)
          V=range(ndays)
          for i in range(ndays):
            U[i]=u[i*n:(i+1)*n,...].mean(0)
            V[i]=v[i*n:(i+1)*n,...].mean(0)

          u=U
          v=V

      elif not time is False and isinstance(time,int) and u.ndim==3:
        u=u[time,...]
        v=v[time,...]

      if isinstance(u,list):
        for i in range(len(u)):
          u[i],v[i]=calc.rot2d(u[i],v[i],-self.grid.angle)
      else:
        if u.ndim==3:
          for i in range(u.shape[0]):
            u[i,...],v[i,...]=calc.rot2d(u[i,...],v[i,...],-self.grid.angle)
        elif u.ndim==2:
          u,v=calc.rot2d(u,v,-self.grid.angle)

      if self.grid:
        x=self.grid.lon
        y=self.grid.lat
        m=self.grid.mask
        return x,y,u,v,m
      else:
        return u,v

    elif what=='data':
      out={}
      args={}
      if not time is False:
        args['bulk_time']=time

      if not quiet:
        print 'loading Blk data time=',time
        print '* * * * * * * * * * * *.'

      if not quiet: print '*',
      out['tair']  = self.use('tair',  **args)
      if not quiet: print '*',
      out['rhum']  = self.use('rhum',  **args)
      if not quiet: print '*',
      out['pres']  = self.use('pres',  **args)
      if not quiet: print '*',
      out['prate'] = self.use('prate', **args)
      if not quiet: print '*',
      out['radsw'] = self.use('radsw', **args)
      if not quiet: print '*',
      out['radlw'] = self.use('radlw', **args)
      if not quiet: print '*',
      out['dlwrf'] = self.use('dlwrf', **args)
      if not quiet: print '*',
      out['uwnd']  = self.use('uwnd',  **args)
      if not quiet: print '*',
      out['vwnd']  = self.use('vwnd',  **args)
      if not quiet: print '*',
      out['sustr'] = self.use('sustr', **args)
      if not quiet: print '*',
      out['svstr'] = self.use('svstr', **args)
      if not quiet: print '*.'
      out['wspd']  = self.use('wspd',  **args)

    blktime=self.use('bulk_time',**args)
    t0=self.atts['t0']

    return out,blktime,t0


