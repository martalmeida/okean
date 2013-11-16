import numpy as np
import pylab
import os
import netcdftime

from okean import calc, netcdf, ticks
import roms_tools as rt
from derived import Derived

try:
  from mpl_toolkits.basemap import Basemap
except:
  Basemap=False

__all__='His','Grid'

class Common:
  def load_vars(self,vars):
    var,nc=netcdf.var(self.name)
    varnames=var.keys()
    for i in vars:
      if not self.quiet: print ':: loading var ',  i

      if isinstance(i,dict):
        iname, iopts = i.items()[0]
      elif isinstance(i,basestring):
        iname=i
        iopts=i
      elif isinstance(i,tuple):
        iname=False
        iopts=i

      if not isinstance(iopts,tuple):
        iopts=(iopts,)

      found=0
      if not hasattr(self,'var_as'): self.var_as={}
      for j in iopts:
        if not iname: iname=j
        if j in varnames:
          setattr(self,iname,var[j][:])
          found=1
          if not self.quiet: print '  - found %s as %s'  % (iname,j)

          # store original name and other info, like units
          try: units=var[j].atts['units'].value
          except: units=''
          self.var_as[iname]={'name':j,'units':units}
          break
      if not found and not self.quiet: print ':: %s not found in %s' % (str(iopts),self.name)

    nc.close()


  def load_dims(self):
    dms=netcdf.fdim(self.name)
    for k in dms.keys():
      if k in ('ocean_time','time','scrum_time'):
        setattr(self,'TIME',dms[k])
      else:
        setattr(self,k.upper(),dms[k])

  def load_atts(self):
    atts=netcdf.fatt(self.name)
    self.atts={}
    for k in atts.keys(): self.atts[k]=atts[k].value

  def load_grid(self,grd=False):
    ats=netcdf.fatt(self.name)
    if not grd:
      if 'grd_file' in ats.keys(): grd=ats['grd_file'].value
      if grd and not os.path.isfile(grd): grd=False

    if grd: self.grid=Grid(grd,self.quiet)
    else:
      try:
        self.grid=Grid(self.name)
      except:
        self.grid=False

  def show(self):
    netcdf.show(self.name)

  def use(self,varname,**kargs):
    return netcdf.use(self.name,varname,**kargs)

  def hasz(self,v):
    '''
    True if variable has depth dimension (s_rho or s_w)
    '''

    dnames=netcdf.vdim(self.name,v)
    return 's_rho' in dnames or 's_w' in dnames

  def hast(self,v):
    '''
    True if variable has *time* dimension
    '''

    dnames=netcdf.vdim(self.name,v)
    for d in dnames:
      if d.find('time')>=0: return True

    return False

  def var_at(self,v,_3d=True):
    '''
    returns location of the variable in the 2d/3d grid
    Possible values are u,v,r(rho), p(psi)
    If _3d, also checks if at w points
    '''
    dims=netcdf.vdim(self.name,v).keys()

    if _3d and 's_w' in dims: return 'w'
    elif 'xi_u'   in dims or   'eta_u'   in dims: return 'u'
    elif 'xi_v'   in dims or   'eta_v'   in dims: return 'v'
    elif 'xi_psi' in dims and  'eta_psi' in dims: return 'p'
    elif 'xi_rho' in dims and  'eta_rho' in dims: return 'r'

    return False


class Grid(Common):
  def __init__(self,grd,quiet=True):
    if grd.startswith('http'): self.name=grd
    else: self.name=os.path.realpath(grd)
    self.type='grid'
    self.quiet=quiet

    self.load()
    self.load_dims()
    self.load_uvp()
    self.load_atts()

  def load(self):
    vars={'lon':('lon_rho','x_rho')}, {'lonu':('lon_u','x_u')},\
         {'lonv':('lon_v','x_v')},    {'lonp':('lon_psi','x_psi')},\
         {'lat':('lat_rho','y_rho')}, {'latu':('lat_u','y_u')},\
         {'latv':('lat_v','y_v')},    {'latp':('lat_psi','y_psi')},\
         {'mask':'mask_rho'},{'masku':'mask_u'},{'maskv':'mask_v'},{'maskp':'mask_psi'},\
         'h','angle', 'pm','pn'
    self.load_vars(vars)

    # some vars may not be present... set defaults:
    if not hasattr(self,'mask'): self.mask=np.ones(self.h.shape)
    if not hasattr(self,'angle'): self.angle=np.zeros(self.h.shape)


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

  def s_levels(self,sparams,zeta=0,h=False,ruvpw='r',i=False,j=False,k=False):
    ruvpw=ruvpw[0]
    isW=ruvpw=='w'

    if h is False:
      h=self.h

    try:
      zeta.shape==h.shape
    except:
      zeta=np.tile(zeta,h.shape).astype(h.dtype)

    if h.ndim==2:
      h=rt.rho2uvp(h,ruvpw)
      zeta=rt.rho2uvp(zeta,ruvpw)

    zr,zw=rt.s_levels(h,zeta,sparams)
    if isW:z=zw
    else: z=zr

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


  def ingrid(self,x,y,retInds=False):
    '''in polygon of region border'''
    xb,yb=self.border()
    i=calc.inpolygon(x,y,xb,yb)
    if retInds:
      return x[i!=0], y[i!=0], np.where(i)[0]
    else: return x[i!=0], y[i!=0]


  def resolution(self):
    dx=1./self.pm
    dy=1./self.pn
    dx=np.ma.masked_where(self.mask==0,dx)
    dy=np.ma.masked_where(self.mask==0,dy)
    return self.lon,self.lat,dx,dy


  def plot(self,**kargs):
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
    cmap       = kargs.get('cmap',pylab.cm.gist_earth_r)
    proj       = kargs.get('proj','merc')
    ax         = kargs.get('ax',False)

    if not ax:
      fig=pylab.figure()
      ax=pylab.axes()
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
      m = Basemap(projection=proj,lat_ts=0.0,
                  #lon_0=self.lon[0].mean(),
                  resolution=resolution,
                  urcrnrlon=Lonlims[1], urcrnrlat=Latlims[1],
                  llcrnrlon=Lonlims[0], llcrnrlat=Latlims[0])

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
        if m.xmax-m.xmin<m.ymax-m.ymin: orientation='horizontal'
        else: orientation='vetical'

        cbh = pylab.colorbar(pch,orientation=orientation,ax=ax)
        cbh.set_label('Depth')
        if title=='auto': ax.set_title(self.name)
        elif title: ax.set_title(title)

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

    # netcdf time:
    try:
      # units must have the format <units> since <date>
      self.datetime=netcdftime.num2date(self.time,self.var_as['time']['units'])
    except: self.datetime=False


    # time is usually seconds, but may be days!!, so:
    if self.var_as['time']['units'].strip().startswith('days'): self.time=self.time*86400.

    if len(self.time)>1:
      self.dt=self.time[1]-self.time[0]
    else: self.dt=0

    self.tdays=self.time/86400.


    # s params:
    # hasattr(obj, '__contains__')
    if isinstance(self.name,basestring):
      self.s_params=rt.s_params(self.name)
    else: self.s_params=rt.s_params(self.name[0])



  def path_s_levels(self,time,x,y,rw='r'):
    xr,yr,h,m=self.grid.vars()
    zeta=self.use('zeta',SEARCHtime=time)
    h    = calc.griddata(xr,yr,h,x,y,extrap=False)
    zeta = calc.griddata(xr,yr,zeta,x,y,extrap=False)

    z=rt.s_levels(h,zeta,self.s_params,rw=rw)
    z=np.squeeze(z)
    return np.ma.masked_where(np.isnan(z),z)


  def s_levels(self,time,ruvpw='r',i=False,j=False,k=False,extrapZeta=False):
    ruvpw=ruvpw[0]

    h=self.grid.h
    zeta=self.use('zeta',SEARCHtime=time)

    if extrapZeta:
      if not calc.ismarray(zeta): zeta=np.ma.masked_where(self.grid.mask==0,zeta)
      zeta=calc.mask_extrap(self.grid.lon,self.grid.lat,zeta)

    h=rt.rho2uvp(h,ruvpw)
    zeta=rt.rho2uvp(zeta,ruvpw)

    z=rt.s_levels(h,zeta,self.s_params,rw=ruvpw)

    if k is False: k=slice(None)
    if j is False: j=slice(None)
    if i is False: i=slice(None)

    return z[k,j,i]

# --------

  def slicei(self,varname,ind,time=0,dist=1,plot=False,**opts):
    savename=False
    if 'savename' in opts.keys(): savename=opts['savename']

    d,x,y,z,v=[[]]*5

    if varname not in netcdf.varnames(self.name):
      print ':: variable %s not found' % varname
      if dist: return  d,z,v
      else:  return  x,y,z,v

    isU=varname=='u'
    isV=varname=='v'

    if isU:
      XI= self.XI_U
      xi='XI_U'
    elif isV:
      if hasattr(self,'XI_V'):
        XI= self.XI_V
        xi='XI_V'
      else: # this happens in current roms-agrif
        XI= self.XI_RHO
        xi='XI_RHO'
    else:
      XI= self.XI_RHO
      xi='XI_RHO'

    if ind >= XI:  print 'j = %d exceeds %s dimension (%d)' % (ind,xi,XI)

    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname),i=ind)
    karg={'SEARCHtime':time,xi.lower(): ind}
    v=self.use(varname,**karg)
    if v.ndim==2:
      z=self.s_levels(time=time,ruvpw=self.var_at(varname),i=ind)
    elif v.ndim==1:
      z=[]

    N=v.shape[0]
    if dist:
      d=calc.distance(x,y)
      if plot:

        pylab.figure()
        if v.ndim==2:
          pylab.pcolor(np.tile(d,(N,1)),z,np.ma.masked_where(np.tile(m,(N,1))==0,v),shading='flat')
          pylab.plot(d,-h)
          pylab.colorbar(shrink=.7)
        elif  v.ndim==1:
           pylab.plot(d,v)

        if savename:
          pylab.savefig(savename,dpi=pylab.gcf().dpi)
          pylab.close(pylab.gcf())

      if v.ndim==2:
        return np.tile(d,(N,1)),z,v #,np.tile(m,(N,1))==0
      elif v.ndim==1:
        return d,z,v #,m==0

    else:
      if plot:
        pylab.figure()
        if v.ndim==2:
          pylab.pcolor(np.tile(y,(N,1)),z,np.ma.masked_where(np.tile(m,(N,1))==0,v),shading='flat')
          pylab.colorbar(shrink=.7)
          pylab.plot(y,-h)
        elif  v.ndim==1:
          pylab.plot(y,v)

        if savename:
          pylab.savefig(savename,dpi=pylab.gcf().dpi)
          pylab.close(pylab.gcf())

      if v.ndim==2:
        return np.tile(x,(N,1)),np.tile(y,(N,1)),z,v #,np.tile(m,(N,1))==0
      elif v.ndim==1:
         return x,y,z,v #,m==0


  def slicej(self,varname,ind,time=0,dist=1,plot=False,**opts):
    savename=False
    if 'savename' in opts.keys(): savename=opts['savename']

    d,x,y,z,v=[[]]*5

    if varname not in netcdf.varnames(self.name):
      print ':: variable %s not found' % varname
      if dist: return  d,z,v
      else:  return  x,y,z,v

    isU=varname=='u'
    isV=varname=='v'

    if isU:
      if hasattr(self,'ETA_U'):
        ETA= self.ETA_U
        eta='ETA_U'
      else: # this happens in current roms-agrif
        ETA= self.ETA_RHO
        eta='ETA_RHO'
    elif isV:
      ETA= self.ETA_V
      eta='ETA_V'
    else:
      ETA= self.ETA_RHO
      eta='ETA_RHO'

    if ind >= ETA:
      print 'j = %d exceeds %s dimension (%d)' % (ind,eta,ETA)
      if dist: return  d,z,v
      else:  return  x,y,z,v

    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname),j=ind)
    karg={'SEARCHtime':time,eta.lower(): ind}
    v=self.use(varname,**karg)
    if v.ndim==2:
      z=self.s_levels(time=time,ruvpw=self.var_at(varname),j=ind)
    elif v.ndim==1:
      z=[]

    N=v.shape[0]
    if dist:
      d=calc.distance(x,y)
      if plot:
        pylab.figure()
        if v.ndim==2:
          pylab.pcolor(np.tile(d,(N,1)),z,np.ma.masked_where(np.tile(m,(N,1))==0,v),shading='flat')
          pylab.plot(d,-h)
          pylab.colorbar(shrink=.7)
        elif v.ndim==1:
          pylab.plot(d,v)

        if savename:
          pylab.savefig(savename,dpi=pylab.gcf().dpi)
          pylab.close(pylab.gcf())

      if v.ndim==2:
        return np.tile(d,(N,1)),z,v
      elif v.ndim==1:
         return d,z,v #np.ma.masked_where(m==0,v)

    else:
      if plot:
        pylab.figure()
        if v.ndim==2:
          pylab.pcolor(tile(x,(N,1)),z,ma.masked_where(tile(m,(N,1))==0,v),shading='flat')
          pylab.colorbar(shrink=.7)
          pylab.plot(x,-h)
        elif v.ndim==1:
          pylab.plot(x,v)

        if savename:
          pylab.savefig(savename,dpi=pylab.gcf().dpi)
          pylab.close(pylab.gcf())

      if v.ndim==2:
        m=np.tile(m,(N,1))==0
        return np.tile(x,(N,1)),np.tile(y,(N,1)),z,v #np.ma.masked_where(m,v)
      elif v.ndim==1:
         return x,y,z,v #np.ma.masked_where(m==0,v)



  def slicek(self,varname,ind,time=0,plot=False,**opts):
    x,y,z,v=[[]]*4

    RetVarOnly=False
    savename=False
    retMsg=False

    if ind in ('s','surf','surface'): ind=-1
    if 'retvar'   in opts.keys(): RetVarOnly = opts['retvar']
    if 'savename' in opts.keys(): savename   = opts['savename']
    if 'msg'      in opts.keys(): retMsg     = opts['msg']


    if varname not in netcdf.varnames(self.name):
      msg=':: variable %s not found' % varname
      if retMsg: return x,y,z,v,msg
      else:
        print msg
        return x,y,z,v

    isW=varname=='w'
    hasZ=self.hasz(varname)

    if hasZ:
      if isW and ind >= self.S_W:
        msg='k = %d exceeds S_W dimension (%d)' % (ind,self.S_W)
        if retMsg: return x,y,z,v,msg
        else:
          print msg
          return x,y,z,v

      elif ind >= self.S_RHO:
        msg='k = %d exceeds S_RHO dimension (%d)' % (ind,self.S_RHO)
        if retMsg: return x,y,z,v,msg
        else:
          print msg
          return x,y,z,v

    if self.hast(varname) and time>=self.TIME:
      msg='t = %d exceeds TIME dimension (%d)' % (time,self.TIME)
      if retMsg: return x,y,z,v,msg
      else:
        print msg
        return x,y,z,v

    if isW: v=self.use(varname,SEARCHtime=time,s_w=ind)
    else:   v=self.use(varname,SEARCHtime=time,s_rho=ind)

    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname))
    if hasZ:
      z=self.s_levels(time,k=ind,ruvpw=self.var_at(varname))
    else:
      z=self.use('zeta',SEARCHtime=time)
      z=rt.rho2uvp(z,varname)

    if not np.ma.isMA(v): # roms agrif ...
      v=np.ma.masked_where(m==0,v)

    if plot:
      p=self.grid.plot(bathy=None,**opts)
      xm, ym = p(x,y)
      pch=pylab.pcolor(xm,ym,v,shading='flat')
      pylab.colorbar(pch,shrink=.7)
      pylab.axis([p.xmin, p.xmax, p.ymin, p.ymax])
      if savename:
        pylab.savefig(savename,dpi=pylab.gcf().dpi)
        pylab.close(pylab.gcf())

    if RetVarOnly:
      return v
    else:
      if retMsg: return x,y,z,v,''
      else: return x,y,z,v


  def slicez(self,varname,ind,time=0,plot=False,**opts):
    x,y,z,v=[[]]*4

    savename  = False
    retMsg    = False
    surf_nans = True

    if 'savename'  in opts.keys(): savename  = opts['savename']
    if 'msg'       in opts.keys(): retMsg    = opts['msg']
    if 'surf_nans' in opts.keys(): surf_nans = opts['surf_nans']


    if varname not in netcdf.varnames(self.name):
      msg=':: variable %s not found' % varname
      if retMsg: return x,y,z,v,msg
      else:
        print msg
        return x,y,z,v

    if self.hast(varname) and time>=self.TIME:
      msg='t = %d exceeds TIME dimension (%d)' % (time,self.TIME)
      if retMsg: return x,y,z,v,msg
      else:
        print msg
        return x,y,z,v

    v=self.use(varname,SEARCHtime=time);
    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname))
    zeta=self.use('zeta',SEARCHtime=time)
    zeta=rt.rho2uvp(zeta,varname)

    if len(v.shape)==2: # no vertical comp
      if retMsg: return x,y,ind+np.zeros(v.shape),np.ma.masked_where(m==0,v),''
      else: return x,y,ind+np.zeros(v.shape),np.ma.masked_where(m==0,v)

    v,mask=rt.slicez(v,m,h,zeta,self.s_params,ind,surf_nans)
    v=np.ma.masked_where(mask,v)

    if plot:
      p=self.grid.plot(bathy=None,**opts)
      xm, ym = p(x,y)
      pch=pylab.pcolor(xm,ym,v,shading='flat')
      pylab.colorbar(pch,shrink=.7)
      pylab.axis([p.xmin, p.xmax, p.ymin, p.ymax])
      if savename:
        pylab.savefig(savename,dpi=pylab.gcf().dpi)
        pylab.close(pylab.gcf())

    if retMsg: return x,y,ind+np.zeros(v.shape),v,''
    else: return x,y,ind+np.zeros(v.shape),v


  def slicell(self,varname,X,Y,time=0,plot=False,**opts):
    x,y,z,v,m=[[]]*5

    data      = opts.get('data',False)
    dist      = opts.get('dist',False)
    extrap    = opts.get('extrap',False)
    varOnly   = opts.get('retv',False)
    maskLimit = opts.get('lmask',0.5) # points where interpolated mask are above
                                      # this value are considered as mask!
                                      # Most strict value is 0


    X=np.array(X)
    Y=np.array(Y)
    if X.ndim>1: X=np.squeeze(X)
    if Y.ndim>1: Y=np.squeeze(X)

    if varname not in netcdf.varnames(self.name):
      print ':: variable %s not found' % varname
      return x,y,z,v,m

    if time>=self.TIME:
      print 't = %d exceeds TIME dimension (%d)' % (time,self.TIME)
      return x,y,z,v,m

    if data is False: v=self.use(varname,SEARCHtime=time)
    else: v=data

    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname))

    if v.ndim==3:
      V=calc.griddata(x,y,v,X,Y,extrap=extrap,mask2d=m==0, keepMaskVal=maskLimit)
    elif v.ndim==2:
      V=calc.griddata(x,y,np.ma.masked_where(m==0,v),X,Y,extrap=extrap, keepMaskVal=maskLimit)

    if varOnly: return V

    # Z:
    if v.ndim==3:
      Z=self.path_s_levels(time,X,Y,rw=varname[0])
    else: Z=0.*X


    # X, Y, dist:
    if dist:
      Dist=calc.distance(X,Y)
      if v.ndim==3:
        Dist=np.tile(Dist,(v.shape[0],1))
    else:
      if v.ndim==3:
        X=np.tile(X,(v.shape[0],1))
        Y=np.tile(Y,(v.shape[0],1))

    if dist:
      return Dist,Z,V
    else:
      return X,Y,Z,V


  def sliceuv(self,ind,time=0,plot=False,**opts):
    savename=False
    if savename in opts.keys(): savename=opts['savename']

    if ind=='bar':
      isK,uname,vname,ind = True, 'ubar','vbar', 9999
    elif ind in ('s','surf','surface'):
      isK,uname,vname,ind = True,  'u','v', -1
    elif ind>=0:
      isK,uname,vname,ind = True,  'u','v', ind
    elif ind <0:
      isK=False
      isK,uname,vname,ind = False, 'u','v', ind


    if isK:
      xu,yu,zu,u=self.slicek(uname,ind,time)
      xv,yv,zv,v=self.slicek(vname,ind,time)
    else:
      xu,yu,zu,u=self.slicez(uname,ind,time)
      xv,yv,zv,v=self.slicez(vname,ind,time)

    z=(zu[1:,:]+zu[:-1,:])/2.
    u=( u[1:,:]+ u[:-1,:])/2.
    v=( v[:,1:]+ v[:,:-1])/2.

    x=self.grid.lonp
    y=self.grid.latp

    # rotate uv:
    ang=rt.rho2uvp(self.grid.angle,'p')
    u,v=calc.rot2d(u,v,-ang)


    if plot:
      p=self.grid.plot(bathy=None)
      xm, ym = p(x,y)
      mm=0*m
      mm[::3,::3]=1
      mm=mm*m==1
      s=np.sqrt(u**2+v**2)
      q=pylab.quiver(xm[mm],ym[mm],u[mm],v[mm],s[mm])
      pylab.colorbar(shrink=.7)
      pylab.axis([p.xmin, p.xmax, p.ymin, p.ymax])
      qk = pylab.quiverkey(q, .05, .01, 0.2, '0.2 m/s',labelpos='E',
                                              coordinates='axes')
      if savename:
        pylab.savefig(savename,dpi=pylab.gcf().dpi)
        pylab.close(pylab.gcf())

    return x,y,z,u,v


  def sliceiso(self,vname,iso,time,plot=False,**opts):
    '''
    Depths where variable, increasing/dec with depth, has some value.

    Output is masked where surface is higher/lower than value or where all
    water column is lower/higher than value.
    '''

    v=self.use(vname,SEARCHtime=time)
    z=self.s_levels(time=time,ruvpw=self.var_at(vname))
    v=rt.depthof(v,z,iso)

    if plot:
      p=self.grid.plot(bathy=None)
      xm, ym = p(x,y)
      pylab.pcolormesh(xm,ym,v)
      pylab.colorbar()

    return x,y,v


  def time_series(self,varname,x,y,times=False,depth=False,**opts):
    t,z,v=[[]]*3

    RetVarOnly=False
    savename=False
    retMsg=False
    nearest=True  # UNIQUE case Currently
    plot=False

    if 'retvar'   in opts.keys(): RetVarOnly = opts['retvar']
    if 'savename' in opts.keys(): savename   = opts['savename']
    if 'msg'      in opts.keys(): retMsg     = opts['msg']
    if 'nearest'  in opts.keys(): nearest    = opts['msg']
    if 'plot'     in opts.keys(): plot       = opts['msg']


    if varname not in netcdf.varnames(self.name):
      msg=':: variable %s not found' % varname
      if retMsg: return t,z,v,msg
      else:
        print msg
        return t,z,v

    isW=varname=='w'
    hasZ=self.hasz(varname)

    if not depth is False and hasZ :
      if isW and depth >= self.S_W:
        msg='k = %d exceeds S_W dimension (%d)' % (depth,self.S_W)
        if retMsg: return t,z,v,msg
        else:
          print msg
          return t,z,v

      elif depth >= self.S_RHO:
        msg='k = %d exceeds S_RHO dimension (%d)' % (depth,self.S_RHO)
        if retMsg: return t,z,v,msg
        else:
          print msg
          return t,z,v

    if times is False:
      times=0,len(self.tdays)

    if self.hast(varname) and times[-1]>self.TIME:
      msg='t = %d exceeds TIME dimension (%d)' % (times[-1],self.TIME)
      if retMsg: return t,z,v,msg
      else:
        print msg
        return t,z,v

    lon,lat,hr,mr=self.grid.vars(ruvp=self.var_at(varname))
    dist=(lon-x)**2+(lat-y)**2
    i,j=np.where(dist==dist.min())
    i,j=i[0],j[0]
    v=self.use(varname,xiSEARCH=j,etaSEARCH=i,SEARCHtime=range(times[0],times[-1])).T

    # time:
    t=self.tdays[range(times[0],times[-1])]
    t=t+0.*v

    # depth:
    if hasZ:
      h=self.grid.h[i,j]
      zeta=self.use('zeta',xiSEARCH=j,etaSEARCH=i,SEARCHtime=range(times[0],times[-1]))
      h=h+0*zeta
      z=rt.s_levels(h,zeta,self.s_params,rw=varname)
      z=np.squeeze(z)

      if not depth is False:
        if depth>=0:
          t=t[0,:]
          z=z[depth,:]
          v=v[depth,:]
        else:
          t0,z0=t,z
          t=t[0,:]
          z=depth+0.*t
          v=calc.griddata(t0,z0,v,t,z,extrap=False,norm_xy=True)

    else: # not hasZ
      z=0.*t


    if plot:
      pass # TODO


    if RetVarOnly:
      return v
    else:
      if retMsg: return t,z,v,''
      else: return t,z,v



  def extrap_at_mask(self,vname,time,quiet=1,**opts):
    if not quiet: print 'Extrap at mask:\n  loading'
    v=self.use(vname,SEARCHtime=time)
    x,y,h,m=self.grid.vars(vname)

    if self.hasz(vname): # 3d
      N=v.shape[0]
      for n in range(N):
        if not quiet: print '  extrap level %d of %d' % (n,N)
        v[n,...]=calc.mask_extrap(x,y,np.ma.masked_where(m==0,v[n,...]))
    else: # 2d
        v=calc.mask_extrap(x,y,np.ma.masked_where(m==0,v))

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
      self.datetime=netcdftime.num2date(self.time,self.var_as['time']['units'])
    except: self.datetime=False

    # time is usually seconds, but may be days!!, so:
    if self.var_as['time']['units'].strip().startswith('days'): self.time=self.time*86400.

    if len(self.time)>1:
      self.dt=self.time[1]-self.time[0]
    else: self.dt=0

    self.tdays=self.time/86400.


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


