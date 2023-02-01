import numpy as np
import os

from okean import calc, netcdf, vis, ticks, cookbook as cb
from okean.roms import roms_tools as rt
from okean.roms.derived import Derived

try:
  from mpl_toolkits.basemap import Basemap
except:
  Basemap=False

__all__='His','Grid','MHis','MGrid'

class propVar():
  def __init__(self,vname): self.vname=vname
  def get(self,gob): return gob._get(self.vname)


class Common():
  def _get(self,vname):
    if not hasattr(self,'_store'): self._store={}
    if not vname in self._store:
      self._store[vname]=self.use(vname)

    return self._store[vname]

  def load_vars(self,vars):
    var=netcdf.var(self.nc)
    for i in vars:
      if not self.quiet: print(':: loading var %s'%i)

      if isinstance(i,dict):
        iname, iopts = list(i.items())[0]
      elif cb.isstr(i):
        iname=i
        if iname.startswith('@'): iopts=iname[1:]
        else: iopts=iname

      if not isinstance(iopts,tuple):
        iopts=(iopts,)

      found=0
      if not hasattr(self,'var_as'): self.var_as={}
      for j in iopts:
        if j in var:
          if iname.startswith('@'): # variable will be loaded only when needed!
            iname=iname[1:]

            # can't use self._get inside lambda function... j is always the last element!
            # so, instead of:
            # setattr(self.__class__, iname, property(fget=lambda self: self._get(j)))
            # a propVar variable is used, storing the riht j value:
            ob=propVar(j)
            setattr(self.__class__,iname,property(fget=ob.get))##lambda self: self.getx(vv)))
          else:
            setattr(self,iname,var[j][:])

          found=1
          if not self.quiet: print('  - found %s as %s'%(iname,j))

          # store original name and other info, like units
          try: units=var[j].atts['units'].value
          except: units=''
          self.var_as[iname]={'name':j,'units':units}
          break
      if not found and not self.quiet: print(':: %s not found in %s'%(str(iopts),self.name))

  def load_dims(self):
    dms=netcdf.fdim(self.nc)
    for k in dms:
      if k in ('ocean_time','time','scrum_time'):
        setattr(self,'TIME',dms[k])
      else:
        setattr(self,k.upper(),dms[k])

  def load_atts(self):
    atts=netcdf.fatt(self.nc)
    self.atts={}
    self.atts.update(atts)

  def load_grid(self,grd=False):
    ats=netcdf.fatt(self.nc)
    if not grd:
      if 'grd_file' in ats: grd=ats['grd_file'].value
      if grd and not grd.startswith('http') and not os.path.isfile(grd): grd=False

    if grd: self.grid=Grid(grd,self.quiet)
    else:
      try:
        self.grid=Grid(self.name)
      except:
        self.grid=False

  def show(self): netcdf.show(self.nc._nc)

  def showvar(self,vname): netcdf.showvar(self.nc,vname)

  def use(self,varname,**kargs):
    return netcdf.use(self.nc,varname,**kargs)

  def vaxes(self,v):
    '''
    returns variables axes as string:
    ex: 'vaxes('temp') is 'tzyx'
    indicating time (t), vertical (z) and horizontal (y,x) dimensions
    '''

    dnames=netcdf.vdim(self.nc,v)
    out=''
    for d in dnames:
      if d.find('time')>=0: out+='t'
      elif d in ('s_rho','s_w'): out+='z'
      if d.startswith('eta_'): out+='y'
      if d.startswith('xi_'): out+='x'

    return out

  def vloc(self,v):
##  def var_at(self,v):##,_3d=True):
    '''
    returns location of the variable in the 2d/3d grid
    Possible values are u,v,r(rho), p(psi)
    If _3d, also checks if at w points
    '''
    dims=netcdf.vdim(self.nc,v)

    hLoc='r' # r,u,v
    vLoc='r' # r,w

    if   's_w' in dims: vLoc='w'
    if   'xi_u'   in dims or   'eta_u'   in dims: hLoc='u' # or instead of and cos of agrif ...
    elif 'xi_v'   in dims or   'eta_v'   in dims: hLoc='v'
    elif 'xi_psi' in dims and  'eta_psi' in dims: hLoc='p'
#    return hLoc,vLoc
    return hLoc+vLoc

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
    self.load_dims()
    self.load()
    # load_uvp was used to create uvp variables based on rho variables. The
    # objective was to load few variables in case of remote opendap url. In
    # the current implementation the variables are loaded only when needed,
    # so, there is no big gain not to load uvp variables from source.
    # User may still call self.load_uvp() after object creation to ensure, ie,
    # to force that less data is downloaded from remote grid
    if 0: self.load_uvp()
    self.load_atts()

    # about proj:
    self.load_proj_info()


  def load_proj_info(self):
    # 1) load from grid and parse attrs needed for proj

    self.proj_info={}
    for i in ['proj4_string','wkt_string','proj_dict']:
      if i in netcdf.fatt(self.nc):
        self.proj_info[i]=netcdf.fatt(self.nc,i)
      else:
        self.proj_info[i]=''

    # parse dict string:
    o={}
    if self.proj_info['proj_dict']:
      for i in self.proj_info['proj_dict'].split():
        name,val=i[1:].split('=')
        try:
          val=float(val)
        except: pass
        o[name]=val

    self.proj_info['proj_dict']=o

    # 2) create a set of reasonable projection options from the data;
    #  user can always change them after. Let us consider only the lcc
    # and merc and tmerc projections, and use lcc by default if none
    # is provided in proj_info. Also let's make sure the projection is
    # well centred (may not be at grid creation)

    if o: prj=o['proj']
    else: prj='lcc'

    pnew={}
    if self.spherical:
      xlim=self.lon.min(),self.lon.max()
      ylim=self.lat.min(),self.lat.max()
      dlon=xlim[1]-xlim[0]
      dlat=ylim[1]-ylim[0]
      lon_c=0.5*(xlim[0]+xlim[1])
      lat_c=0.5*(ylim[0]+ylim[1])

      if prj=='lcc':
        lat_1=lat_c-dlat/4.
        lat_2=lat_c+dlat/4.
        pnew=dict(proj=prj,lon_0=lon_c,lat_0=lat_c,lat_1=lat_1,lat_2=lat_2)
      elif prj=='merc':
        pnew=dict(proj=prj,lon_0=lon_c,lat_ts=lat_c)
      elif prj=='tmerc':
        pnew=dict(proj=prj,lon_0=lon_c,lat_0=lat_c)

      self.proj_info['proj_dict_nice']=pnew

      # 3) keep data extent, with some margin (needed for plots)
      marginX=dlon/10
      marginY=dlat/10
      self.proj_info['extent']=xlim[0]-marginX,xlim[1]+marginX, ylim[0]-marginY,ylim[1]+marginY
      self.proj_info['extent_no_margin']=xlim[0],xlim[1],ylim[0],ylim[1]


  def get_projection(self,prjtype='nice',cartopy=False,**ckargs):
    '''
    prjtype, proj type, if 'nice', a set of reasonable options are used,
      else if 'original', the original basemap_opts0 is used
    if cartopy, returns cartopy object instead of pyproj
    ckargs=cartopy extra arguments (e.g. cutoff for lcc projection)
    '''

    if prjtype=='nice': opts=self.proj_info['proj_dict_nice']
    elif prjtype=='original':
      # use WKT if present, otherwise use proj_dict
      if self.proj_info['wkt_string']: opts=self.proj_info['wkt_string']
      else: opts=self.proj_info['proj_dict']
    else: return 'unknown prjtype'

    if cartopy:
      import cartopy.crs as ccrs
      if opts['proj']=='lcc':
        #if opts['lat_0']<0: cutoff=30
        #else: cutoff=-30
        ylim=self.lat.min(),self.lat.max()
        yc=(ylim[0]+ylim[1])/2
        if yc>0: cutoff=ylim[0]-1
        else: cutoff=ylim[1]+1

        cutoff=ckargs.pop('cutoff',cutoff)
        p=ccrs.LambertConformal(central_longitude=opts['lon_0'],
                                central_latitude=opts['lat_0'],
                                #secant_latitudes=(opts['lat_1'],opts['lat_2']),
                                # secant latitudes deprecated, use standard paralels
                                # (https://scitools.org.uk/cartopy/docs/v0.17/crs/projections.html)
                                standard_parallels=(opts['lat_1'],opts['lat_2']),
                                cutoff=cutoff,**ckargs)
      elif opts['proj']=='merc':
        p=ccrs.Mercator(central_longitude=opts['lon_0'],
                        latitude_true_scale=opts['lat_ts'],**ckargs)
      elif opts['proj']=='tmerc':
        p=ccrs.TransverseMercator(central_longitude=opts['lon_0'],
                                  central_latitude=opts['lat_0'],**ckargs)
      else: return 'not implemented yet'
      return p
    else:
      import pyproj
      return pyproj.Proj(opts)


#  def load_proj(self):
#    # load and parse attrs needed for basemap projection
#
#    self.proj_info={}
#    for i in ['proj4string','basemap_opts']:
#      if i in netcdf.fatt(self.nc):
#        self.proj_info[i]=netcdf.fatt(self.nc,i)
#      else:
#        self.proj_info[i]=''
#
#    o={}
#    if self.proj_info['basemap_opts']: # parse it
#      for i in self.proj_info['basemap_opts'].split():
#        name,val=i[1:].split('=')
#        try:
#          val=float(val)
#        except: pass
#        o[name]=val
#
#    self.proj_info['basemap_opts0']=o
#
#    # Now let ys make sure the projection is well centered (may not be
#    # at grid creation):
#    # - build a set of reasonable options;  user can always change them
#    # after. Let us consider only the lcc and merc projections, and use
#    # lcc by default if none is provided at grid atts
#    if o: prj=o['projection']
#    else: prj='lcc'
#
#    if self.spherical:
#      xlim=self.lon.min(),self.lon.max()
#      ylim=self.lat.min(),self.lat.max()
#      dlon=xlim[1]-xlim[0]
#      dlat=ylim[1]-ylim[0]
#
#      if prj=='lcc':
#        lon_0=0.5*(xlim[0]+xlim[1])
#        lat_0=0.5*(ylim[0]+ylim[1])
#        lat_1=lat_0-dlat/4.
#        lat_2=lat_0+dlat/4.
#
#        W=110*np.cos(lat_0*np.pi/180)*(xlim[1]-xlim[0])*1e3
#        H=110*(ylim[1]-ylim[0])*1e3
#        W*=1.2
#        H*=1.2
#
#        o=dict(resolution='i',projection=prj,
#               lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0,
#               width=W,height=H)
#
#      elif prj=='merc':
#        dlon_=dlon/10
#        dlat_=dlat/10
#        o=dict(resolution='c',projection=prj,
#               llcrnrlon=xlim[0]-dlon_, llcrnrlat=ylim[0]-dlat_,
#               urcrnrlon=xlim[1]+dlon_, urcrnrlat=ylim[1]+dlat_,
#               lon_0=0.5*(xlim[0]+xlim[1]),
#               lat_0=0.5*(ylim[0]+ylim[1]))
#
#      self.proj_info['basemap_opts']=o
#      self.proj_info['extent']=xlim[0]-dlon/10,xlim[1]+dlon/10, ylim[0]-dlat/10,ylim[1]+dlat/10
#
#
#  def get_projection(self,prjtype='nice'):
#    '''
#    prjtype, proj type, if 'nice', a set of reasonable options are used,
#      else if 'original', the original basemap_opts0 is used
#    '''
#
#    if prjtype=='nice': opts=self.proj_info['basemap_opts']
#    elif prjtype=='original': opts=self.proj_info['basemap_opts0']
#    else: return 'unknown prjtype'
#
#    if opts: return Basemap(**opts)


  def load(self):
    # if grid file is remote (opendap) we could load less variables. But since
    # @varname will be loaded only when used, there is no need to do it because
    # by default almost no variables are loaded
    vars={'@lon':('lon_rho','x_rho')}, {'@lonu':('lon_u','x_u')},\
         {'@lonv':('lon_v','x_v')},    {'@lonp':('lon_psi','x_psi')},\
         {'@lat':('lat_rho','y_rho')}, {'@latu':('lat_u','y_u')},\
         {'@latv':('lat_v','y_v')},    {'@latp':('lat_psi','y_psi')},\
         {'mask':'mask_rho'},{'@masku':'mask_u'},{'@maskv':'mask_v'},{'@maskp':'mask_psi'},\
         '@h','angle', '@pn','@pm','spherical'

    self.load_vars(vars)

    # some vars may not be present... set defaults:
    if not hasattr(self,'mask'): self.mask=np.ones((self.ETA_RHO,self.XI_RHO))
    if not hasattr(self,'angle'): self.angle=np.zeros((self.ETA_RHO,self.XI_RHO))

    # about spherical variable:
    # may be 1, T, an array with 'T','','',...
    # let's also assume that empty or None sph means True!! (may not be a very smart idea)
    sph=self.spherical
    trueVals=[1,'T','',None]

    # check if numeric:
    if np.issubdtype(sph.dtype,np.number):
      #if sph.size: sph=sph[0] # may have size and no shape !
      if np.sum(sph.shape): sph=sph[0]
    else:

      if sph.shape==(): # another netcdf4 version, another issue!!
        try:
          sph=sph[()].decode()
        except:
          sph=sph[()]
      else:
        try:
          try:
            sph=''.join(sph).strip()
          except TypeError:
            if sph.shape==(): # size is 1 though !
              sph=sph[()].decode()
            else:
              sph=b''.join(sph).strip().decode() # bytes type needed for python3

          if len(sph): sph=sph[0]
        except: pass

    self.spherical=sph in trueVals

  def load_grid(): self.load()

  def load_uvp(self):
    '''
    Set uvp lon, lat and mask based on variables at rho, instead of
    loading them from file
    '''
    [setattr(self,'lon'+i,rt.rho2uvp(self.lon,i)) for i in 'uvp']
    [setattr(self,'lat'+i,rt.rho2uvp(self.lat,i)) for i in 'uvp']
    self.masku,self.maskv,self.maskp=rt.uvp_mask(self.mask)

    if not hasattr(self,'XI_U'):  self.XI_U  = self.XI_RHO-1
    if not hasattr(self,'XI_V'):  self.XI_V  = self.XI_RHO
    if not hasattr(self,'ETA_U'): self.ETA_U = self.ETA_RHO
    if not hasattr(self,'ETA_V'): self.ETA_V = self.ETA_RHO-1


  def border_multi_corner(self):
    mc=self.use('maskc')
    if np.unique(mc).size==1: return self.border() # it is not multi-corner!

    l,n=mc.shape

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
    '''
    return lon, lat, h and mask at r, u, v or p
    i,j slices can also be used
    '''

    ruvp=ruvp[0]
    if ruvp=='u':
      x = self.lonu
      y = self.latu
      h = rt.rho2uvp(self.h,'u')
      m = self.masku
    elif ruvp=='v':
      x = self.lonv
      y = self.latv
      h = rt.rho2uvp(self.h,'v')
      m = self.maskv
    elif ruvp=='p':
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

    #x,y,h,m=[i.data for i in [x,y,h,m] if np.ma.isMA(i)]
    x,y,h,m=[k.data if  np.ma.isMA(k) else k for k in [x,y,h,m]]

    return x[j,i],y[j,i],h[j,i],m[j,i]


  def s_levels(self,sparams,zeta=0,h=False,loc='rr',i=False,j=False,k=False):
    hLoc,vLoc=loc
    if h is False: h=self.h

    try:
      zeta.shape==h.shape
    except:
      zeta=np.tile(zeta,h.shape).astype(h.dtype)

    if h.ndim==2:
      h=rt.rho2uvp(h,hLoc)
      zeta=rt.rho2uvp(zeta,hLoc)

    z=rt.s_levels(h,zeta,sparams,rw=vLoc)

    if k is False: k=slice(None)
    if j is False: j=slice(None)
    if i is False: i=slice(None)

########################    return np.squeeze(z[k,j,i])
    return z[k,j,i]


  def ingrid(self,x,y,**kargs):#,retInds=False):
    '''in polygon of region border'''
    xb,yb=self.border(**kargs)
    return calc.inpolygon(x,y,xb,yb)

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

      try:
        from mayavi import mlab
      except:
        print('error: mayavi is required')
        return

      import sys
      pver=eval(sys.version[:3]) 
      qtapi=os.environ.get('QT_API','')

      if pver<3 and qtapi!='pyqt':
        print('warning: environmen variable QT_API may need to be set as pyqt '\
              'in python <3 so that mayavi can be used after pylab')

        a=raw_input('Wanna continue ([n],y) ?')
        if a!='y': return


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

  def plot_data(self,**kargs):
    prjtype=kargs.pop('prjtype','nice')
    from matplotlib.cm import gist_earth_r
    h=np.ma.masked_where(self.mask==0,self.h)
    a=vis.Data(x=self.lon,y=self.lat,v=h)
    ht=ticks.loose_label_n(self.h.min(),self.h.max(),7)
    a.set_param(field__plot='contourf',field__cvals=ht,field__cmap=gist_earth_r)
    a.set_param(**kargs)

    # also show domain boundary:
    xb,yb=self.border()
    b=vis.Data(x=xb,v=yb)
    b.set_param(d1_line__options=dict(lw=0.5,color='k',ls='-',zorder=1e3))
    a.extra=[b]

##    # about projection:
##    if prjtype=='original' and self.proj_info['basemap_opts0']: d=self.proj_info['basemap_opts0']
##    elif prjtype=='nice': d=self.proj_info['basemap_opts']
##    a.set_projection(d,extent=self.proj_info['extent'])

    p=self.get_projection(prjtype,cartopy=1)
    a.set_projection(p,extent=self.proj_info['extent'])

    return a


  def plot(self,**kargs):
    a=self.plot_data(**kargs)
    labels=kargs.pop('labels',0)
    kargs.pop('prjtype',0)
    a.plot(labels=labels,**kargs)
    return a


class MGrid():
  def __init__(self,name):
    if cb.isstr(name):
      import glob
      f=os.path.splitext(name)[0][:-1]+'*'+os.path.splitext(name)[1]
      files=glob.glob(f)
      files.sort()
    else:
      files=name

    if cb.isstr(files[0]): # list of filenames
      self._grids=[Grid(i) for i in files]
    else: self._grids=name # list of Grid, used by MHis only

  def __getitem__(self,i):
    return self._grids[i]

  def __len__(self): return len(self._grids)

  def __repr__(self):
    return '<MGrid at 0x%x>'%id(self)+\
    '\n'+''.join([' -- '+os.path.basename(i.name)+'\n' for i in self._grids])

  def _set_ingrid(self,ruvp='r'):
    for i in range(len(self)):
      aname='ingrd_'+ruvp
      setattr(self[i],aname,[])
      ob=getattr(self[i],aname)
      if i<len(self)-1:
        x,y=self[i].vars(ruvp)[:2]
        for j in range(i+1,len(self)):
          ob+=[self[j].ingrid(x,y)]

  def border(self,**kargs):
    xb=[]
    yb=[]
    n=0
    for i in self._grids:
      x,y=i.border(**kargs)
      xb+=[x]
      yb+=[y]
      n=max(n,x.size)

    X=np.full((n,len(self._grids)),np.nan)
    Y=np.full((n,len(self._grids)),np.nan)

    for i in range(len(self._grids)):
      X[:xb[i].size,i]=xb[i]
      Y[:yb[i].size,i]=yb[i]

    return X,Y

  def plot_data(self,**kargs):
    prjtype=kargs.pop('prjtype','nice')

    from matplotlib.cm import gist_earth_r

    if 0:
      j=np.ceil(np.median(range(len(self)))).astype('i')
      ht=ticks.loose_label_n(self[j].h.min(),self[j].h.max(),7)
    else:
      ht=ticks.loose_label_n(self[0].h.min(),self[0].h.max(),7)

    if len(self)>1:
      ht2=ticks.loose_label_n(ht[0],ht[1]/2,5)[1:]
      ht=np.concatenate(([ht[0]],ht2,ht[1:]))
      ht=np.unique(np.sort(ht))

    ht=kargs.get('field__cvals',ht)
    ht=kargs.get('cvals',ht)

    # contour colors:
    from matplotlib import rcParams
    colors=rcParams['axes.prop_cycle']()

    for c,i in enumerate(self):
      h=np.ma.masked_where(i.mask==0,i.h)
      b=vis.Data(x=i.lon,y=i.lat,v=h)
      b.set_param(field__plot='contourf',field__cvals=ht,field__cmap=gist_earth_r)
      b.set_param(**kargs)
      if c==0:
        a=b
      else:
        a.extra+=[b]

      c=vis.Data(x=i.lon,y=i.lat,v=h)
      c.set_param(field__plot='contour',field__cvals=ht,field__cmap=next(colors)['color'],
                  field__linewidths=0.5,plot__zorder=2) # show above everything else
      c.label='bathy'
      a.extra+=[c]

    # show all borders:
    for i in self:
      xb,yb=i.border()
      c=vis.Data(x=xb,v=yb)
      c.set_param(d1_line__options=dict(lw=0.5,color='k',ls='-'))
      a.extra+=[c]


    # about projection:
##    if prjtype=='original' and self[0].proj_info['basemap_opts0']: d=self[0].proj_info['basemap_opts0']
##    elif prjtype=='nice': d=self[0].proj_info['basemap_opts']
##    a.set_projection(d,extent=self[0].proj_info['extent'])

    p=self[0].get_projection(prjtype,cartopy=1)
    a.set_projection(p,extent=self[0].proj_info['extent'])

    return vis.MData([a])


  def plot(self,**kargs):
    a=self.plot_data(**kargs)
    kargs.pop('prjtype',0)
    a.plot(**kargs)


class His(Common,Derived):
  def __init__(self,his,grd=False,quiet=True):
    self.name=his
    self.isremote=False
    if cb.isstr(his):
      if his.startswith('http'): self.isremote=True
      else: self.name=os.path.realpath(his)

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

    if np.ma.is_masked(zeta):
      # h=np.ma.masked_where(zeta.mask,h)
      # comment above in order to keep right depth inside the domain
      # on masked regions
      pass

    else:# probably a clm/ini file. Mask maskc point (pygridgen complex grid):
      if 'maskc' in netcdf.varnames(self.grid.nc):
        mc=self.grid.use('maskc')
        if inds: mc=mc[j0:j1,i0:i1]
        h=np.ma.masked_where(mc==0,h)
        zeta=np.ma.masked_where(mc==0,zeta)

    h    = calc.griddata(xr,yr,h,x,y,extrap=False)
    zeta = calc.griddata(xr,yr,zeta,x,y,extrap=False)

    # better interpolate the gaps, otherwise s_levels will use default
    # values for zeta and h:
    X=np.arange(h.size)
    if np.ma.is_masked(h): h=np.interp(X,X[~h.mask],h[~h.mask])
    if np.ma.is_masked(zeta): zeta=np.interp(X,X[~zeta.mask],zeta[~zeta.mask])
    

###    return h,zeta,self.s_params,rw
    return rt.s_levels(h,zeta,self.s_params,rw=rw)
###    z=np.squeeze(z)
###    return np.ma.masked_where(np.isnan(z),z)


  def s_levels(self,time,loc='rr',i=False,j=False,k=False,extrapZeta=False):
    try:
      hLoc,vLoc=loc
    except:
      hLoc,vLoc=loc,'r'

    h=self.grid.h
    zeta=self.use('zeta',SEARCHtime=time)

    if extrapZeta:
      #if not calc.ismarray(zeta): zeta=np.ma.masked_where(self.grid.mask==0,zeta)
      if not np.ma.is_masked(zeta): zeta=np.ma.masked_where(self.grid.mask==0,zeta)
      zeta=calc.mask_extrap(self.grid.lon,self.grid.lat,zeta)

    h=rt.rho2uvp(h,hLoc)
    zeta=rt.rho2uvp(zeta,hLoc)
    z=rt.s_levels(h,zeta,self.s_params,rw=vLoc)

    if k is False: k=slice(None)
    if j is False: j=slice(None)
    if i is False: i=slice(None)

    return z[k,j,i]


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
        if 't' in self.vaxes(varname) and ind>=self.TIME:
          msg='t = %d exceeds TIME dimension (%d)' % (ind,self.TIME)

      # check dim k:
      elif dim=='k':
        if 'z' in self.vaxes(varname):
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
    return ','.join(sorted(res.split(',')))

  def slicei(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicei')).split(',')

    out=vis.Data()
    out.label='slicei'
    out.msg=self.check_slice(varname,t=time,i=ind)
    if out.msg: return out

    v=self.use(varname,SEARCHtime=time,xi_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v):
      m=self.grid.vars(ruvp=self.vloc(varname)[0],i=ind)[-1]
      if v.ndim==2: m=np.tile(m,(v.shape[0],1))
      v=np.ma.masked_where(m==0,v)

    out.v=v
    out.info['v']['name']=varname
    out.info['v']['slice']='i=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass

    # coords:
    if 'z' in coords and v.ndim==2:
      out.z=self.s_levels(time=time,loc=self.vloc(varname),i=ind)
      out.info['z']=dict(name='Depth',units='m')

    if any([i in coords for i in 'xyd']):
      x,y,h,m=self.grid.vars(ruvp=self.vloc(varname)[0],i=ind)

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

    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[vis.Data()]
      if 'd' in coords: out.extra[0].x=out.d[0]
      if 'y' in coords: out.extra[0].x=out.y[0]
      if 'x' in coords: out.extra[0].y=out.x[0]
      out.extra[0].v=-h
      out.extra[0].config['d1.plot']='fill_between'
      out.extra[0].config['d1.y0']=-h.max()-(h.max()-h.min())/20.
      out.extra[0].label='bottom'

    out.coordsReq=','.join(sorted(coords))
    return out


  def slicej(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicej')).split(',')

    out=vis.Data()
    out.label='slicej'
    out.msg=self.check_slice(varname,t=time,j=ind)
    if out.msg: return out

    v=self.use(varname,SEARCHtime=time,eta_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v):
###      m=self.grid.vars(ruvp=self.var_at(varname)[0],j=ind)[-1]
      m=self.grid.vars(ruvp=self.vloc(varname)[0],j=ind)[-1]
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
###      out.z=self.s_levels(time=time,loc=self.var_at(varname),j=ind)
      out.z=self.s_levels(time=time,loc=self.vloc(varname),j=ind)
      out.info['z']=dict(name='Depth',units='m')

    if any([i in coords for i in 'xyd']):
###      x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0],j=ind)
      x,y,h,m=self.grid.vars(ruvp=self.vloc(varname)[0],j=ind)

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

########    if 't' in coords and self.hast(varname): out.t=self.time[time]
    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if v.ndim==2:
      out.extra=[vis.Data()]
      if 'd' in coords: out.extra[0].x=out.d[0]
      if 'x' in coords: out.extra[0].y=out.x[0]
      if 'y' in coords: out.extra[0].x=out.y[0]
      out.extra[0].v=-h
      out.extra[0].config['d1.plot']='fill_between'
      out.extra[0].config['d1.y0']=-h.max()-(h.max()-h.min())/20.
      out.extra[0].label='bottom'

    out.coordsReq=','.join(sorted(coords))
    return out


  def slicek(self,varname,ind,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicek')).split(',')

    out=vis.Data()
    out.label='slicek'
    out.msg=self.check_slice(varname,t=time,k=ind)
    if out.msg: return out

    v=self.use(varname,SEARCHtime=time,s_SEARCH=ind)

    # add mask if not masked:
    if not np.ma.isMA(v): 
###      m=self.grid.vars(ruvp=self.var_at(varname)[0])[-1]
      m=self.grid.vars(ruvp=self.vloc(varname)[0])[-1]
      v=np.ma.masked_where(m==0,v)

    out.v=v
    out.info['v']['name']=varname
#####    if self.hasz(varname): out.info['v']['slice']='k=%d'%ind
    if 'z' in self.vaxes(varname): out.info['v']['slice']='k=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
    if 'z' in coords and 'z' in self.vaxes(varname):
      out.z=self.s_levels(time,k=ind,loc=self.vloc(varname))
      out.info['z']=dict(name='Depth',units='m')


    if any([i in coords for i in 'xy']):
      x,y,h,m=self.grid.vars(ruvp=self.vloc(varname)[0])

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

    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if out.v.ndim==2: # always?
      out.extra=[vis.Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.cmap']='k'
      out.extra[0].label='bathy'


    out.coordsReq=','.join(sorted(coords))

    # about projection:
##    out.set_projection(self.grid.proj_info['basemap_opts'],extent=self.grid.proj_info['extent'])
    p=self.grid.get_projection(cartopy=1)
    out.set_projection(p,extent=self.grid.proj_info['extent'])

    return out


  def slicez(self,varname,ind,time=0,**opts):
    surf_mask=opts.get('surf_mask',True)
    spline=opts.get('spline',True)
    coords=opts.get('coords',self._default_coords('slicez')).split(',')

    out=vis.Data()
    out.label='slicez'
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out

    if not 'z' in self.vaxes(varname):
      return self.slicek(varname,ind,time,**opts)

    v=self.use(varname,SEARCHtime=time)
###    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0])
    x,y,h,m=self.grid.vars(ruvp=self.vloc(varname)[0])
    zeta=self.use('zeta',SEARCHtime=time)
    zeta=rt.rho2uvp(zeta,self.vloc(varname)[0])

    out.v=rt.slicez(v,m,h,zeta,self.s_params,ind,surf_mask,spline)

    out.info['v']['name']=varname
    if calc.isarray(ind):
      out.info['v']['slice']='z= array %.2f to %.2f'%(ind.min(),ind.max())
    else:
      out.info['v']['slice']='z=%d'%ind
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass


    # coords:
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

    if 'z' in coords:
      out.z=ind+np.zeros(out.v.shape)
      out.info['z']=dict(name='Depth',units='m')

    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if out.v.ndim==2:
      out.extra=[vis.Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.cmap']='k'
      out.extra[0].label='bathy'


    out.coordsReq=','.join(sorted(coords))

    # about projection:
##    out.set_projection(self.grid.proj_info['basemap_opts'],extent=self.grid.proj_info['extent'])
    p=self.grid.get_projection(cartopy=1)
    out.set_projection(p,extent=self.grid.proj_info['extent'])

    return out


  def slicell(self,varname,X,Y,time=0,**opts):
    coords=opts.get('coords',self._default_coords('slicell')).split(',')

    data      = opts.get('data',False)
    extrap    = opts.get('extrap',False)
    maskLimit = opts.get('lmask',0.5) # points where interpolated mask are above
                                      # this value are considered as mask!
                                      # Most strict value is 0

    out=vis.Data()
    out.label='slicell'
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out#None,aux

    X=np.asarray(X)
    Y=np.asarray(Y)
    if X.ndim>1: X=np.squeeze(X)
    if Y.ndim>1: Y=np.squeeze(X)

###    x,y,h,m=self.grid.vars(ruvp=self.var_at(varname)[0])
    x,y,h,m=self.grid.vars(ruvp=self.vloc(varname)[0])
    if True: # extrat only portion of data needed:
      i0,i1,j0,j1=calc.ij_limits(x, y, (X.min(),X.max()),(Y.min(),Y.max()), margin=1)
      xi='%d:%d'%(i0,i1)
      eta='%d:%d'%(j0,j1)

      if data is False: V=self.use(varname,SEARCHtime=time,xi_SEARCH=xi,eta_SEARCH=eta)
      #else: v=data[...,j0:j1,i0:i1]
      else: V=data[...,j0:j1,i0:i1]

      x=x[j0:j1,i0:i1]
      y=y[j0:j1,i0:i1]
      h=h[j0:j1,i0:i1]
      m=m[j0:j1,i0:i1]

    else:
      if data is False: V=self.use(varname,SEARCHtime=time)
      #else: v=data
      else: V=data

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
###      out.z=self.path_s_levels(time,X,Y,rw=self.var_at(varname)[1],inds=inds)
###
#######      out.z,zw=self.path_s_levels(time,X,Y,rw=False,inds=inds)
#######      if self.vloc(varname)[1]=='w': out.z=zw
      out.z=self.path_s_levels(time,X,Y,rw=self.vloc(varname)[1],inds=inds)
      out.info['z']=dict(name='Depth',units='m')

    if 'd' in coords:
      d=calc.distance(X,Y)
      if d[-1]-d[0]>1e4:
        d=d/1000.
        dunits='km'
      else: dunits='m'

      if v.ndim==2: d=np.tile(d,(v.shape[0],1))
      out.d=d
      out.info['d']=dict(name='Distance',units=dunits)

    if 'x' in coords:
      if v.ndim==2: X=np.tile(X,(v.shape[0],1))
      out.x=X
      out.info['x']=dict(name='Longitude',units=r'$\^o$E')

    if 'y' in coords:
      if v.ndim==2: Y=np.tile(Y,(v.shape[0],1))
      out.y=Y
      out.info['y']=dict(name='Latitude',units=r'$\^o$N')

#######    if 't' in coords and self.hast(varname): out.t=self.time[time]
    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if v.ndim==2: ################3 and not out.z is None: # zeta and bottom already calculated
      out.extra=[vis.Data()]
      if 'd' in coords: out.extra[0].x=out.d[0]
      if 'x' in coords: out.extra[0].y=out.x[0]
      if 'y' in coords: out.extra[0].x=out.y[0]
####      #h=-zw[0]
      h    = calc.griddata(x,y,h,X,Y,extrap=False)
      out.extra[0].v=-h # bottom
      out.extra[0].config['d1.plot']='fill_between'
      out.extra[0].config['d1.y0']=-h.max()-(h.max()-h.min())/20.
      out.extra[0].label='bottom'



    out.coordsReq=','.join(sorted(coords))
    return out


  def sliceuv(self,ind,time=0,**opts):#plot=False,**opts):
    coords=opts.get('coords',self._default_coords('sliceuv')).split(',')

    if ind=='bar':
      slc,uname,vname,ind = self.slicek, 'ubar','vbar', 9999
    elif ind in ('s','surf','surface'):
      slc,uname,vname,ind = self.slicek,  'u','v', -1
    elif ind>=0:
      slc,uname,vname,ind = self.slicek,  'u','v', ind
    elif ind <0:
      slc,uname,vname,ind = self.slicez, 'u','v', ind

    outu=slc(uname,ind,time,**opts)
    outv=slc(vname,ind,time,**opts)

    if   outu.msg: return outu
    elif outv.msg: return outv

    # at psi:
    u,v=outu.v,outv.v
    u=( u[1:,:]+ u[:-1,:])/2.
    v=( v[:,1:]+ v[:,:-1])/2.

    # rotate uv:
    if 0:
      ang=rt.rho2uvp(self.grid.angle,'p')
    else:
      ang=np.cos(self.grid.angle)+1j*np.sin(self.grid.angle)
      ang=np.angle(rt.rho2uvp(ang,'p'))

    u,v=calc.rot2d(u,v,-ang)

    out=outu
    out.label=out.label+'_uv'
    out.v=u+1j*v
    out.info['v']['name']=uname+vname

    if 'z' in coords: out.z=(out.z[1:,:]+out.z[:-1,:])/2.
    if 'x' in coords: out.x=(out.x[1:,:]+out.x[:-1,:])/2.
    if 'y' in coords: out.y=(out.y[1:,:]+out.y[:-1,:])/2.

    out.coordsReq=','.join(sorted(coords))
    return out


  def sliceiso(self,varname,iso,time,**opts):
    '''
    Depths where variable, increasing/dec with depth, has some value.

    Output is masked where surface is higher/lower than value or where all
    water column is lower/higher than value.
    '''


    coords=opts.get('coords',self._default_coords('sliceiso')).split(',')

    out=vis.Data()
    out.label='slice_iso'
    out.msg=self.check_slice(varname,t=time)
    if out.msg: return out

###    if not self.hasz(varname):
    if not 'z' in self.vaxes(varname):
      out.msg='a depth dependent variable is needed for sliceiso!'
      return out

    v=self.use(varname,SEARCHtime=time)
####    z=self.s_levels(time=time,ruvpw=self.var_at(varname))
###    z=self.s_levels(time=time,loc=self.var_at(varname))
    z=self.s_levels(time=time,loc=self.vloc(varname))
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
###      x,y,h=self.grid.vars(ruvp=self.var_at(varname))[:3]
      x,y,h=self.grid.vars(ruvp=self.vloc(varname))[:3]

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

    if 'z' in coords:
      # makes no sense... v is the depth!
      coords.remove('z')

    if 't' in coords and 't' in self.vaxes(varname): out.t=self.time[time]

    if out.v.ndim==2:
      out.extra=[vis.Data()]
      if 'x' in coords: out.extra[0].x=out.x
      if 'y' in coords: out.extra[0].y=out.y
      out.extra[0].v=h
      if h.max()>1000: cvals=200.,1000.
      elif h.max()>200: cvals=50.,100.,200.
      else: cvals=3
      out.extra[0].config['field.plot']='contour'
      out.extra[0].config['field.cvals']=cvals
      out.extra[0].config['field.cmap']='k'
      out.extra[0].label='bathy'

    out.coordsReq=','.join(sorted(coords))

    # about projection:
###    out.set_projection(self.grid.proj_info['basemap_opts'],extent=self.grid.proj_info['extent'])
    p=self.grid.get_projection(cartopy=1)
    out.set_projection(p,extent=self.grid.proj_info['extent'])

    return out


  def time_series(self,varname,x,y,times=None,depth=None,**opts):
    coords=opts.get('coords',self._default_coords('time_series')).split(',')

    if times is None: times=range(0,self.time.size)

    # depth or s_level: check if is float or if is negative!
    isDepth=False
    if not depth is None:
       if calc.isiterable(depth): depth=np.asarray(depth)
       if calc.isarray(depth):
         isDepth=np.any(depth<0) or depth.kind!='i' 
       else: isDepth=depth<0 or np.asarray(depth).dtype.kind!='i'

    out=vis.Data()
    out.label='time_series'
    if not depth is None and not isDepth:
      out.msg=self.check_slice(varname,t=np.max(times),k=depth) 
    else:
      out.msg=self.check_slice(varname,t=np.max(times)) 

    if out.msg: return out

    # find nearest point:
###    lon,lat,hr,mr=self.grid.vars(ruvp=self.var_at(varname))
    lon,lat,hr,mr=self.grid.vars(ruvp=self.vloc(varname))
    dist=(lon-x)**2+(lat-y)**2
    i,j=np.where(dist==dist.min())
    i,j=i[0],j[0]

    if not depth is None and not isDepth: arg={'s_SEARCH':depth}
    else: arg={}
    v=self.use(varname,xiSEARCH=j,etaSEARCH=i,SEARCHtime=times,**arg).T

    # calculate depths:
###    if self.hasz(varname):
    if 'z' in self.vaxes(varname):
      h=self.grid.h[i,j]
      zeta=self.use('zeta',xiSEARCH=j,etaSEARCH=i,SEARCHtime=times)
      h=h+0*zeta
####      z=rt.s_levels(h,zeta,self.s_params,rw=varname)
###      z=rt.s_levels(h,zeta,self.s_params,rw=self.var_at(varname)[1])
      z=rt.s_levels(h,zeta,self.s_params,rw=self.vloc(varname)[1])
      z=np.squeeze(z)

    # depth slice:
###    if isDepth and self.hasz(varname):
    if isDepth and 'z' in self.vaxes(varname):
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

          v=rt.slicez(v[...,np.newaxis],land_mask,
               self.grid.h[i,j]*np.ones((nt,1),dtype=v.dtype), # bottom depth
               zeta[:,np.newaxis],self.s_params,depth,
               surface_masked=opts.get('surf_mask',True),
               spline=opts.get('spline',True))[...,0]

      else: # one time only
        v=np.interp(depth,z,v,left=np.nan,right=np.nan)
        v=np.ma.masked_where(np.isnan(v),v) 


    out.v=v
    out.info['v']['name']=varname
    out.info['v']['slice']='time series'
    try: out.info['v']['units']=netcdf.vatt(self.nc,varname,'units')
    except: pass
     
    # coords
#########    if 't' in coords and self.hast(varname):
    if 't' in coords and 't' in self.vaxes(varname):
      if v.ndim==2:
        out.t=np.tile(self.time[times],(v.shape[0],1))
        from matplotlib.dates import date2num
        out.tnum=np.tile(date2num(self.time[times]),(v.shape[0],1))
      else: out.t=self.time[times]
      out.info['t']['name']='Time'
      out.info['tnum']=dict(name='Time',units=self.var_as['time']['units'])

###    if 'z' in coords and self.hasz(varname):
    if 'z' in coords and 'z' in self.vaxes(varname):
      if not depth is None:
        if not isDepth: out.z=z[depth,...]
        else: out.z=depth+0*v
      else: out.z=z
      out.info['z']=dict(name='Depth',units='m')

    if 'x' in coords:
      out.x=lon[i,j]
      if self.grid.spherical:
         out.info['x']=dict(name='Longitude',units=r'$\^o$E')
      else:
        out.x=x/1000.
        out.info['x']=dict(name='X-position',units='km')

    if 'y' in coords:
      out.y=lat[i,j]
      if self.grid.spherical:
        out.info['y']=dict(name='Latitude',units=r'$\^o$N')
      else:
        out.y=y/1000.
        out.info['y']=dict(name='Y-position',units='km')


    out.coordsReq=','.join(sorted(coords))
    return out


#  def extrap_at_mask(self,vname,time,**opts): 
#    method=opts.get('method','delaunay') # also available 'easy'
#    quiet=opts.get('quiet',True)
#
#    if method=='delaunay': extrap=calc.mask_extrap
#    elif method=='easy': extrap=calc.easy_extrap
#
#    if not quiet: print('Extrap at mask:\n  loading')
#    v=self.use(vname,SEARCHtime=time)
#    x,y,h,m=self.grid.vars(vname)
#
#    if 'z' in self.vaxes(vname): # 3d
#      N=v.shape[0]
#      for n in range(N):
#        if not quiet: print('  extrap level %d of %d' % (n,N))
#        v[n,...]=extrap(x,y,np.ma.masked_where(m==0,v[n,...]))
#    else: # 2d
#        v=extrap(x,y,np.ma.masked_where(m==0,v))
#
#    if np.ma.isMA(v) and v.count()==v.size: v=v.data # no need for masked array
#    return v

class MHis():
  def __init__(self,name,grd=False):
    import glob
    if cb.isstr(name):
      f=os.path.splitext(name)[0][:-1]+'*'+os.path.splitext(name)[1]
      files=glob.glob(f)
      files.sort()
    else:
      files=name

    if cb.isstr(grd):
      f=os.path.splitext(grd)[0][:-1]+'*'+os.path.splitext(grd)[1]
      gfiles=glob.glob(f)
      gfiles.sort()

    if grd:
      self._his=[]
      for i in range(len(files)):
        self._his+=[His(files[i],gfiles[i])]
    else:
      self._his=[His(i) for i in files]

    ###self.grid=MGrid([i.grid.name for i in self._his])
    self.grid=MGrid([i.grid for i in self._his])

  def __getitem__(self,i):
    return self._his[i]

  def __len__(self): return len(self._his)

  def __repr__(self):
    return '<MHis at 0x%x>'%id(self)+\
    '\n'+''.join([' -- '+os.path.basename(i.name)+'\n' for i in self._his])

  def tinds(self,it0):
    if it0>self[0].time.size-1:
      return [it0]*len(self),['bad it0','','']

    it=[it0]
    msg=['']
    for c,i in enumerate(self[1:]):
      if it0<i.time.size and i.time[it0]==self[0].time[it0]:
        it+=[it0]
        msg+=['']
      else:
        # find closest time index:
        it2 = np.argmin(np.abs(i.time -self[0].time[it0]))
        it+=[it2]
        if i.time[it2]!=self[0].time[it0]:
          msg+=['time %d differs from slice 0: %s'%(c+1,i.time[it2].isoformat())]

    return it,msg

  def _slice(self,slc,varname,ind,it0,**opts):
##    border=opts.get('border',1)

    if slc=='ll': x,y=ind

    o=[]
    It,twarn=self.tinds(it0)
    for c,i in enumerate(self):
      meth=getattr(i,'slice'+slc)
      if slc=='uv':  o+=[meth(ind,It[c],**opts)]
      elif slc=='ll':o+=[meth(varname,x,y,It[c],**opts)]
      else:          o+=[meth(varname,ind,It[c],**opts)]
      #if msg[c]:
      #  if o[-1].msg: o[-1].msg+='\n'
      #  o[-1].msg+=msg[c]


    # add filled border for next domain:
    if slc in ['i','j','k','z','iso']:#,'uv']:
      for i in range(len(self)-1):
        xb,yb=self[i+1].grid.border()
        o[i].extra+=[vis.Data(xb,yb)]
        o[i].extra[-1].label='border grid %d'%(i+1)
        #if i>0:
        o[i].extra[-1].config['d1.plot']='fill'
        o[i].extra[-1].config['d1_fill.options']['lw']=0
        o[i].extra[-1].config['d1_fill.options']['facecolor']='w'
        #o[i].extra[-1].config['plot.zorder']=1
        #else:
        #  o[i].extra[-1].config['d1.plot']='plot'
    elif slc=='uv': # for uv, better to mask arrows inside child domains
      if not hasattr(self.grid,'ingrd'): self.grid._set_ingrid('p')
      from functools import reduce
      for c,i in enumerate(o):
        if not i.v is None:
          if len(self.grid[c].ingrd_p):
            mask=reduce(lambda i,j: i&j,self.grid[c].ingrd_p)
            mask=mask&(~i.v.mask)
            i.nmask=mask # nested domains mask
            i.v.mask=i.v.mask|i.nmask
          else: i.nmask=None



    # vfield settings:
    if not o[0].v is None:
      vfield=o[0].get_param('vfield')
      if not vfield['options']['scale']:
        vfield['options']={'units':'width','scale':10}

      # use same vfield settings for all slices:
      for i in o:
        for k in vfield: i.config['vfield.'+k]=vfield[k]

      # no key for slices[1:]:
      for i in o[1:]:
        i.config['vfield.key_XYU']=0,0,0


    # field settings:
    if not o[0].v is None:
      field=o[0].get_param('field')

      if field['clim'] is False:
        if np.iscomplexobj(o[0].v):
          field['clim']=np.abs(o[0].v).min(),np.abs(o[0].v).max()
        else:
          field['clim']=o[0].v.min(),o[0].v.max()

        if field['clim'][0]==field['clim'][1]:
          field['clim']=field['clim'][0]-1,field['clim'][0]+1

      if field['cvals'] is False:
        tk=ticks.loose_label_n(field['clim'][0],field['clim'][1],7)
        field['cvals']=tk

      # use same field settings for all slices:
      for i in o:
        for k in field: i.config['field.'+k]=field[k]


    # also same settings for extras like bathy contours:
    Lab=['bathy']
    for lab in Lab:
      for e in o[0].extra:
        if e.label==lab and not e.v is None:
          field=e.get_param('field')
          if field['clim'] is False: field['clim']=e.v.min(),e.v.max()
          if field['cvals'] is False or  not calc.isiterable(field['cvals']):
            # the isiterable here is because it can be an integer (nof fixed set)
            tk=ticks.loose_label_n(field['clim'][0],field['clim'][1],3)
            field['cvals']=tk

      for i in o:
        for e in i.extra:
          if e.label==lab:
            for k in field: e.config['field.'+k]=field[k]

    # add border for all domains:
#    xb,yb=self.grid.border()
#    o[0].extra+=[vis.Data(x=xb,v=yb)]
    borders=[]
    for i in self:
      xb,yb=i.grid.border()
      c=vis.Data(x=xb,v=yb)
      c.set_param(d1_line__options=dict(lw=0.5,color='k',ls='-'))
      borders+=[c]

    o[-1].extra+=borders

    # set zorder:
    c=1
    for i in o:
      c+=0.1
      i.config['plot.zorder']=c
      for j in i.extra:
        c+=0.1
        j.config['plot.zorder']=c

    # place some stuff above continents:
    z=o[0].config['proj.continents']['zorder']
    z+=0.1
    # vfield:
    if slc=='uv':
      for i in o: i.config['plot.zorder']=z

    # place borders above continents:
    for b in borders: b.config['plot.zorder']=z


    return vis.MData(o,warnings=twarn)


  def slicek(self,varname,ind,time,**opts):
    return self._slice('k',varname,ind,time,**opts)
  def slicez(self,varname,ind,time,**opts):
    return self._slice('z',varname,ind,time,**opts)
  def sliceiso(self,varname,ind,time,**opts):
    return self._slice('iso',varname,ind,time,**opts)

#  def slicei(self,varname,ind,time,**opts):
#    opts['border']=0
#    return self._slice('i',varname,ind,time,**opts)

  def slicei(self,varname,ind,time,**opts):
#    opts['border']=0
    return self._slice('i',varname,ind,time,**opts)
  def slicej(self,varname,ind,time,**opts):
#    opts['border']=0
    return self._slice('j',varname,ind,time,**opts)

  def slicell(self,varname,X,Y,time,**opts):
#    opts['border']=0
    return self._slice('ll',varname,[X,Y],time,**opts)

  def sliceuv(self,ind,time,**opts):
#    opts['border']=0
    return self._slice('uv','',ind,time,**opts)

  def slice_derived(self,varname,ind,time,**opts):
    return self._slice('_derived',varname,ind,time,**opts)
