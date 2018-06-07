'''
tools for NOAA GNOME oil spill model
See also okean.gshhs for bna coast file generation

Martinho MA
'''

import numpy as np
import datetime
from okean import netcdf, calc

tunits='seconds since 2000-01-01  0:00:00 0:00'


def create_uv(fname,xi,eta):
  '''Surface currents file for GNOME'''

  nc=netcdf.ncopen(fname,'w',version=3)
  nc.add_dim('lon',xi)
  nc.add_dim('lat',eta)
  nc.add_dim('time',0)
  nc.close()
  nc=netcdf.ncopen(fname,'a')

  # lon
  v=nc.add_var('lon','double',('lat','lon'))
  v.add_att('long_name','Longitude')
  v.add_att('units','degrees_east')
  v.add_att('standard_name','longitude')

  # lat
  v=nc.add_var('lat','double',('lat','lon'))
  v.add_att('long_name','Latitude')
  v.add_att('units','degrees_north')
  v.add_att('standard_name','latitude')

  # depth
  v=nc.add_var('depth','double',('lat','lon'))
  v.add_att('long_name','Bathymetry')
  v.add_att('units','meters')
  v.add_att('positive','down')
  v.add_att('standard_name','depth')

  # mask
  v=nc.add_var('mask','double',('lat','lon'))
  v.add_att('long_name','Land Mask')
  v.add_att('units','nondimensional')
  v.add_att('standard_name','land_binary_mask')

  # u
  v=nc.add_var('u','double',('time','lat','lon'),fill_value=-99999.)
  v.add_att('long_name','Eastward Water Velocity')
  v.add_att('units','m/s')
  v.add_att('missing_value',-99999.)
  v.add_att('scale_factor',1.)
  v.add_att('add_offset',0.)
  v.add_att('standard_name','surface_eastward_sea_water_velocity')

  # v
  v=nc.add_var('v','double',('time','lat','lon'),fill_value=-99999.)
  v.add_att('long_name','Northward Water Velocity')
  v.add_att('units','m/s')
  v.add_att('missing_value',-99999.)
  v.add_att('scale_factor',1.)
  v.add_att('add_offset',0.)
  v.add_att('standard_name','surface_northward_sea_water_velocity')

  # time
  v=nc.add_var('time','double',('time',))
  v.add_att('long_name','Time')
  v.add_att('units',tunits)
  v.add_att('standard_name','time')

  nc.add_att('file_type','Full_Grid')
  nc.add_att('Conventions','COARDS')
  nc.add_att('grid_type','curvilinear')
  nc.add_att('z-type','s-coordinate')
  nc.add_att('model','ROMS')
  nc.add_att('title', 'forecast')
  nc.close()


def create_wind(fname,xi,eta):
  '''Surface wind file for GNOME'''

  nc=netcdf.ncopen(fname,'w',version=3)
  nc.add_dim('lon',xi)
  nc.add_dim('lat',eta)
  nc.add_dim('time',0)
  nc.close()
  nc=netcdf.ncopen(fname,'a')

  # lon
  v=nc.add_var('lon','double',('lat','lon'))
  v.add_att('long_name','Longitude')
  v.add_att('units','degrees_east')
  v.add_att('standard_name','longitude')
  v.add_att('point_spacing','even')

  # lat
  v=nc.add_var('lat','double',('lat','lon'))
  v.add_att('long_name','Latitude')
  v.add_att('units','degrees_north')
  v.add_att('standard_name','latitude')
  v.add_att('point_spacing','even')

  # depth
  v=nc.add_var('depth','double',('lat','lon'))
  v.add_att('long_name','Bathymetry')
  v.add_att('units','meters')
  v.add_att('positive','down')
  v.add_att('standard_name','depth')

  # u
  v=nc.add_var('air_u','double',('time','lat','lon'),fill_value=-9.9999e+32)
  v.add_att('long_name','Eastward Wind')
  v.add_att('units','m/s')
  v.add_att('scale_factor',1.)
  v.add_att('add_offset',0.)
  v.add_att('standard_name','eastward_wind')

  # v
  v=nc.add_var('air_v','double',('time','lat','lon'),fill_value=-9.9999e+32)
  v.add_att('long_name','Northward Wind')
  v.add_att('units','m/s')
  v.add_att('scale_factor',1.)
  v.add_att('add_offset',0.)
  v.add_att('standard_name','northward_wind')

  # time
  v=nc.add_var('time','double',('time',))
  v.add_att('long_name','Valid Time')
  v.add_att('units',tunits)
  v.add_att('standard_name','time')

  nc.add_att('Convetions','CF-1.3')
  nc.add_att('grid_type','curvilinear')
  nc.add_att('title', 'forecast')
  nc.close()


def his2gnome(fname,his,grd=False,nomask=False,gshhsMask=True,xylim=False,dates=False,ij=(1,1)):
  '''
  Creates GNOME wind file
  Ex:
    his2gnome(out,his,grd,dates=dates,ij=(2,2))

  if gshhsMask, the high res mask file mask_gshhs.npy will be created at 1st usage.
  Mask is based on high (h) resolution gshhs data which must be available (env variable
  GSHHS_MASK must be set). 
  '''

  if not grd: grd=his
  deta,dxi=ij

  dims=netcdf.fdim(his)
  xi,eta=dims['xi_rho'],dims['eta_rho']
  xi0,eta0=xi,eta

  nc0=netcdf.ncopen(his)
  time=netcdf.nctime(nc0,'ocean_time')
  # for roms agrif:
  #t=netcdf.use(nc0,'scrum_time')
  #time=netcdf.num2date(t,'seconds since %d-01-01' % year0)

  x0=netcdf.use(grd,'lon_rho')
  y0=netcdf.use(grd,'lat_rho')
  ang=netcdf.use(grd,'angle')

  if not xylim is False:
    xlim=xylim[:2]
    ylim=xylim[2:]
    i1,i2,j1,j2=calc.ij_limits(x0,y0,xlim,ylim)
    print(i1,i2,j1,j2)
    xi=i2-i1
    eta=j2-j1
  else:
    i1,i2=0,xi
    j1,j2=0,eta

  XI  ='%d:%d:%d' %(i1,i2,dxi)
  ETA ='%d:%d:%d' %(j1,j2,deta)

  xi=len(range(i1,i2,dxi))
  eta=len(range(j1,j2,deta))
  # create file:
  create_uv(fname,xi,eta)

  nc=netcdf.ncopen(fname,'a')
  for v0,v in ('lon_rho','lon'),('lat_rho','lat'),('mask_rho','mask'),('h','depth'):
    print('filling %s with %s' % (v,v0))
    nc.vars[v][:]=netcdf.use(grd,v0,xi_rho=XI,eta_rho=ETA)

  if nomask:
    print('NO MASK !!!')
    nc.vars['mask'][:]=1

  if gshhsMask:
    try:
     mask=np.load('mask_gshhs.npy')
    except:
      mask=1+0*netcdf.use(nc0,'mask_rho',xi_rho=XI,eta_rho=ETA)
      mask=mask.astype('bool')
      x=netcdf.use(grd,'lon_rho',xi_rho=XI,eta_rho=ETA)
      y=netcdf.use(grd,'lat_rho',xi_rho=XI,eta_rho=ETA)

      from okean import gshhs
      axis=x.min(),x.max(),y.min(),y.max()
      g=gshhs.gshhs(axis, resolution='h',area_thresh=0., max_level=2,clip=True)
      for lon, lat, level in zip(g.lon, g.lat, g.level):
        if level == 1: # land
          print('mask ',lon.shape)
          i=calc.inpolygon(x,y,lon,lat)
          mask=mask & ~i

      mask.dump('mask_gshhs.npy')


    nc.vars['mask'][:]=mask


  x=x0[j1:j2:deta,i1:i2:dxi]
  y=y0[j1:j2:deta,i1:i2:dxi]
  ang=ang[j1:j2:deta,i1:i2:dxi]

  n=-1
  for it in range(len(time)):
    if not dates is False:
      d0,d1=dates
      if time[it]<d0 or time[it]>=d1: continue

    n+=1
    U=np.zeros((eta0,xi0),'f')
    V=np.zeros((eta0,xi0),'f')

    nc.vars['time'][n]=netcdf.date2num(time[it],tunits)

    # for roms agrif:
    #u=netcdf.use(nc0,'u',time=it,s_rho=-1)
    #v=netcdf.use(nc0,'v',time=it,s_rho=-1)
    u=netcdf.use(nc0,'u',ocean_time=it,s_rho=-1)
    v=netcdf.use(nc0,'v',ocean_time=it,s_rho=-1)

    # mask extrap:
    print('mask extrap...')

    u=calc.mask_extrap(x0,y0,np.ma.masked_where(u==0,u))
    v=calc.mask_extrap(x0,y0,np.ma.masked_where(v==0,v))

    U[:,1:-1]=0.5*(u[:,:-1]+u[:,1:])
    U[:,0]=u[:,0]
    U[:,-1]=u[:,-1]

    V[1:-1,:]=0.5*(v[:-1.:]+v[1:,:])
    V[0,:]=v[0,:]
    V[-1,:]=v[-1,:]

    U=U[j1:j2,i1:i2]
    V=V[j1:j2,i1:i2]
  
    U=U[j1:j2:deta,i1:i2:dxi]
    V=V[j1:j2:deta,i1:i2:dxi]

    # rotate uv:
    print('rotating ...')
    U,V=calc.rot2d(U,V,-ang)

    print('filling uv', n, time[it])
    nc.vars['u'][n,...]=U
    nc.vars['v'][n,...]=V

  nc.close()
  nc0.close()


def frc2gnome(fname,frc,grd,xylim=False,dates=False,ij=(1,1),**kargs):
  '''
  Creates GNOME wind file
  kargs:
    t[u,v]var
    t[u,v]dim
    x[y,ang]var

  Ex:
    .frc2gnome(out,frc,grd,ij=(10,10),dates=dates,**{'tdim':'Time'})
  '''

  deta,dxi=ij

  tvar='time'
  uvar='Uwind'
  vvar='Vwind'
  #tvar='bulk_time'
  #uvar='uwnd'
  #vvar='vwnd'

  tdim='time'
  #tdim='bulk_time'
  xdim='xi_rho'
  ydim='eta_rho'

  xvar='lon_rho'
  yvar='lat_rho'
  angvar='angle'

  if 'tvar' in kargs.keys(): tvar=kargs['tvar']
  if 'uvar' in kargs.keys(): uvar=kargs['uvar']
  if 'vvar' in kargs.keys(): vvar=kargs['vvar']

  if 'tdim' in kargs.keys(): tdim=kargs['tdim']
  if 'xdim' in kargs.keys(): xdim=kargs['xdim']
  if 'ydim' in kargs.keys(): ydim=kargs['ydim']

  if 'xvar' in kargs.keys(): xvar=kargs['xvar']
  if 'yvar' in kargs.keys(): yvar=kargs['yvar']
  if 'angvar' in kargs.keys(): angvar=kargs['angvar']


  dims=netcdf.fdim(grd)
  xi,eta=dims[xdim],dims[ydim]
  xi0,eta0=xi,eta

  ncg=netcdf.ncopen(grd)

  nc0=netcdf.ncopen(frc)
  try:
   t=netcdf.nctime(nc0,tvar)
  except:
    t=netcdf.use(nc0,tvar)
    t=netcdf.num2date(t,'days since %d-01-01' % year0)

  time=netcdf.date2num(t,tunits)

  x0=netcdf.use(grd,xvar)
  y0=netcdf.use(grd,yvar)
  if x0.ndim==1: x0,y0=np.meshgrid(x0,y0)

  if angvar:
    ang=netcdf.use(grd,angvar)

  if not xylim is False:
    xlim=xylim[:2]
    ylim=xylim[2:]
    i1,i2,j1,j2=calc.ij_limits(x0,y0,xlim,ylim)
    xi=i2-i1
    eta=j2-j1
  else:
    i1,i2=0,xi
    j1,j2=0,eta

  XI  ='%d:%d:%d' %(i1,i2,dxi)
  ETA ='%d:%d:%d' %(j1,j2,deta)

  xi=len(range(i1,i2,dxi))
  eta=len(range(j1,j2,deta))

  # create file:
  create_wind(fname,xi,eta)

  nc=netcdf.ncopen(fname,'a')

  x=x0[j1:j2:deta,i1:i2:dxi]
  y=y0[j1:j2:deta,i1:i2:dxi]

  nc.vars['lon'][:]=x
  nc.vars['lat'][:]=y
  if angvar: ang=ang[j1:j2:deta,i1:i2:dxi]

  n=-1
  for it in range(len(time)):

    if not dates is False:
      d0,d1=dates
      if t[it]<d0 or t[it]>=d1: continue

    n+=1
    u=netcdf.use(nc0,uvar,**{xdim:XI,ydim:ETA,tdim:it})
    v=netcdf.use(nc0,vvar,**{xdim:XI,ydim:ETA,tdim:it})

    # rotate uv:
    if angvar:
      print('rotating ...')
      u,v=calc.rot2d(u,v,-ang)


    nc.vars['time'][n]=time[it]
    print('filling uv',n,t[it])
    nc.vars['air_u'][n,...]=u
    nc.vars['air_v'][n,...]=v


  nc.close()
  nc0.close()
  ncg.close()

class GnomeOut:
  '''
  Reads GNOME output for each time
  Ex:
    GnomeOut('moss005')
  will use the files moss005.ms3..6
  '''

  def __init__(self,name):
    self.ms3=name+'.ms3'
    self.ms4=name+'.ms4'
    self.ms5=name+'.ms5'
    self.ms6=name+'.ms6'

    print(' -- loading from %s' % name)

    self.get_date()
    self.get_pos()

  def get_date(self):
   tmp=open(self.ms3).readlines()[4].strip().split('VALIDFOR:')[1].strip()
   self.date=datetime.datetime.strptime(tmp,'%H:%M, %m/%d/%y')

  def get_pos(self):
    # best estimate solution:
    tmp=open(self.ms4).readlines()[1::2]

    self.x=np.zeros(len(tmp),'f')
    self.y=np.zeros(len(tmp),'f')
    self.z=np.zeros(len(tmp),'f')

    for i in range(len(tmp)):
      self.x[i],self.y[i],self.z[i]=[float(k) for k in tmp[i].split()]

    # minimum regret solution:
    if os.path.isfile(self.ms6):
      tmp=open(self.ms6).readlines()[1::2]

      self.xr=np.zeros(len(tmp),'f')
      self.yr=np.zeros(len(tmp),'f')
      self.zr=np.zeros(len(tmp),'f')

      for i in range(len(tmp)):
        self.xr[i],self.yr[i],self.zr[i]=[float(k) for k in tmp[i].split()]

    else: self.xr,self.yr,self.zr=[],[],[]


class Gnome:
  '''
  Get data from gnome forecasts

  Ex: 
  Gnome('moss%04d',ndays=5)
  will use the files moss000.ms3..6, etc
  '''

  def __init__(self,name,ndays,dt=1):
    NT=ndays*24/dt+1
    for i in range(NT):
      F=name%i

      g=GnomeOut(F)
      if i==0:
        self.x=np.zeros((NT,len(g.x)),'f')
        self.y=np.zeros((NT,len(g.y)),'f')
        self.dates=np.zeros(NT,dtype=datetime.datetime)

        self.xr=np.zeros((NT,len(g.xr)),'f')
        self.yr=np.zeros((NT,len(g.yr)),'f')

      self.x[i,:g.x.size]=g.x
      self.y[i,:g.y.size]=g.y
      self.dates[i]=g.date
      if len(g.xr):
        self.xr[i,:g.xr.size]=g.xr
        self.yr[i,:g.yr.size]=g.yr
