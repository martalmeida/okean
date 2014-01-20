'''
http://podaac.jpl.nasa.gov/dataset/CCMP_MEASURES_ATLAS_L4_OW_L3_0_WIND_VECTORS_FLK
'''
import os
import numpy as np
from okean import netcdf, roms, calc, dateu, cookbook as cb
import datetime
import netcdftime


def source(date):
  url0='http://podaac-opendap.jpl.nasa.gov/opendap/allData/ccmp/L3.0/flk/'
  return os.path.join(url0,'%d'%date.year,'%02d'%date.month,'analysis_%s_v11l30flk.nc.gz'%date.strftime('%Y%m%d'))

def read_wind(grd,date,ij=False):
  f=source(date)
  print '-- reading from %s'%f
  time=netcdf.nctime(f,'time')
  try:
    i=np.where(time==date)[0][0]
  except:
    return 'date %s not found'%d.isoformat(' ')

  returnXY=False
  if ij is False:
    returnXY=True
    lon=netcdf.use(f,'lon')
    lat=netcdf.use(f,'lat')
    g=roms.Grid(grd)
    xl0=np.asarray((g.lon.min(),g.lon.max()))
    xl=np.asarray((g.lon.min(),g.lon.max()))
    if np.all(xl<0): xl=xl+360
    elif np.any(xl<0) and np.any(xl>0):
      print 'ERROR: zero crossing not implemented !!!'
      print 'can be done with mpl_toolkits.basemap.shiftgrid ... TODO'
      print '(http://matplotlib.org/basemap/api/basemap_api.html)'
      return

    yl=g.lat.min(),g.lat.max()
    ij=calc.ij_limits(lon,lat,xl,yl,margin=1)

  i0,i1,j0,j1=ij
  u=netcdf.use(f,'uwnd',lon='%d:%d'%(i0,i1),lat='%d:%d'%(j0,j1),time=i)
  v=netcdf.use(f,'vwnd',lon='%d:%d'%(i0,i1),lat='%d:%d'%(j0,j1),time=i)
  if returnXY:
    lon=netcdf.use(f,'lon',lon='%d:%d'%(i0,i1),lat='%d:%d'%(j0,j1))
    lat=netcdf.use(f,'lat',lon='%d:%d'%(i0,i1),lat='%d:%d'%(j0,j1))
    lon,lat=np.meshgrid(lon,lat)
    if np.all(xl0<0): lon=lon-360 # this may be wrong ... if xl is near 0, lon ay have pos and neg values !!! fix this one day ...
    return lon,lat,u,v, ij
  else: return u,v


def get_wind(grd,date0,date1):
  dates=[date0]

  while dates[-1]<date1:
    dates+=[dates[-1]+datetime.timedelta(1/4.)]

  n=-1
  for d in dates:
    n+=1
    if d==dates[0]:
      x,y,u,v,ij=read_wind(grd,d,ij=False)
      U=np.zeros((len(dates))+u.shape,u.dtype)
      V=np.zeros((len(dates))+v.shape,v.dtype)
    else:
      u,v=read_wind(grd,d,ij)

    U[n,...]=u
    V[n,...]=v


def gen_frc(fname,grd,tag='_ccmp'):

  nc=netcdf.Pync(fname,'t',version=3)

  # dims:
  grd_dims=netcdf.fdim(grd)
  gdims='xi_rho','xi_u','xi_v','eta_rho','eta_u','eta_v'
  for name in gdims: nc.add_dim(name,grd_dims[name])
  # time dim:
  nc.add_dim('wind_time',0)

  v=nc.add_var('wind_time',np.dtype('d'),('wind_time',))
  v.add_att('long_name','wind forcing time')
  tunits='days since 1970-01-01'
  v.add_att('units',tunits)

  v=nc.add_var('Uwind'+tag,np.dtype('d'),('wind_time','eta_rho', 'xi_rho'))
  v.add_att('long_name','u-wind')
  v.add_att('units','metre second-1')
  v.add_att('time','time')

  v=nc.add_var('Vwind'+tag,np.dtype('d'),('wind_time','eta_rho', 'xi_rho'))
  v.add_att('long_name','v-wind')
  v.add_att('units','metre second-1')
  v.add_att('time','time')

  # Global Attributes:
  nc.add_att('type','Wind forcing file')
  nc.add_att('title','CCMP wind')
  nc.add_att('grd_file',os.path.realpath(grd))
  from time import ctime
  nc.add_att('history','ROMS  wind  file, '+ctime())
  nc.add_att('author',cb.username()[1]+', '+cb.machinename())

  nc.close()


def fill_frc(fname,time,u,v):

  nc=netcdf.Pync(fname,'w')

  tunits=netcdf.vatt(nc,'wind_time','units')
  time=netcdftime.date2num(time,tunits)
  tind=nc.dims['wind_time']

  nc.vars['wind_time'][tind]=time
  nc.vars['Uwind_ccmp'][tind,...]=u
  nc.vars['Vwind_ccmp'][tind,...]=v

  nc.close()

def make_frc(frcname,grd,date0,date1):
  g=roms.Grid(grd)

  if not os.path.isfile(frcname):
    # create file:
    gen_frc(frcname,grd)
  else:
    last=netcdf.nctime(frcname,'wind_time',wind_time=-1)
    print '-->found file %s with last time %s'%(frcname,last.isoformat())
    date0=last+datetime.timedelta(1/4.)

  dates=[date0]
  while dates[-1]<date1:
    dates+=[dates[-1]+datetime.timedelta(1/4.)]

  n=-1
  for d in dates:
    n+=1
    if d==dates[0]:
      x,y,u,v,ij=read_wind(grd,d,ij=False)
    else:
      u,v=read_wind(grd,d,ij)

    U=calc.griddata(x,y,u,g.lon,g.lat)
    V=calc.griddata(x,y,v,g.lon,g.lat)

     # rotate wind,
    print ' --> rot U,V'
    angle=g.use('angle')
    U,V=calc.rot2d(U,V,angle)

    print '  filling %s'%d.isoformat(' ')
    fill_frc(frcname,d,U,V)

if __name__=='__main__':
  import sys
  if len(sys.argv)<5:
    print 'USAGE: python ccmp.py roms_grd.nc  wind_frc.nc 20100101 20110101'
  else:
    grd=sys.argv[1]
    fname=sys.argv[2]
    date0=dateu.parse_date(sys.argv[3])
    date1=dateu.parse_date(sys.argv[4])
    make_frc(fname,grd,date0,date1)





