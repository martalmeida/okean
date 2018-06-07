import os
import glob
from okean import dateu, netcdf, roms, calc, cookbook as cb
import numpy as np
import datetime



def source(date):
  src=glob.glob('%d/aoscat_%d_*.nc'%(date.year,date.year))
  src.sort()
  # cant return all files cos MFDataset wont work with files
  # without unlimited dimension !
  for s in src:
    t=netcdf.nctime(s,'time')
    if date>=t[0] and date<=t[-1]: return s



def read_wind(grd,date,ij=False):
  f=source(date)
  print('-- reading from %s'%f)
  time=netcdf.nctime(f,'time')

  try:
    i=np.where(time==date)[0][0]
  except:
    return 'date %s not found'%date.isoformat(' ')

  returnXY=False
  if ij is False:
    returnXY=True
    lon=netcdf.use(f,'longitude') # -180..180
    lat=netcdf.use(f,'latitude')
    g=roms.Grid(grd)
    xl0=np.asarray((g.lon.min(),g.lon.max()))
    xl=np.asarray((g.lon.min(),g.lon.max()))
    if np.any(xl>180) or np.any(xl<-180):
      print('ERROR: grid is supposed to be -180<x<180')
      print('Can be implemented with mpl_toolkits.basemap.shiftgrid ... TODO')
      print('(http://matplotlib.org/basemap/api/basemap_api.html)')
      return

    yl=g.lat.min(),g.lat.max()
    ij=calc.ij_limits(lon,lat,xl,yl,margin=1)

  i0,i1,j0,j1=ij
  u=netcdf.use(f,'eastward_wind',longitude='%d:%d'%(i0,i1),latitude='%d:%d'%(j0,j1),time=i)
  v=netcdf.use(f,'northward_wind',longitude='%d:%d'%(i0,i1),latitude='%d:%d'%(j0,j1),time=i)
  if returnXY:
    lon=netcdf.use(f,'longitude',longitude='%d:%d'%(i0,i1),latitude='%d:%d'%(j0,j1))
    lat=netcdf.use(f,'latitude',longitude='%d:%d'%(i0,i1),latitude='%d:%d'%(j0,j1))

    lon,lat=np.meshgrid(lon,lat)
    #if np.all(xl0<0): lon=lon-360 # this may be wrong ... if xl is near 0, lon ay have pos and neg values !!! fix this one day ...
    return lon,lat,u,v, ij
  else: return u,v


def get_wind(grd,date0,date1):
  dates=[date0]

  while dates[-1]<date1:
    dates+=[dates[-1]+datetime.timedelta(1/4.)] # 6h data

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


def gen_frc(fname,grd,tag='_aoscat06'):

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
  nc.add_att('title','ASCAT/OSCAT 06h wind')
  nc.add_att('info','myocean GLO-BLENDED_WIND_L4-V3-OBS_FULL_TIME_SERIE')
  nc.add_att('grd_file',os.path.realpath(grd))
  from time import ctime
  nc.add_att('history','ROMS  wind  file, '+ctime())
  nc.add_att('author',cb.username()[1]+', '+cb.machinename())

  nc.close()


def fill_frc(fname,time,u,v):
  nc=netcdf.Pync(fname,'w')

  tunits=netcdf.vatt(nc,'wind_time','units')
  time=netcdf.date2num(time,tunits)
  tind=nc.dims['wind_time']

  nc.vars['wind_time'][tind]=time
  nc.vars['Uwind_aoscat06'][tind,...]=u
  nc.vars['Vwind_aoscat06'][tind,...]=v

  nc.close()


def make_frc(frcname,grd,date0,date1):
  g=roms.Grid(grd)

  if not os.path.isfile(frcname):
    # create file:
    gen_frc(frcname,grd)
  else:
    last=netcdf.nctime(frcname,'wind_time',wind_time=-1)
    print('-->found file %s with last time %s'%(frcname,last.isoformat()))
    date0=last+datetime.timedelta(1/4.) # 6h data

  dates=[date0]
  while dates[-1]<date1:
    dates+=[dates[-1]+datetime.timedelta(1/4.)] # 6h data

  n=-1
  ij=False
  for d in dates:
    n+=1
    if ij is False:
      tmp=read_wind(grd,d,ij=False)
      if isinstance(tmp,basestring): continue # some error, like date not found!
      else: x,y,u,v,ij=tmp
    else:
      tmp=read_wind(grd,d,ij)
      if isinstance(tmp,basestring): continue
      else: u,v=tmp

    U=calc.griddata(x,y,u,g.lon,g.lat,extrap=True)
    V=calc.griddata(x,y,v,g.lon,g.lat,extrap=True)

     # rotate wind,
    print(' --> rot U,V')
    angle=g.use('angle')
    U,V=calc.rot2d(U,V,angle)

    print('  filling %s'%d.isoformat(' '))
    fill_frc(frcname,d,U,V)


if __name__=='__main__':
  import sys
  if len(sys.argv)==5:
    grd=sys.argv[1]
    fname=sys.argv[2]
    date0=dateu.parse_date(sys.argv[3])
    date1=dateu.parse_date(sys.argv[4])
    make_frc(fname,grd,date0,date1)
  else:
     print('USAGE: python aoscat.py roms_grd.nc  wind_frc.nc 20100101 20110101')

