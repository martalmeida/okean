import pylab as pl
import numpy as np
from okean import air_sea, dateu, netcdf


class Data:
  '''
  Very simple class for x,y,v data.
  Includes the attributes x,y and data.
  By default, also includes units and info.
  '''

  def __init__(self,x,y,data,units,info=''):
    self.x     = x
    self.y     = y
    self.data  = data
    self.units = units
    self.info  = info


class CORDEXData:
  '''
  CORDEX clim data extraction
  '''
  def __init__(self,climFile):
    '''
    climFile is a climatology of CORDEX data and includes:
      - lon, lat time
      - clouds (clt), may not be present in some CORDEX datasets
      - humid (huss)
      - sw_down (rsds)
      - sw_up (rsus)
      - lw_down (rlds)
      - lw_up (rlus); in some datasets this variable is negative, but
          in climFile is is always positive! All lw and sw are positive!
      - rain (pr)
      - tair (tas)
      - u (uas); 10m u wind
      - v (vas); 10m v wind

    The file is expected to have daily (but may not be) data and the
    extremes are repeated by default; ie, the last values are added
    before the 1st and the 1st is added after the last. This is done
    to ensure data at time <= 00h 1-jan and at time>=00h 1-jan year+1.
    '''

    self.f=climFile

  def data(self,date0=False,date1=False,quiet=True,grd=False):
    '''
    Returns atm data date0 (>=) and date1 (<=)
    '''
    if grd:
      from okean import roms
      r=roms.Grid(grd)
      xlim=r.lon.min(),r.lon.max()
      ylim=r.lat.min(),r.lat.max()
      lims=xlim,ylim
    else:
      lims=False

    # get data for all times in file:
    res=cordex_file_data(self.f,lims=lims,quiet=quiet)
    time=res['time']
    if date0:
      date0=dateu.parse_date(date0)
      i,=np.where(time>=date0)
      i=i[0]
    else: i=0

    if date1:
      date1=dateu.parse_date(date1)
      j,=np.where(time<=date1)
      j=j[-1]
    else: j=len(time)

    if date0 or date1:
      for k in res:
        if k =='time':
          res[k]=res[k][i:j+1]
        else:
          res[k].data=res[k].data[i:j+1,...]

    return res


def fill_extremes(data,quiet=False):
  for k in data:
    if not quiet: print('-- repeating extremes for %s'%k)
    if k=='time':
      t=data[k]
      t=np.hstack((t[-1],t,t[0]))
      t[0]=t[0].replace(year=t[0].year-1)
      t[-1]=t[-1].replace(year=t[-1].year+1)
      data[k]=t
    else:
      v=data[k].data
      data[k].data=np.vstack((v[-1][np.newaxis],v,v[0][np.newaxis]))

  return data


def cordex_file_data(f,lims=False,quiet=False):
  '''
  CORDEX data for ROMS

  No accumulated variables are considered
  '''

  out={}

  # time, x, y:
  if not quiet: print(' reading time,x,y')
  out['time']=netcdf.nctime(f,'time')
  x=netcdf.use(f,'lon')
  y=netcdf.use(f,'lat')
  x[x>180]=x[x>180]-360
  if x.ndim==1 and y.ndim==1:
    x,y=np.meshgrid(x,y)

  if np.ma.isMA(x): x=x.data
  if np.ma.isMA(y): y=y.data

  if lims:
    from okean import calc
    xlim,ylim=lims
    i1,i2,j1,j2=calc.ij_limits(x,y,xlim,ylim,margin=3)
  else:
    i0=0
    j0=0
    j1,i1=x.shape

  I=range(i1,i2)
  J=range(j1,j2)
  x=x[j1:j2,i1:i2]
  y=y[j1:j2,i1:i2]

  # tair [K-->C]
  if not quiet: print(' --> T air')
  vname='tair'
  tair=netcdf.use(f,vname,lon=I,lat=J)-273.15
  out['tair']=Data(x,y,tair,'Celsius')

  # R humidity [0--1]
  if not quiet: print(' --> R humidity (from specific humidity)')
  vname='humid'
  q=netcdf.use(f,vname,lon=I,lat=J) # specific humidity
  rhum=q/air_sea.qsat(tair)
  rhum[rhum>1]=1
  out['rhum']=Data(x,y,rhum,'0--1')

  # surface pressure [Pa]
  if not quiet: print(' --> Surface pressure')
  vname='press'
  pres=netcdf.use(f,vname,lon=I,lat=J)
  out['pres']=Data(x,y,pres,'Pa')

  # P rate [kg m-2 s-1 -> cm/d]
  if not quiet: print(' --> P rate')
  vname='rain'
  prate=netcdf.use(f,vname,lon=I,lat=J)
  prate=prate*86400*100/1000.
  prate[prate<0]=0
  out['prate']=Data(x,y,prate,'cm/d')

  # Net shortwave flux  [W m-2]
  if not quiet: print(' --> Net shortwave flux')
  if not quiet: print('       SW down')
  sw_down=netcdf.use(f,'sw_down',lon=I,lat=J)
  if not quiet: print('       SW up')
  sw_up=netcdf.use(f,'sw_up',lon=I,lat=J)
  sw_net=sw_down-sw_up
  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')

  # Net longwave flux  [W m-2]
  if not quiet: print(' --> Net longwave flux')
  if not quiet: print('       LW down')
  lw_down=netcdf.use(f,'lw_down',lon=I,lat=J)
  if not quiet: print('       LW up')
  lw_up=netcdf.use(f,'lw_up',lon=I,lat=J)
  lw_net=lw_down-lw_up
  out['radlw']=Data(x,y,-lw_net,'W m-2',info='positive upward')
  # downward lw:
  out['dlwrf']=Data(x,y,-lw_down,'W m-2',info='negative... downward')
  # signs convention is better explained in wrf.py

  # U and V wind speed 10m
  if not quiet: print(' --> U and V wind')
  uwnd=netcdf.use(f,'u',lon=I,lat=J)
  vwnd=netcdf.use(f,'v',lon=I,lat=J)
  if not quiet: print(' --> calc wind speed and stress')
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=Data(x,y,speed,'m s-1')
  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  out['sustr']=Data(x,y,taux,'Pa')
  out['svstr']=Data(x,y,tauy,'Pa')

  # Cloud cover [0--100 --> 0--1]:
  if not quiet: print(' --> Cloud cover')
  if 'clouds' in netcdf.varnames(f):
    clouds=netcdf.use(f,'clouds',lon=I,lat=J)/100.
    out['cloud']=Data(x,y,clouds,'fraction (0--1)')
  else:
    print('==> clouds not present!')

  return fill_extremes(out,quiet)
