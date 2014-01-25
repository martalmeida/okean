import numpy as np
from ompy.nc import netcdf
from ompy.roms import roms
import os
from ompy.tools import date_tools as dt
from ompy.tools import cookbook as cb
from ompy.tools import num_tools as nt
from ompy.tools import air_sea
import datetime

grd='http://nomads.ncdc.noaa.gov/thredds/dodsC/narr-a/AWIP32-fixed.grb'

# cant find grid in opendap anymore !!
grd='AWIP32-fixed.grb'
grdTxt='latlon.g221.txt'
from okean.cookbook import run
if not os.path.isfile(grd):
  run('wget','http://nomads.ncdc.noaa.gov/data/narr/AWIP32-fixed.grb')

# pygrib cannot read the grid vars !!! so, will use the lon x lat txt
# file form the http server!
if not os.path.isfile(grdTxt):
  run('wget','http://nomads.ncdc.noaa.gov/data/narr/latlon.g221.txt')

#baseurl = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/'
baseurl = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr-a/'

narrdt  = 3 #hours

def load_grid():
  sname='latlon.g221.npz'
  if os.path.isfile(sname):
    a=np.load(sname)
    return a['lon'],a['lat']
  else:
    lines=open(grdTxt).readlines()[5:]
    lon=np.zeros(len(lines),'f')
    lat=np.zeros(len(lines),'f')
    N=0
    c=-1
    for l in lines:
      c+=1
      l=l.split()
      lat[c]=float(l[2])
      lon[c]=float(l[3])
      if N==0 and int(l[1])==2: N=c

    M=lon.size/N
    lon.shape=M,N
    lat.shape=M,N
    np.savez(sname,lon=lon,lat=lat)

  return lon,lat

class NARRData:

  def __init__(self): pass

  def files(self,date0,date1=False):
    date0=dt.parse_date(date0)
    if not date1: date1=date0+datetime.timedelta(1)
    else: date1=dt.parse_date(date1)

    dates=dt.drange(date0,date1,True)

    files=[]
    time=[]
    for d in dates:
      if d==dates[-1]: hours=[0]
      else: hours=range(0,24,narrdt)

      for hour in hours:
        #fname='narr-a_%d_%s_%02d00_000.grb' % (int(d.strftime('%j'))-1,d.strftime('%Y%m%d'),hour)
        fname='narr-a_221_%s_%02d00_000.grb' % (d.strftime('%Y%m%d'),hour)
        files+=[os.path.join(baseurl,d.strftime('%Y%m'),d.strftime('%Y%m%d'),fname)]
        time+=[d+datetime.timedelta(hour/24.)]

    return files, np.array(time)


  def data(self,date0,date1=False,xlim=False,ylim=False,quiet=True):
    files,time=self.files(date0,date1)
    res=cb.odict()

    for i in range(len(files)):
      if not quiet: print '|-> getting from file %s' % files[i]

      try:
        data=narr_file_data(files[i],xlim,ylim,quiet=quiet)
      except:
        if not quiet: print 'CANNOT use data for %s' % time[i].isoformat()
        continue

      if data:
        data['INFO_file']   = files[i]
        res[time[i]]=data
      else:
        if not quiet: print 'NO DATA for %s' % time[i].isoformat()

    miss={}
    return res,miss

def narr_file_data(fname,xlim=False,ylim=False,quiet=False):
  '''
  Returns bulk data from one NARR file
  '''

  out={}

  # loading grid:
  if 0:
    if not quiet: print ' reading lon,lat from file %s' % grd
    nc=netcdf.ncopen(grd)
    x=nc.vars['East_longitude_0-360'][0,...]-360.
    y=nc.vars['Latitude_-90_to_+90'][0,...] # time always 1 !!
    nc.close()
  else:
    if not quiet: print ' reading lon,lat from file %s' % grdTxt
    x,y=load_grid()
    #x=x-360.
    x=-x

  ny,nx=x.shape


  if (xlim,ylim)==(False,False):i0,i1,j0,j1=0,nx,0,ny
  else:
    i0,i1,j0,j1=nt.ij_limits(x, y, xlim, ylim, margin=0)
    x=x[j0:j1,i0:i1]
    y=y[j0:j1,i0:i1]

  try:
    nc=netcdf.ncopen(fname)
  except:
    return {}

  xx=str(i0)+':'+str(i1)
  yy=str(j0)+':'+str(j1)

  tdim=netcdf.fdim(nc,'time1')
  if tdim!=1: print 'WARNING: tdim !=1  !!!!!!'

  # T surface [K->C]
  if not quiet: print ' --> T air'
  tair=netcdf.use(nc,'Temperature_surface',time1=0,x=xx,y=yy)
  tair=tair-273.15
  out['tair']=cb.Data(x,y,tair,'C')

  # R humidity [% -> 0--1]
  if not quiet: print ' --> R humidity'
  rhum=netcdf.use(nc,'Relative_humidity',time1=0,x=xx,y=yy)
  out['rhum']=cb.Data(x,y,rhum/100.,'0--1')

  # surface pressure [Pa]
  if not quiet: print ' --> Surface pressure'
  pres=netcdf.use(nc,'Pressure_surface',time1=0,x=xx,y=yy)
  out['pres']=cb.Data(x,y,pres,'Pa')

  # P rate [kg m-2 s-1 -> cm/d]
  if not quiet: print ' --> P rate'
  prate=netcdf.use(nc,'Precipitation_rate',time1=0,x=xx,y=yy)
  prate=prate*86400*100/1000.
  out['prate']=cb.Data(x,y,prate,'cm/d')

  # Net shortwave flux  [ W m-2]
  if not quiet: print ' --> Net shortwave flux'
  if not quiet: print '       SW down'
  sw_down=netcdf.use(nc,'Downward_shortwave_radiation_flux',time1=0,x=xx,y=yy)
  if not quiet: print '       SW up'
  sw_up=netcdf.use(nc,'Upward_short_wave_radiation_flux_surface',time1=0,x=xx,y=yy)
  sw_net=sw_down-sw_up
  out['radsw']=cb.Data(x,y,sw_net,'W m-2',info='positive downward')

  # Net longwave flux  [W/m^2]
  if not quiet: print ' --> Net longwave flux'
  if not quiet: print '       LW down'
  lw_down=netcdf.use(nc,'Downward_longwave_radiation_flux',time1=0,x=xx,y=yy)
  if not quiet: print '       LW up'
  lw_up=netcdf.use(nc,'Upward_long_wave_radiation_flux_surface',time1=0,x=xx,y=yy)
  lw_net=lw_down-lw_up
  out['radlw']=cb.Data(x,y,-lw_net,'W m-2',info='positive upward')

  # downward lw:
  out['dlwrf']=cb.Data(x,y,-lw_down,'W m-2',info='negative... downward')

  # U and V wind speed 10m
  if not quiet: print ' --> U and V wind'
  # vertical dim is height_above_ground1: 10 and 30 m
  uwnd=netcdf.use(nc,'u_wind_height_above_ground',height_above_ground1=0,time1=0,x=xx,y=yy)
  vwnd=netcdf.use(nc,'v_wind_height_above_ground',height_above_ground1=0,time1=0,x=xx,y=yy)

  if not quiet: print ' --> calc wind speed and stress'
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=cb.Data(x,y,speed,'m s-1')
  out['uwnd']=cb.Data(x,y,uwnd,'m s-1')
  out['vwnd']=cb.Data(x,y,vwnd,'m s-1')
  out['sustr']=cb.Data(x,y,taux,'Pa')
  out['svstr']=cb.Data(x,y,tauy,'Pa')

  # Cloud cover [0--100 --> 0--1]:
  if not quiet: print ' --> Cloud cover'
  clouds=netcdf.use(nc,'Total_cloud_cover',time1=0,x=xx,y=yy)
  out['cloud']=cb.Data(x,y,clouds/100.,'fraction (0--1)')

  nc.close()
  return  out

