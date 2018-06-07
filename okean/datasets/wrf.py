import numpy as np
from okean import air_sea, dateu, netcdf
import datetime

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


class WRFData:
  '''GFS data extraction'''

  def __init__(self,basefolder,files='wrfout*'):
    self.basefolder=basefolder

    import glob
    import os
    self.files=glob.glob(os.path.join(self.basefolder,files))
    self.files.sort()


  def data(self,date0=False,date1=False,quiet=True):
    '''
    Returns atm data form all times in basefolder files or between
    date0 (>=) and date1 (<=)
   '''
    if date0 is False: date0=datetime.datetime(1,1,1)
    else: date0=dateu.parse_date(date0)

    if date1 is False: date1=datetime.datetime(9999,1,1)
    else: date1=dateu.parse_date(date1)

    out=None
    for f in self.files:
      time=read_time(f)
      print(time[0],time[-1])
      if not np.any((time>=date0)&(time<=date1)): continue

      if not quiet: print('-> extracting file %s'%f)
      res=wrf_file_data(f,quiet=quiet)
      time=res['time']
      i,=np.where(time>=date0)
      j,=np.where(time<=date1)
      i=i[0]
      j=j[-1]

      for k in res.keys():
        if k =='time':
          res[k]=res[k][i:j+1]
        else:
          res[k].data=res[k].data[i:j+1]

      if not out:
        out=res
      else:
        i0=np.where(res['time']>out['time'][-1])[0][0]
        for k in res.keys():
          if k =='time':
            out[k]=np.hstack((out[k],res[k][i0:]))
          else:
            out[k].data=np.concatenate((out[k].data,res[k].data[i0:]))

    return out


def parse_time(t):
  '''Parse wrf time as 'y','y','y','y','-','... ect or as yyyymmdd.xx
  '''
  from dateutil import parser

  if t.dtype.char=='S' and t.ndim==1: # one time of the type 'y','y','y','y','-','... ect
    t.shape=1,t.size

  nt=t.shape[0]
  tout=np.zeros(nt,datetime.datetime)

  if t.dim==2: # 'y','y','y','y','-','... ect
    for i in range(nt):
       tmp=''.join(t[i])
       tout[i]=parser.parse(tmp.replace('_',' '))
  else: # yyyymmdd.xx
    for i in range(nt):
      L=t[i]
      ano,mes,dia,addDay=int(str(L)[:4]),int(str(L)[4:6]),int(str(L)[6:8]),float(str(L)[8:])
      tout[i]=datetime.datetime(ano,mes,dia)+datetime.timedelta(days=addDay)

  return tout


def accum2avg(v,dt): # accum from initial reccord!
  
  nt,nj,ni=v.shape
  v=(v[1:]-v[:-1])/dt.seconds # Y s-1
  vc=(v[1:]+v[:-1])/2.
  V=np.zeros((nt,nj,ni),v.dtype)
  V[1:-1]=vc
  V[0]=v[0] # --> closest one (t diff is dt/2)
  V[-1]=v[-1]
  return V

#  dt=(t[-1]-t[0])
#  ndays=dt.days+dt.seconds/86400.
#  ndays=np.ceil(ndays)
#  Day0=datetime.datetime(t[0].year,t[0].month,t[0].day)
#  for n in range(ndays):
#    day0=Day0+datetime.datetime(


def read_time(file):
  return parse_time(netcdf.use(file,'Times'))
 

def wrf_file_data(file,quiet=False):
  '''
  WRF data for ROMS

  '''

  out={}

  # time:
  if not quiet: print(' --> get time')
  time=read_time(file)

  out['time']=time

  # lon,lat:
  if not quiet: print(' --> reading x,y')
  x=netcdf.use(file,'XLONG',Time=0)#**{'0': 0})
  y=netcdf.use(file,'XLAT',Time=0)#**{'0': 0})

  # tair [K-->C]
  if not quiet: print(' --> T air')
  tair=netcdf.use(file,'T2')-273.15
  out['tair']=Data(x,y,tair,'Celsius')

  # R humidity [kg/kg --> 0--1]
  if not quiet: print(' --> R humidity from QV at 2m')
  wv=netcdf.use(file,'Q2') # water vapor mixing ratio at 2m
  rhum=wv/air_sea.qsat(tair)
  rhum[rhum>1]=1
  out['rhum']=Data(x,y,rhum,'0--1')

  # surface pressure [Pa]
  if not quiet: print(' --> Surface pressure')
  pres=netcdf.use(file,'PSFC')
  out['pres']=Data(x,y,pres,'Pa')

  # P rate [mm --> cm day-1]
  if not quiet: print(' --> P rate (rainc+rainnc)')
  rainc  = netcdf.use(file,'RAINC')
  rainnc = netcdf.use(file,'RAINNC')
  prate=rainc+rainnc
  if not quiet: print('      accum2avg...')
  prate=accum2avg(prate,dt=time[1]-time[0]) # mm s-1
  conv= 0.1*86400       # from mm s-1      --> cm day-1
  prate=prate*conv # cm day-1
  prate[prate<0]=0 # interpolation errors may result in negative rain!
  out['prate']=Data(x,y,prate,'cm day-1')

  # LW, SW, latent, sensible signs:
  # positive (downward flux, heating) or negative (upward flux, cooling)
  #https://www.myroms.org/forum/viewtopic.php?f=1&t=2621

  # Net shortwave flux  [W m-2]
  if not quiet: print(' --> Net shortwave flux')
  sw_down=netcdf.use(file,'SWDOWN')
  albedo=netcdf.use(file,'ALBEDO')
  sw_net=sw_down*(1-albedo)
  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')

  # Net longwave flux  [W m-2]
  if not quiet: print(' --> Net longwave flux')
  lw_down=netcdf.use(file,'GLW') # positive
  # sst needed:
  if not quiet: print('     --> SST for LW up')
  sst=netcdf.use(file,'SST') # K
  lw_net = air_sea.lwhf(sst,lw_down) # positive down
  # here vars have roms-agrif signs --> radlw is positive upward!
  #conversion to ROMS is done in surface.py
  out['radlw']=Data(x,y,-lw_net,'W m-2',info='positive upward')
  out['dlwrf']=Data(x,y,-lw_down,'W m-2',info='positive upward')

  # U and V wind speed 10m
  if not quiet: print(' --> U and V wind')
  uwnd=netcdf.use(file,'U10')
  vwnd=netcdf.use(file,'V10')
  if not quiet: print(' --> calc wind speed and stress')
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=Data(x,y,speed,'m s-1')
  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  out['sustr']=Data(x,y,taux,'Pa')
  out['svstr']=Data(x,y,tauy,'Pa')

  # Cloud cover [0--1]:
  if not quiet: print(' --> Cloud cover for LONGWAVE. Use LONGWAVE_OUT instead...')
  if 0:
    pass
# next code is wrong! If cloud cover is really needed, it needs to be calculated using wrfpost.
# See http://www2.mmm.ucar.edu/wrf/users/docs/user_guide/users_guide_chap8.html#_ARWpost_1
#
#  if 'CLDFRA' in netcdf.varnames(file):
#    clouds=netcdf.use(file,'CLDFRA').sum(-3)
#    clouds=np.where(clouds>1,1,clouds)
  else:
    if not quiet: print('CLDFRA not found!! Using SST and air_sea.cloud_fraction')
    sst=netcdf.use(file,'SST')-273.15
    clouds=air_sea.cloud_fraction(lw_net,sst,tair,rhum,Wtype='net')
    clouds[clouds<0]=0
    clouds[clouds>1]=1

  out['cloud']=Data(x,y,clouds,'fraction (0--1)')

  return out

