import pylab as pl
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


class ERA5Data:
  '''
  ECMWF ERA5 data extraction
  '''
  def __init__(self,basefolder):
    '''
    Input:
      basefolder: folder with the data files(s)
    '''
    self.basefolder=basefolder

    import glob
    import os
    self.files=glob.glob(os.path.join(self.basefolder,'*.nc'))


  def data(self,date0=False,date1=False,quiet=True):
    '''
    Returns atm data form all times in basefolder files or between
    date0 (>=) and date1 (<=)
    '''
    # get data for all times in file:
    res=era5_file_data(self.files,quiet=quiet)
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
      for k in res.keys():
        if k =='time':
          res[k]=res[k][i:j+1]
        else:
          res[k].data=res[k].data[i:j+1,...]

    return res


#def accum2avg(v,nforec):
#  '''
#  Accumulation fields to average
#  See http://www.ecmwf.int/products/data/archive/data_faq.html#accumulations
#3  for details
#  '''
#  N_DAILY_SIMS  = 2 # 00h and 12 h
#  N_FOREC_STEPS = nforec # 2 or 4 (3,6,9,12 or 6,12)
#  DT_FOREC=24/(N_DAILY_SIMS*N_FOREC_STEPS)
#
#  # accum to avg (intermediate times):
#  v.shape=v.shape[0]/N_FOREC_STEPS,N_FOREC_STEPS,v.shape[1],v.shape[-1]
#  v[:,1:,...]=(v[:,1:,...]-v[:,:-1,...])/(3600*DT_FOREC)
#  v[:,0,...]=v[:,0,...]/(3600*DT_FOREC)
#
#  # data at original time, except last one:
#  v[:,:-1,...]=(v[:,1:,...]+v[:,:-1,...])/2.
#
#  # data at last forec time is still unknown, ie,
#  # v[:,-1,...] is currently wrong!
#  # average with next first elements:
#  v[:-1,-1,...]=(v[:-1,-1,...]+v[1:,0,...])/2.
#
#  # return to original shape:
#  v.shape=v.shape[0]*v.shape[1],v.shape[2],v.shape[-1]
#
#  # now, only last value is wrong !
#  # let us consider is the same as previous
#  v[-1,...]=v[-2,...]
#
#  return v


#def fill_t0(v,fill='next'):
#  '''
#  Forecast only variables don't have t=00h (for the first day only!),
#  so, fill that time with data from next time, usually 3h, fill='next',
#  or fill with one value, fill=0, for instance
#  This should be done with forecast variables (like radiation SW, LW
#  and precipitation)
#  The output will have forst dimension increase by 1
#  '''
#
#  if fill=='next': fill=v[0,...]
#
#  sh=list(v.shape)
#  sh[0]+=1
#  vv=np.zeros(sh,dtype=v.dtype)
#  vv[0,...]=fill
#  vv[1:,...]=v
#  return vv


def fill_tend(v,fill='prev'):
  '''
  Downloaded data may not include last requested date, eg, 24h of last day.
  So we can fill that time
  with data from previous time (fill='prev') or with some value, fill=20,
  for instance (not a good idea, better use previous value!!)
  '''
  if fill=='prev': fill=v[-1]
  sh=list(v.shape)
  sh[0]+=1
  vv=np.zeros(sh,dtype=v.dtype)
  vv[-1]=fill
  vv[:-1]=v
  return vv


def era5_file_data(files,quiet=False):
  '''
  ECMWF ERA5 data for ROMS

  variables:
    # radiation:
      msdwlwrf  -  mean_surface_downward_long_wave_radiation_flux
      #msdwswrf -  mean_surface_downward_short_wave_radiation_flux -- not needed
      msnlwrf   -  mean_surface_net_long_wave_radiation_flux
      msnswrf   -  mean_surface_net_short_wave_radiation_flux
    # rain:
      mtpr      - mean_total_precipitation_rate
    # wind:
      u10       - 10m_u_component_of_wind
      v10       - 10m_v_component_of_wind
    # temp:
      t2m       - 2m_temperature
      d2m       - 2m_dewpoint_temperature (for relative humidity)
    # pres:
      sp        - surface_pressure
    # clouds:
      tcc       - 'total_cloud_cover
  '''

###  # some variables may have different names! 
###  Vars={}
###  Vars['v10u']='v10u','u10'
###  Vars['v10v']='v10v','v10'
###  Vars['v2t']='v2t','t2m'
###  Vars['v2d']='v2d','d2m'
###
####  def find_v(name):
###    if name in Vars.keys():
###      for v in Vars[name]:
###        if varfile(v): return v
###    else: return name


  def varfile(var):
    for f in files:
      if var in netcdf.varnames(f): return f


##3  def check_var_type(var):
###    # new interim dataserver provides forec+analysis vars with extra dim
###    # 'type', 0 or 1
###    if var.ndim==4:
###      if not quiet: print('      dealing with var type... '),
###      v=np.zeros(var.shape[1:],var.dtype)
###      v[::2]=var[0,::2,...]
###      v[1::2]=var[1,1::2,...]
###      var=v
###      if not quiet: print('done.')
###
###    return var


  out={}

  # time:
  File=varfile('t2m') # air temp, for instance
  if not quiet: print(' reading  time from file %s' % File)
  time=netcdf.nctime(File,'time')
  # check if last time at 00h:
  if  time[-1].hour!=0:
    dateEnd=datetime.datetime(time[-1].year,time[-1].month,time[-1].day)+datetime.timedelta(days=1)
    time=np.append(time,dateEnd)
  else: dateEnd=0
  out['time']=time

###  # all times from analysis file, except last ind which will be
######  # the last time of forecast file
###  aFile=varfile(find_v('v2t')) # air temp, for instance
######  fFile=varfile(find_v('ssr')) # sw rad, for instance
###  if not quiet: print(' reading "analysis" time from file %s' % aFile)
###  aTime=netcdf.nctime(aFile,'time')
###  aTime.sort() # analysis+forecast files may not have time sorted!!
###  if not quiet: print(' reading "forecast" time from file %s' % fFile)
###  fTime=netcdf.nctime(fFile,'time')
###  fTime.sort() # this one should be sorted...
###  time=np.append(aTime,fTime[-1])
###  out['time']=time
##
#  # calc number of forecast steps stored,nforec (used by accum2avg)
#  if [fTime[i].hour for i in range(8)]==range(3,22,3)+[0]: nforec=4
#  elif [fTime[i].hour for i in range(4)]==range(6,19,6)+[0]: nforec=2
#  else:
#    if not quiet: print('INTERIM WRONG TIME: cannot n forec steps')
#    return
#
###  if not quiet: print(' ==> n forecast steps = %d' % nforec)

  # x,y:
  if not quiet: print(' reading x,y from file %s' % File)
  x=netcdf.use(File,'longitude')
  y=netcdf.use(File,'latitude')
  x[x>180]=x[x>180]-360
  if x.ndim==1 and y.ndim==1:
    x,y=np.meshgrid(x,y)

  # tair [K-->C]
  if not quiet: print(' --> T air')
  vname='t2m'
  f=varfile(vname)
  tair=netcdf.use(f,vname)-273.15
  if dateEnd:
    if not quiet: print('      fill_tend...')
    tair=fill_tend(tair)

  out['tair']=Data(x,y,tair,'Celsius')


  # R humidity [0--1]
  if not quiet: print(' --> R humidity (from T dew)')
  vname='d2m'
  f=varfile(vname)
  Td=netcdf.use(f,vname)-273.15
  if dateEnd:
    if not quiet: print('      fill_tend... (T dew)')
    Td=fill_tend(Td)

  T=tair
  rhum=air_sea.relative_humidity(T,Td)
  rhum[rhum>1]=1
  out['rhum']=Data(x,y,rhum,'0--1')


  # surface pressure [Pa]
  if not quiet: print(' --> Surface pressure')
  vname='sp'
  f=varfile(vname)
  pres=netcdf.use(f,vname)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    pres=fill_tend(pres)

  out['pres']=Data(x,y,pres,'Pa')


  def avg_fix_time(v,DT):
    '''Fix data to right time in avg rate fields (ie, prev half hour to now)
       See:
       https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Meanratesandaccumulations
    '''
    DTstep=1
    a=DTstep/2.
    b=DT-a
    u=np.zeros_like(v)
    u[:-1]=(v[:-1]*b+v[1:]*a)/(a+b)
    # last one lost, use prev value (from 30 min before)
    u[-1]=v[-1]
    return u


  DT=(time[1]-time[0]).total_seconds()/3600. # hours

  # P rate [kg m-2 s-1 --> cm day-1]
  if not quiet: print(' --> P rate')
  vname='mtpr'
  f=varfile(vname)
  prate=netcdf.use(f,vname)
  if not quiet: print('      avg_fix_time - DTstep=1h - DT=%.2f h'%DT)
  prate=avg_fix_time(prate,DT)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    prate=fill_tend(prate)

  conv= 100*86400/1000. # from kg m-2 s-1 --> cm day-1
  prate=prate*conv # cm day-1
  prate[prate<0]=0
  out['prate']=Data(x,y,prate,'cm day-1')

#  dt=(TIMe[1]-time[0]).total_seconds()/3600.
#  if not quiet: print('      accum to avg at correct time - DTstep=1h - DT=%.2f h'%DT)
#  prate2=accum2avg(prate,DT)
#  prate2=prate2*1000
#  return prate,prate2
#  prate=check_var_type(prate)
#  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
#  if not quiet: print('      accum2avg...')
#  prate=accum2avg(prate,nforec)
#  conv= 100*86400       # from m s-1      --> cm day-1
#  #conv= 100*86400/1000. # from kg m-2 s-1 --> cm day-1
#  prate=prate*conv # cm day-1
#  if not quiet: print('      fill_t0...')
#  prate=fill_t0(prate)
#  prate[prate<0]=0
#  out['prate']=Data(x,y,prate,'cm day-1')
#
#  return out

  # Net shortwave flux  [W m-2 --> W m-2]
  if not quiet: print(' --> Net shortwave flux')
  vname='msnswrf'
  f=varfile(vname)
  sw_net=netcdf.use(f,vname)
  if not quiet: print('      avg_fix_time - DTstep=1h - DT=%.2f h'%DT)
  sw_net=avg_fix_time(sw_net,DT)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    sw_net=fill_tend(sw_net)

  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')


  # Net longwave flux  [W m-2 --> W m-2]
  if not quiet: print(' --> Net longwave flux')
  vname='msnlwrf'
  f=varfile(vname)
  lw_net=netcdf.use(f,vname)*-1 # positive upward (*-1)
  # here vars have roms-agrif signs --> radlw is positive upward!
  # conversion to ROMS is done in surface.py
  if not quiet: print('      avg_fix_time - DTstep=1h - DT=%.2f h'%DT)
  lw_net=avg_fix_time(lw_net,DT)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    lw_net=fill_tend(lw_net)

  out['radlw']=Data(x,y,lw_net,'W m-2',info='positive upward')


  # longwave down:
  if not quiet: print(' --> Down longwave flux')
  vname='msdwlwrf'
  f=varfile(vname)
  lw_down=netcdf.use(f,vname)*-1 # positive upward (*-1)
  if not quiet: print('      avg_fix_time - DTstep=1h - DT=%.2f h'%DT)
  lw_down=avg_fix_time(lw_down,DT)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    lw_down=fill_tend(lw_down)

  out['dlwrf']=Data(x,y,lw_down,'W m-2',info='positive upward')


  # U and V wind speed 10m
  if not quiet: print(' --> U and V wind')
  vname='u10'
  f=varfile(vname)
  uwnd=netcdf.use(f,vname)
  vname='v10'
  f=varfile(vname)
  vwnd=netcdf.use(f,vname)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    uwnd=fill_tend(uwnd)
    vwnd=fill_tend(vwnd)

  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  # speed and stress:
  if 0:
    if not quiet: print(' --> calc wind speed and stress')
    speed = np.sqrt(uwnd**2+vwnd**2)
    taux,tauy=air_sea.wind_stress(uwnd,vwnd)
    out['wspd']=Data(x,y,speed,'m s-1')
    out['sustr']=Data(x,y,taux,'Pa')
    out['svstr']=Data(x,y,tauy,'Pa')


  # Cloud cover [0--1]:
  if not quiet: print(' --> Cloud cover')
  vname='tcc'
  f=varfile(vname)
  clouds=netcdf.use(f,vname)
  if dateEnd:
    if not quiet: print('      fill_tend...')
    clouds=fill_tend(clouds)

  out['cloud']=Data(x,y,clouds,'fraction (0--1)')

  return out

