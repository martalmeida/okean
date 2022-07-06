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


class INTERIMData:
  '''
  ECMWF ERA INTERIM data extraction
  '''
  def __init__(self,basefolder):
    '''
    Input:
      basefolder: folder with all the analysis+forecast files
        for a common period !! For instance .../ecmwf_data/2008

    The number of forecast steps  in forecast only data is
    usually 4 (forecast step 3,6,9 and 12). If using all
    these steps, the analysis data must include the forecast
    steps 3 and 9, so that we end up with data every 3 hours.

    n forecast  2 means analysis data includes only analysis (step 0)
    and forecast data inlcudes the steps 6 and 12 only

    See interim_file_data for more info
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
    res=interim_file_data(self.files,quiet=quiet)
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


def accum2avg(v,nforec):
  '''
  Accumulation fields to average
  See http://www.ecmwf.int/products/data/archive/data_faq.html#accumulations
  for details
  '''
  N_DAILY_SIMS  = 2 # 00h and 12 h
  N_FOREC_STEPS = nforec # 2 or 4 (3,6,9,12 or 6,12)
  DT_FOREC=24/(N_DAILY_SIMS*N_FOREC_STEPS)

  # accum to avg (intermediate times):
  v.shape=v.shape[0]//N_FOREC_STEPS,N_FOREC_STEPS,v.shape[1],v.shape[-1]
  v[:,1:,...]=(v[:,1:,...]-v[:,:-1,...])/(3600*DT_FOREC)
  v[:,0,...]=v[:,0,...]/(3600*DT_FOREC)

  # data at original time, except last one:
  v[:,:-1,...]=(v[:,1:,...]+v[:,:-1,...])/2.

  # data at last forec time is still unknown, ie,
  # v[:,-1,...] is currently wrong!
  # average with next first elements:
  v[:-1,-1,...]=(v[:-1,-1,...]+v[1:,0,...])/2.

  # return to original shape:
  v.shape=v.shape[0]*v.shape[1],v.shape[2],v.shape[-1]

  # now, only last value is wrong !
  # let us consider is the same as previous
  v[-1,...]=v[-2,...]

  return v


def fill_t0(v,fill='next'):
  '''
  Forecast only variables don't have t=00h (for the first day only!),
  so, fill that time with data from next time, usually 3h, fill='next',
  or fill with one value, fill=0, for instance
  This should be done with forecast variables (like radiation SW, LW
  and precipitation)
  The output will have forst dimension increase by 1
  '''

  if fill=='next': fill=v[0,...]

  sh=list(v.shape)
  sh[0]+=1
  vv=np.zeros(sh,dtype=v.dtype)
  vv[0,...]=fill
  vv[1:,...]=v
  return vv


def fill_tend(v,fill='prev'):
  '''
  Analysis (maybe with forecast data at 3 and 9h) vars should not include
  the data for last 24h (last day only). So we can fill that time
  with data from previous time (fill='prev') or with some value, fill=20,
  for instance (not a good idea, better use previous value!!)
  '''
  if fill=='prev': fill=v[-1,...]
  sh=list(v.shape)
  sh[0]+=1
  vv=np.zeros(sh,dtype=v.dtype)
  vv[-1,...]=fill
  vv[:-1,...]=v
  return vv


def interim_file_data(files,quiet=False):
  '''
  ECMWF ERA INTERIM data for ROMS

  To be used with data obtained from the new server:
  http://apps.ecmwf.int/datasets/

  and not the old one (http://data-portal.ecmwf.int/data/d/interim_daily/)
  To deal with data from old server use module interim_past

  => forecast vars:
  time=00h, 12h
  step=3,6,9,12 --> n forec steps=4

      needed (suggestion):
      - Surface net solar radiation (ssr)
      - Surface thermal radiation (str)
      - Total precipitation (tp)
      others:
      - Surface thermal radiation downwards (strd)
      - Evaporation (e)
      - ...

  =>analysis+forec vars: (interim analysis starts every 6h; forecast starts every 12h !!)
  time=00h,6h,12h,18h
  step=0,3,9 (3 and 9 are for the forecast 00h and 12h)

      needed (suggestion):
      - Surface pressure (sp)
      - Total cloud cover    (tcc)
      - 10 metre U wind component (v10u or u10)
      - 10 metre V wind component  (v10v or v10)
      - 2 metre temperature (v2t or t2m)
      - 2 metre dewpoint temperature (v2d or d2m)

  Accumulated vars (rad SW, LW and precipitation are converted to averages
  by acum2avg
  '''

  # some variables may have different names! 
  Vars={}
  Vars['v10u']='v10u','u10'
  Vars['v10v']='v10v','v10'
  Vars['v2t']='v2t','t2m'
  Vars['v2d']='v2d','d2m'

  def find_v(name):
    if name in Vars.keys():
      for v in Vars[name]:
        if varfile(v): return v
    else: return name


  def varfile(var):
    for f in files:
      if var in netcdf.varnames(f): return f


  def check_var_type(var):
    # new interim dataserver provides forec+analysis vars with extra dim
    # 'type', 0 or 1
    if var.ndim==4:
      if not quiet: print('      dealing with var type... '),
      v=np.zeros(var.shape[1:],var.dtype)
      v[::2]=var[0,::2,...]
      v[1::2]=var[1,1::2,...]
      var=v
      if not quiet: print('done.')

    return var


  out={}

  # time:
  # all times from analysis file, except last ind which will be
  # the last time of forecast file
  aFile=varfile(find_v('v2t')) # air temp, for instance
  fFile=varfile(find_v('ssr')) # sw rad, for instance
  if not quiet: print(' reading "analysis" time from file %s' % aFile)
  aTime=netcdf.nctime(aFile,'time')
  aTime.sort() # analysis+forecast files may not have time sorted!!
  if not quiet: print(' reading "forecast" time from file %s' % fFile)
  fTime=netcdf.nctime(fFile,'time')
  fTime.sort() # this one should be sorted...
  time=np.append(aTime,fTime[-1])
  out['time']=time

  # calc number of forecast steps stored,nforec (used by accum2avg)
  if [fTime[i].hour for i in range(8)]==list(range(3,22,3))+[0]: nforec=4
  elif [fTime[i].hour for i in range(4)]==list(range(6,19,6))+[0]: nforec=2
  else:
    if not quiet: print('INTERIM WRONG TIME: cannot n forec steps')
    return

  if not quiet: print(' ==> n forecast steps = %d' % nforec)

  # x,y:
  if not quiet: print(' reading x,y from file %s' % files[0])
  x=netcdf.use(files[0],'longitude')
  y=netcdf.use(files[0],'latitude')
  x[x>180]=x[x>180]-360
  if x.ndim==1 and y.ndim==1:
    x,y=np.meshgrid(x,y)


  # tair [K-->C]
  if not quiet: print(' --> T air')
  vname=find_v('v2t')
  f=varfile(vname)
  # time may not be monotonically increasing !!
  # when using mix of analysis and forecast variables and steps
  sortInds=np.argsort(netcdf.use(f,'time'))
  tair=netcdf.use(f,vname,time=sortInds)-273.15
  tair=check_var_type(tair)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend...')
  tair=fill_tend(tair)
  out['tair']=Data(x,y,tair,'Celsius')

  # R humidity [0--1]
  if not quiet: print(' --> R humidity (from T dew)')
  vname=find_v('v2d')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  Td=netcdf.use(f,vname,time=sortInds)-273.15
  Td=check_var_type(Td)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend... (T dew)')
  Td=fill_tend(Td)
  T=tair
  rhum=air_sea.relative_humidity(T,Td)
  ##  rhum=((112-0.1*T+Td)/(112+0.9*T))**8
  rhum[rhum>1]=1
  out['rhum']=Data(x,y,rhum,'0--1')

  # surface pressure [Pa]
  if not quiet: print(' --> Surface pressure')
  vname=find_v('sp')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  pres=netcdf.use(f,vname,time=sortInds)
  pres=check_var_type(pres)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend...')
  pres=fill_tend(pres)
  out['pres']=Data(x,y,pres,'Pa')

  # P rate [m --> cm day-1]
  if not quiet: print(' --> P rate')
  vname=find_v('tp')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  prate=netcdf.use(f,vname,time=sortInds)
  prate=check_var_type(prate)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      accum2avg...')
  prate=accum2avg(prate,nforec)
  conv= 100*86400       # from m s-1      --> cm day-1
  #conv= 100*86400/1000. # from kg m-2 s-1 --> cm day-1
  prate=prate*conv # cm day-1
  if not quiet: print('      fill_t0...')
  prate=fill_t0(prate)
  prate[prate<0]=0
  out['prate']=Data(x,y,prate,'cm day-1')

  # Net shortwave flux  [W m-2 s+1 --> W m-2]
  if not quiet: print(' --> Net shortwave flux')
  vname=find_v('ssr')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  sw_net=netcdf.use(f,vname,time=sortInds)
  sw_net=check_var_type(sw_net)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      accum2avg...')
  sw_net=accum2avg(sw_net,nforec)
  if not quiet: print('      fill_t0...')
  sw_net=fill_t0(sw_net)
  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')


  # Net longwave flux  [W m-2 s+1 --> W m-2]
  if not quiet: print(' --> Net longwave flux')
  vname=find_v('str')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  lw_net=netcdf.use(f,vname,time=sortInds)*-1 # let us consider positive upward (*-1)
  lw_net=check_var_type(lw_net)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      accum2avg...')
  lw_net=accum2avg(lw_net,nforec)
  if not quiet: print('      fill_t0...')
  lw_net=fill_t0(lw_net)
  out['radlw']=Data(x,y,lw_net,'W m-2',info='positive upward')

  # longwave down:
  # can be obtained from clouds!!
  if not quiet: print(' --> Down longwave flux')
  vname=find_v('strd')
  f=varfile(vname)
  if f:
    sortInds=np.argsort(netcdf.use(f,'time'))
    lw_down=netcdf.use(f,vname,time=sortInds)*-1 # let us consider positive upward (*-1)
    lw_down=check_var_type(lw_down)
    if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
    if not quiet: print('      accum2avg...')
    lw_down=accum2avg(lw_down,nforec)
    if not quiet: print('      fill_t0...')
    lw_down=fill_t0(lw_down)
    out['dlwrf']=Data(x,y,lw_down,'W m-2',info='negative... downward')
  else:  print('down long wave CANNOT BE USED')

  # U and V wind speed 10m
  if not quiet: print(' --> U and V wind')
  vname=find_v('v10u')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  uwnd=netcdf.use(f,vname,time=sortInds)
  uwnd=check_var_type(uwnd)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend...')
  uwnd=fill_tend(uwnd)
  vname=find_v('v10v')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  vwnd=netcdf.use(f,vname,time=sortInds)
  vwnd=check_var_type(vwnd)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend...')
  vwnd=fill_tend(vwnd)

  if not quiet: print(' --> calc wind speed and stress')
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=Data(x,y,speed,'m s-1')
  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  out['sustr']=Data(x,y,taux,'Pa')
  out['svstr']=Data(x,y,tauy,'Pa')

  # Cloud cover [0--1]:
  if not quiet: print(' --> Cloud cover')
  vname=find_v('tcc')
  f=varfile(vname)
  sortInds=np.argsort(netcdf.use(f,'time'))
  clouds=netcdf.use(f,vname,time=sortInds)
  clouds=check_var_type(clouds)
  if not quiet and np.any(sortInds!=range(len(sortInds))): print('      sort DONE')
  if not quiet: print('      fill_tend...')
  clouds=fill_tend(clouds)
  out['cloud']=Data(x,y,clouds,'fraction (0--1)')

  return out

