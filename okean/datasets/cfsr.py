'''
Tools to extract the fields required for ocean model atmospheric
bulk forcing, from the atm dataset CFSR

Martinho MA, UFBA, Sep 2013
'''

import os
import datetime
import numpy as np
from okean import netcdf, dateu, air_sea
import glob


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


class CFSRData:
  '''
  CFSR data extraction
  '''
  def __init__(self,basefolder):
    '''
    Input:
      basefolder: folder with all data for a common period. Data
      should be in folders as:
        - basefolder/ST   Surface temperature
        - basefolder/RH   Relative humidity
        - basefolder/SP   Surface pressure
        - basefolder/PR   Precipitation rate
        - basefolder/RAD  Short&long wave radiation
        - basefolder/UV   Winf u, v 10 m
        - basefolder/CC   Cloud cover


    See cfsr_file_data for more info
    '''

    self.basefolder=basefolder

    files={}
    files['st']  = os.path.join(basefolder,'ST/*.nc')
    files['rh']  = os.path.join(basefolder,'RH/*.nc')
    files['sp']  = os.path.join(basefolder,'SP/*.nc')
    files['pr']  = os.path.join(basefolder,'PR/*.nc')
    files['rad'] = os.path.join(basefolder,'RAD/*.nc')
    files['uv']  = os.path.join(basefolder,'UV/*.nc')
    files['cc']  = os.path.join(basefolder,'CC/*.nc')
    self.files=files

  def data(self,date0=False,date1=False,quiet=True):
    '''
    Returns atm data form all times in basefolder files or between
    date0 (>=) and date1 (<=)
    '''
    # get data for all times in file:
    res=cfsr_file_data(self.files,quiet=quiet)
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


def cfsr_file_data(files,quiet=False):
  '''
  Returns bulk data from one CFRS files
  '''


  def load_time(f):
    time=np.array((),datetime.datetime)
    ff=glob.glob(f)
    ff.sort()
    for f in ff: time=np.append(time,netcdf.nctime(f,'time'))
    return time


  def fix_time(t,var,t0,t1):
    # convert 1h, 7h, ... to 0h, 6h, ...
    if t[0].hour in [1,7,13,19]: # not all! sp analysis starts at 0, 6,...!
      print '     1,7,... to 0,6,...'
      var=(var[1:]*5+var[:-1]*1)/6.
      t=t[1:]-datetime.timedelta(hours=1)

    cond=(t>=t0)&(t<=t1)
    t=t[cond]
    var=var[cond]

    if t[0]>t0:
      dt=t[0]-t0
      dt=dt.days*24+dt.seconds/3600. # hours
      print 'missind data at start: %.2d h missing --> repeating 1st data'%dt
      v=np.zeros((var.shape[0]+1,)+var.shape[1:],var.dtype)
      v[1:]=var
      v[0]=var[0]
      var=v

    if t[-1]<t1:
      dt=t1-t[-1]
      dt=dt.days*24+dt.seconds/3600. # hours
      print 'missind data at end: %.2d h missing --> repeating last data'%dt
      v=np.zeros((var.shape[0]+1,)+var.shape[1:],var.dtype)
      v[:-1]=var
      v[-1]=var[-1]
      var=v

    return var


  out={}

  # time:
  if 0:
    time=netcdf.nctime(files['cc'],'time')
    # files have diff units !! so, cannot load all times at once!
    # thse result will use only units of 1st file!!
  else:
    time=load_time(files['cc'])


  out['time']=time

  # T air [K->C]
  if not quiet: print ' --> T air'
  f=files['st']
  tair=netcdf.use(f,'TMP_L103')
  tair=tair-273.15
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    tair=fix_time(ttmp,tair,time[0],time[-1])
    print ' ...fixed!'
  out['tair']=Data(x,y,tair,'C')


  # R humidity [%-->0--1]
  if not quiet: print ' --> R humidity'
  f=files['rh']
  rhum=netcdf.use(f,'R_H_L103')
  rhum=rhum/100.
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    rhum=fix_time(ttmp,rhum,time[0],time[-1])
    print ' ...fixed!'
  out['rhum']=Data(x,y,rhum,'0--1')


  # surface pressure [Pa]
  if not quiet: print ' --> Surface pressure'
  f=files['sp']
  pres=netcdf.use(f,'PRES_L1')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    pres=fix_time(ttmp,pres,time[0],time[-1])
    print ' ...fixed!'
  out['pres']=Data(x,y,pres,'Pa')


  # P rate [kg m-2 s-1 -> cm/d]
  if not quiet: print ' --> P rate'
  f=files['pr']
  prate=netcdf.use(f,'PRATE_L1')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  # Conversion kg m^-2 s^-1  to cm/day
  prate=prate*86400*100/1000.
  prate=np.where(abs(prate)<1.e-4,0,prate)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    prate=fix_time(ttmp,prate,time[0],time[-1])
    print ' ...fixed!'
  out['prate']=Data(x,y,prate,'cm/d')


  # Net shortwave flux  [W/m^2]
  if not quiet: print ' --> Net shortwave flux'
  if not quiet: print '       SW down'
  f=files['rad']
  sw_down=netcdf.use(f,'DSWRF_L1_Avg_1')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  if not quiet: print '       SW up'
  sw_up=netcdf.use(f,'USWRF_L1_Avg_1')
  sw_net=sw_down-sw_up
  sw_net=np.where(sw_net<1.e-10,0,sw_net)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    sw_net=fix_time(ttmp,sw_net,time[0],time[-1])
    print ' ...fixed!'
  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')


  # Net longwave flux  [W/m^2]
  if not quiet: print ' --> Net longwave flux'
  if not quiet: print '       LW down'
  f=files['rad']
  lw_down=netcdf.use(f,'DLWRF_L1_Avg_1')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  if not quiet: print '       LW up'
  lw_up=netcdf.use(f,'ULWRF_L1_Avg_1')
  lw_net=lw_down-lw_up
  lw_net=np.where(np.abs(lw_net)<1.e-10,0,lw_net)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    lw_net=fix_time(ttmp,lw_net,time[0],time[-1])
    lw_down=fix_time(ttmp,lw_down,time[0],time[-1])
    print ' ...fixed!'
  # ROMS convention: positive upward
  out['radlw']=Data(x,y,-lw_net,'W m-2',info='positive upward')
  # downward lw:
  out['dlwrf']=Data(x,y,-lw_down,'W m-2',info='negative... downward')


  # U and V wind speed 10m
  if not quiet: print ' --> U and V wind'
  f=files['uv']
  uwnd=netcdf.use(f,'U_GRD_L103')
  vwnd=netcdf.use(f,'V_GRD_L103')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    uwnd=fix_time(ttmp,uwnd,time[0],time[-1])
    vwnd=fix_time(ttmp,vwnd,time[0],time[-1])
    print ' ...fixed!'
  #
  if not quiet: print ' --> calc wind speed and stress'
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=Data(x,y,speed,'m s-1')
  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  out['sustr']=Data(x,y,taux,'Pa')
  out['svstr']=Data(x,y,tauy,'Pa')


  # Cloud cover [0--100 --> 0--1]:
  if not quiet: print ' --> Cloud cover'
  f=files['cc']
  clouds=netcdf.use(f,'T_CDC_L200')
  x=netcdf.use(f,'lon'); x[x>180]=x[x>180]-360
  y=netcdf.use(f,'lat')
  x,y=np.meshgrid(x,y)
  clouds=clouds/100.
  # check time:
  ttmp=load_time(f)
  if ttmp.size==time.size and np.all(ttmp==time): print '    time ok'
  else:
    print '   time differs !!!!',
    clouds=fix_time(ttmp,clouds,time[0],time[-1]) # cc is currently used to get
                                                  # time, so this step is not needed!
    print ' ...fixed!'
  out['cloud']=Data(x,y,clouds,'fraction (0--1)')

  # rhum has different resolution (0.5, just like dew point!)
  # so, i can edit surface.py or just interpolate here rhum to
  # other vars resolution:
  if out['rhum'].data.shape!=out['uwnd'].data.shape:
    from okean import calc
    print 'rhum shape differs!! --> interp:'
    nt,ny,nx=out['uwnd'].data.shape
    x,y=out['uwnd'].x,out['uwnd'].y
    rhum=np.zeros((nt,ny,nx), out['rhum'].data.dtype)
    for it in range(nt):
      if it%100==0: print '  %d of %d'%(it,nt)
      rhum[it]=calc.griddata(out['rhum'].x,out['rhum'].y,out['rhum'].data[it],x,y)

    out['rhum']=Data(x,y,rhum,'0--1')


  return out
