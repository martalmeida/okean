import numpy as np
import datetime

from okean import calc, cookbook as cb, dateu as dts, netcdf, air_sea
from okean.roms import roms, roms_tools as rt
import gennc


def load_blkdata_gfs(gfspath,date0,date1=False,nforec=0,quiet=True):
  '''
  nforec=0 means not forecast, auto means default from config file
  '''
  from okean.datasets import gfs
  a=gfs.GFSData(gfspath)
  if nforec==0:
    data,miss=a.data_analysis(date0,date1,quiet=quiet)
  else:
    data,miss=a.data_forecast(date0,nforec,quiet=quiet)

  out=cb.odict()
  time=data.keys()
  for t in time:
    out[t]={}
    for k in data[t].keys():

      if  k.startswith('INFO'):
        out[t][k]=data[t][k]
      else:
        out[t][k]=data[t][k].data

        # x and y may change with time, but should be the same for each time (gfs file)
        if not out[t].has_key('x'):
          out[t]['x']=data[t][k].x
          out[t]['y']=data[t][k].y

  return out, miss


def load_blkdata_wrf(wrfpath,date0=False,date1=False,quiet=True):
  from okean.datasets import wrf
  a=wrf.WRFData(wrfpath)
  data=a.data(date0,date1,quiet)

  out=cb.odict()
  time=data['time']
  for it in range(len(time)):
    out[time[it]]={}
    for k in data.keys():

      if k in ('time',): continue
      elif  k.startswith('INFO'):
        out[time[it]][k]=data[k]
      else:
        out[time[it]][k]=data[k].data[it,...]

        # x and y should be the same, so:
        if not out[time[it]].has_key('x'):
          out[time[it]]['x']=data[k].x
          out[time[it]]['y']=data[k].y


  return out


def load_blkdata_ncep(): pass


def load_blkdata_narr(date0,date1,quiet=True):
  from okean.datasets import narr
  a=narr.NARRData()
  data,miss=a.data(date0,date1,quiet=quiet)

  out=cb.odict()
  time=data.keys()
  for t in time:
    out[t]={}
    for k in data[t].keys():

      if  k.startswith('INFO'):
        out[t][k]=data[t][k]
      else:
        out[t][k]=data[t][k].data
        if not out[t].has_key('x'):
          out[t]['x']=data[t][k].x
          out[t]['y']=data[t][k].y

  return out, miss



def load_blkdata_interim(interimpath,date0=False,date1=False,quiet=True,past=False):
  '''
  If past, then module interim_past will be used to extract data from
  old interim server. Otherwise interim module will be used.
  '''
  if past:
    from okean.datasets import interim_past as interim
  else:
    from okean.datasets import interim

  a=interim.INTERIMData(interimpath)
  data=a.data(date0,date1,quiet)

  out=cb.odict()
  time=data['time']
  for it in range(len(time)):
    out[time[it]]={}
    for k in data.keys():

      if k in ('time',): continue
      elif  k.startswith('INFO'):
        out[time[it]][k]=data[k]
      else:
        out[time[it]][k]=data[k].data[it,...]

        # x and y should be the same, so:
        if not out[time[it]].has_key('x'):
          out[time[it]]['x']=data[k].x
          out[time[it]]['y']=data[k].y


  return out

def load_blkdata_cfsr(cfsrpath,date0=False,date1=False,quiet=True):
  from okean.datasets import cfsr
  a=cfsr.CFSRData(cfsrpath)
  data=a.data(date0,date1,quiet)

  out=cb.odict()
  time=data['time']
  for it in range(len(time)):
    out[time[it]]={}
    for k in data.keys():

      if k in ('time',): continue
      elif  k.startswith('INFO'):
        out[time[it]][k]=data[k]
      else:
        out[time[it]][k]=data[k].data[it,...]

        # x and y should be the same, so:
        if not out[time[it]].has_key('x'):
          out[time[it]]['x']=data[k].x
          out[time[it]]['y']=data[k].y


  return out




def conv_units(data,model,quiet=True):
  '''
  Some bulk variables have different units in ROMS Rutgers and ROMS-AGRIF
  Type can thus be 'roms' or 'roms-agrif'

  variable   units      roms     roms-agrif    description
  ----------------------------------------------------------------
    tair      C           C          C         surface temperature
    rhum     0--1       0--100      0--1       humidity
    pres      Pa         mb          Pa        surface pressure
    prate    cm/day   kg m-2 s-1    cm/d       total precipitation
    radsw    W m-2       W m-2      W m-2      sw radiation
    radlw    W m-2      (-)W m-2    W m-2      lw radiation
    dlwrf    W m-2      (-)W m-2    W m-2      down lw radiation
    uwnd     m s-1       m s-1      m s-1      u 10 m
    vwnd     m s-1       m s-1      m s-1      v 10 m
    ustress   Pa          Pa         Pa        u stress
    vstress   Pa          Pa         Pa        v stress
    cloud    0--1        0--1       ------     cloud cover
  ----------------------------------------------------------------
  '''

  if model=='roms':
    if not quiet: print '    rhum',
    data['rhum']  = data['rhum']*100                # 0--1  --> 0--100
    if not quiet: print ' pres',
    data['pres']  = data['pres']/100.               # Pa --> mb
    if not quiet: print ' prate',
    data['prate'] = data['prate']*1000/(86400*100)  # cm day-1 --> kg m-2 s-1
    if not quiet: print ' radlw'
    data['radlw'] = data['radlw']*-1 # positive down
    if not quiet: print ' dlwrf',
    data['dlwrf'] = data['dlwrf']*-1 # positive down


def data2romsblk(data,grd,**kargs):
  quiet        = True
  keepOriginal = 2 # if 2, will keep data only inside rt.grid_vicinity rectange
                   # if 1, all original data is kept
  Margin       = 4 # for griddata and for original data to keep

  for k in kargs.keys():
    if k=='quiet': quiet=kargs[k]
    elif k.lower().startswith('keepor'): keepOriginal=kargs[k]
    elif k=='margin': Margin=kargs[k]


  g=roms.Grid(grd)

  cond=False
  out={}

  for vname in data.keys():
    if vname.startswith('INFO') or vname in 'xy': continue
    # vnames starting with INFO are an info key, without data

    if not quiet: print '  interp %s' % vname

    if cond is False:
      cond,inds=rt.grid_vicinity(grd,data['x'],data['y'],
                                 margin=Margin,rect=True,retinds=True)
      i1,i2,j1,j2=inds

    out[vname]=calc.griddata(data['x'][cond],data['y'][cond],
                           data[vname][cond],g.lon,g.lat,extrap=True)

    if keepOriginal:
      # keep original data:
      if keepOriginal==1:
        out[vname+'_original']=data[vname] # keep all original data
      elif keepOriginal==2:
        out[vname+'_original']=data[vname][j1:j2,i1:i2] # keep data only inside vicinity rectangle

      if not out.has_key('x_original'):
        if keepOriginal==1:
          out['x_original']=data['x']
          out['y_original']=data['y']
        elif keepOriginal==2:
          out['x_original']=data['x'][j1:j2,i1:i2]
          out['y_original']=data['y'][j1:j2,i1:i2]

  # about wind:
  if not quiet: print ' --> rot U,V wind and U,V wind stress'
  angle=g.use('angle')
  out['uwnd'],out['vwnd']=calc.rot2d(out['uwnd'],out['vwnd'],angle)
  out['sustr'],out['svstr']=calc.rot2d(out['sustr'],out['svstr'],angle)

  out['sustr']=rt.rho2uvp(out['sustr'],'u')
  out['svstr']=rt.rho2uvp(out['svstr'],'v')


  return out


def make_blk_interim(interimpath,grd,bulk,date0=False,date1=False,**kargs):
  '''
  kargs:
    tunits, ex: 'seconds since yyyy-mm-dd'
    quiet
    create
    model, 'roms' or 'roms-agrif'
    title
    keeporiginal, save original data in nc file
    margin, for original data to be saved)
    attr, additional global arguments, ex: attr={'some_info': 'abc','what': 123}
    past,  use True for old data server data, default False
    ...
  '''

  quiet  = kargs.get('quiet',0)
  create = kargs.get('create',1)
  model  = kargs.get('model','roms') # or roms-agrif
  past   = kargs.get('past',False) # use True for old data server data

  data=load_blkdata_interim(interimpath,date0,date1,quiet,past)

  # about original data, run data2romsblk once to test for x_original:
  tmp=data2romsblk(data[data.keys()[0]],grd,**kargs)
  if 'x_original' in tmp.keys(): original=tmp['x_original'].shape
  else: original=False

  q=gennc.GenBlk(bulk,grd,**kargs)
  if create: q.create(model,original)

  for d in data.keys():
    if model=='roms':
       if not quiet: print '  converting units:',
       conv_units(data[d],model,quiet)

    D=data2romsblk(data[d],grd,**kargs)
    D['date']=d

    if not quiet: print '  =>filling date=%s' % d.isoformat(' ')
    q.fill(D,quiet=quiet)


def make_blk_gfs(gfspath,grd,bulk,date0,date1=False,nforec=0,**kargs):
  '''
  see make_blk_interim
  '''

  quiet  = kargs.get('quiet',0)
  create = kargs.get('create',1)
  model  = kargs.get('model','roms') # or roms-agrif

  data,miss=load_blkdata_gfs(gfspath,date0,date1,nforec=nforec,quiet=quiet)

  # about original data, run data2romsblk once to test for x_original:
  tmp=data2romsblk(data[data.keys()[0]],grd,**kargs)
  if 'x_original' in tmp.keys(): original=tmp['x_original'].shape
  else: original=False

  # add gfs INFO to file global attributes:
  if not 'attr' in kargs.keys(): kargs['attr']={}
  s=''
  for k in data.keys(): s+=' # '+k.isoformat(' ')+' '+data[k]['INFO_file']+' isbest '+str(data[k]['INFO_isbest'])
  kargs['attr']['sources']=s

  # common to interim----------------------
  q=gennc.GenBlk(bulk,grd,**kargs)
  if create: q.create(model,original)

  for d in data.keys():
    if model=='roms':
       if not quiet: print '  converting units:',
       conv_units(data[d],model,quiet)

    D=data2romsblk(data[d],grd,**kargs)
    D['date']=d

    if not quiet: print '  =>filling date=%s' % d.isoformat(' ')
    q.fill(D,quiet=quiet)


def make_blk_narr(grd,bulk,date0,date1=False,**kargs):
  '''
  see make_blk_interim

  ps: to not include original data in file, use kargg keepor=False
  '''

  quiet  = kargs.get('quiet',0)
  create = kargs.get('create',1)
  model  = kargs.get('model','roms') # or roms-agrif

  data,miss=load_blkdata_narr(date0,date1,quiet=quiet)

  # about original data, run data2romsblk once to test for x_original:
  tmp=data2romsblk(data[data.keys()[0]],grd,**kargs)
  if 'x_original' in tmp.keys(): original=tmp['x_original'].shape
  else: original=False

  # add narr INFO to file global attributes:
  if not 'attr' in kargs.keys(): kargs['attr']={}
  try:
    s=''
    for k in data.keys(): s+=' # '+k.isoformat(' ')+' '+data[k]['INFO_file']+' isbest '+str(data[k]['INFO_isbest'])
    kargs['attr']['sources']=s
  except: pass

  # common to interim----------------------
  q=gennc.GenBlk(bulk,grd,**kargs)
  if create: q.create(model,original)

  for d in data.keys():
    if model=='roms':
       if not quiet: print '  converting units:',
       conv_units(data[d],model,quiet)

    D=data2romsblk(data[d],grd,**kargs)
    D['date']=d

    if not quiet: print '  =>filling date=%s' % d.isoformat(' ')
    q.fill(D,quiet=quiet)


def make_blk_cfsr(cfsrpath,grd,bulk,date0=False,date1=False,**kargs):
  '''
  see make_blk_interim

  ps: to not include original data in file, use kargg keepor=False
  '''

  quiet  = kargs.get('quiet',0)
  create = kargs.get('create',1)
  model  = kargs.get('model','roms') # or roms-agrif

  data=load_blkdata_cfsr(cfsrpath,date0,date1,quiet) # unique diff from make_blk_interim !!

  # about original data, run data2romsblk once to test for x_original:
  tmp=data2romsblk(data[data.keys()[0]],grd,**kargs)
  if 'x_original' in tmp.keys(): original=tmp['x_original'].shape
  else: original=False

  q=gennc.GenBlk(bulk,grd,**kargs)
  if create: q.create(model,original)

  for d in data.keys():
    if model=='roms':
       if not quiet: print '  converting units:',
       conv_units(data[d],model,quiet)

    D=data2romsblk(data[d],grd,**kargs)
    D['date']=d

    if not quiet: print '  =>filling date=%s' % d.isoformat(' ')
    q.fill(D,quiet=quiet)

def make_blk_wrf(wrfpath,grd,bulk,date0=False,date1=False,**kargs):
  '''
  see make_blk_interim

  '''

  quiet  = kargs.get('quiet',0)
  create = kargs.get('create',1)
  model  = kargs.get('model','roms') # or roms-agrif

  data=load_blkdata_wrf(wrfpath,date0,date1,quiet) # unique diff from make_blk_interim !!

  # about original data, run data2romsblk once to test for x_original:
  tmp=data2romsblk(data[data.keys()[0]],grd,**kargs)
  if 'x_original' in tmp.keys(): original=tmp['x_original'].shape
  else: original=False

  q=gennc.GenBlk(bulk,grd,**kargs)
  if create: q.create(model,original)

  for d in data.keys():
    if model=='roms':
       if not quiet: print '  converting units:',
       conv_units(data[d],model,quiet)

    D=data2romsblk(data[d],grd,**kargs)
    D['date']=d

    if not quiet: print '  =>filling date=%s' % d.isoformat(' ')
    q.fill(D,quiet=quiet)


# --------------------------- OLD stuff from here... to be removed soon --------
def __update_wind(fname,datapath,source,**kargs):
  if source=='quikscat':
    new_wind_info='wind from quikscat'
    from okean.datasets import quikscat
    a=quikscat.WINDData(datapath)
  elif source=='blended':
    new_wind_info='wind from myocean blended'
    from okean.datasets import blended_wind
    a=blended_wind.WINDData(datapath)

  time=netcdf.nctime(fname,'time')
  date0=dts.next_date(time[0],-1)
  date1=dts.next_date(time[-1],+2)
  data=a.data(date0,date1)

  update_wind(fname,data,new_wind_info,**kargs)

def update_wind_quikscat(fname,datapath,**kargs):
  __update_wind(fname,datapath,'quikscat',**kargs)

def update_wind_blended(fname,datapath,**kargs):
  __update_wind(fname,datapath,'blended',**kargs)

def update_wind_blended2(fname,datapaths,**kargs):
  '''
  In days without blended data will try to use quikscat data
  '''
  from okean.datasets import quikscat
  from okean.datasets import blended_wind
  a=blended_wind.WINDData(datapaths[0])
  b=quikscat.WINDData(datapaths[1])

  time=netcdf.nctime(fname,'time')
  date0=dts.next_date(time[0],-1)
  date1=dts.next_date(time[-1],+2)

  data=a.data(date0,date1)

  # limit are... otherwise, quikscat interp will be very slow!
  grd=netcdf.fatt(fname,'grd_file')
  import os
  if not os.path.isfile(grd): grd=kargs['grd']
  cond,inds=rt.grid_vicinity(grd,data['x'],data['y'],margin=5,rect=True,retinds=True)
  i1,i2,j1,j2=inds
  for d in data.keys():
    if   d == 'x': data[d]=data[d][i1:i2]
    elif d == 'y': data[d]=data[d][j1:j2]
    else: data[d]=data[d][j1:j2,i1:i2]


  # check for missing days:
  time0=data.keys()
  x0=data['x']
  y0=data['y']
  x0,y0=np.meshgrid(x0,y0)
  time0.remove('x')
  time0.remove('y')

  out=cb.odict()
  out['x']=x0
  out['y']=y0
  info=''
  qs_ij_limits_done=False
  for d in dts.drange(date0,date1):
    found=0
    for t in time0:
      if (t.year,t.month,t.day)==(d.year,d.month,d.day):
        print '==> blended : ',t
        out[t]=data[t]
        found=1

    if not found: # use quikscat:
      print '==> quikscat : ',d.strftime('%Y-%m-%d')
      tmp= b.data(d,dts.next_date(d))
      if not tmp.has_key('x'): continue
      x,y=tmp['x'],tmp['y']
      x,y=np.meshgrid(x,y)

      # reduce qs data:
      if not qs_ij_limits_done:
        i1,i2,j1,j2=calc.ij_limits(x,y,[x0.min(),x0.max()],[y0.min(),y0.max()])
        qs_ij_limits_done=True

      x=x[j1:j2,i1:i2]
      y=y[j1:j2,i1:i2]
      tmp[tmp.keys()[0]]=tmp[tmp.keys()[0]][j1:j2,i1:i2]


      print '  griddata u'
      u=calc.griddata(x,y,tmp[tmp.keys()[0]].real,x0,y0)
      print '  griddata v'
      v=calc.griddata(x,y,tmp[tmp.keys()[0]].imag,x0,y0)
      out[tmp.keys()[0]]=u+1.j*v
      info+='#'+d.strftime('%Y%m%d')


  new_wind_info='blended+quikscat at days: '+info
  update_wind(fname,out,new_wind_info,**kargs)


def update_wind(fname,data,new_wind_info,**kargs):
  '''
  new_wind_info will be added to fname as global attribute
  (ex: 'new wind from database xxx')

  '''
  quiet  = False
  Margin = 4 # for griddata and for original data to keep
  grd    = False
  keepOriginal = 2

  for k in kargs.keys():
    if   k=='quiet':  quiet=kargs[k]
    elif k=='margin': Margin=kargs[k]
    elif k=='grid':   grd=kargs[k]
    elif k.lower().startswith('keepor'): keepOriginal=kargs[k]

  if not grd: grd=netcdf.fatt(fname,'grd_file')
  g=roms.Grid(grd)

  x0=data['x']
  y0=data['y']
  if x0.ndim==1: x0,y0=np.meshgrid(x0,y0)
  dates=data.keys()[:]
  dates.remove('x')
  dates.remove('y')
  dates=np.array(dates)


  # interp in time:
  time=netcdf.nctime(fname,'time')
  cond=False
  tind=-1
  fob=gennc.GenBlk(fname,grd)
  for t in time:
    newWind={}
    tind+=1

    I,=np.where(dates==t)
    if I.size:
      I=I[0]
      uv=data[dates[I]]
      if not quiet: print t,dates[I]
    else:
      i1,=np.where(dates>t)
      i0,=np.where(dates<t)
      if i0.size and i1.size:
        i1=i1[0]
        i0=i0[-1]
        d0=t-dates[i0]
        d1=dates[i1]-t
        d0=d0.days+d0.seconds/86400.
        d1=d1.days+d1.seconds/86400.
        uv=(data[dates[i0]]*d1+data[dates[i1]]*d0)/(d0+d1)
        if not quiet: print t,dates[i0],dates[i1], d0,d1
      elif not i1.size:
        uv=data[dates[-1]]
        if not quiet: print t,dates[-1]
      elif not i0.size:
        uv=data[dates[0]]
        if not quiet: print t,dates[0]

    # interp to grid:
    if cond is False:
      cond,inds=rt.grid_vicinity(grd,x0,y0,margin=Margin,rect=True,retinds=True)
      i1,i2,j1,j2=inds

    if not quiet: print ' --> inter uv %s' % t.isoformat(' ')
    u=calc.griddata(x0[cond],y0[cond],uv.real[cond],g.lon,g.lat,extrap=True)
    v=calc.griddata(x0[cond],y0[cond],uv.imag[cond],g.lon,g.lat,extrap=True)


    # rotate wind, calc stress:
    if not quiet: print ' --> rot U,V wind and U,V wind stress'
    wspd=np.sqrt(u**2+v**2)
    sustr,svstr=air_sea.wind_stress(u,v)
    angle=g.use('angle')
    uwnd,vwnd=calc.rot2d(u,v,angle)
    sustr,svstr=calc.rot2d(sustr,svstr,angle)
    sustr=rt.rho2uvp(sustr,'u')
    svstr=rt.rho2uvp(svstr,'v')

    # update wind data:
    newWind['date']  = t
    newWind['uwnd']  = uwnd
    newWind['vwnd']  = vwnd
    newWind['sustr'] = sustr
    newWind['svstr'] = svstr
    newWind['wspd']  = wspd
    # original xy:
    if tind==0:
      newWind['attr']={'new_wind_info':new_wind_info}
      if not quiet: print ' --> original xy'
      if keepOriginal==1:
        newWind['x_wind']=x0
        newWind['y_wind']=y0
      elif keepOriginal==2:
        newWind['x_wind']=x0[j1:j2,i1:i2]
        newWind['y_wind']=y0[j1:j2,i1:i2]


    # add to file:
    if not quiet: print ' --> adding to file'
    fob.update_wind(newWind,quiet=quiet)
    if not quiet: print ''
