from os.path import isfile
import numpy as np
import datetime

from okean import calc, cookbook as cb, netcdf
from okean.roms import roms_tools as rt, roms
import gennc

class DataAccess:
  # dims:
  xdim = 'lon'
  ydim = 'lat'
  zdim = 'depth'
  tdim = 'time'

  # vars xyz:
  x_name = 'lon'
  y_name = 'lat'
  z_name = 'depth'

  # TSUV
  temp_name = 'temp'
  salt_name = 'salt'
  u_name    = 'u'
  v_name    = 'v'
  ssh_name  = 'ssh'

  # time:
  time_name  = 'time'

  def __init__(self,**kargs):
    for k in kargs.keys():
      if not hasattr(self,k): print('new key %s' % k )
      setattr(self,k,kargs[k])

def load_data(f,quiet=0,**kargs):
  '''
  Loads prognostic variables (temp,salt,u,v,ubar,vbar,zeta) from
  netcdf file or opendap server. Also loads lon,lat, depth, and time.

  If f is a file, it must include the 1d variables lon,lat and depth;
  the 2d variable ssh (zeta) and the 3d variables temp, salt, u and v;
  ie, each file must contain data for a simgle time. The file must also
  contain the variable time.

  If f is a opendap address, it must contain also all these variables
  or the ones defined in the input karg settings (DataAccess object)

  To deal with the case of variables in different files/opendap addresses,
  f can also be a dictionary with keys the variables and values the files
  or opendap addresses. In this case, the keys must be:
    - temp
    - salt
    - u
    - v
    - ssh
    - misc, for lon, lat, depth, time and dimensions
      or xy for lon,lat and x,ydim; z for depth and zdim, time for time

  The output data (dict) is suitable to be used by data2roms, which
  interpolates the data to ROMS 3d grid.
  Also outputs an error/status string.

  kargs:
    inds, dict with dimension names/values (where time dim can be integer
          or datetime)
    settings, DataAccess object
    extra, extra misc vars to load [(outKey0,fileVar0),...]
    t_units, units of variable time, by default the att  units is used
  '''

  sett=DataAccess()
  inds={}
  extra=[]
  t_units=[]
  if 'settings' in kargs.keys(): sett    = kargs['settings']
  if 'inds'     in kargs.keys(): inds    = kargs['inds']
  if 'extra'    in kargs.keys(): extra   = kargs['extra']
  if 't_units'  in kargs.keys(): t_units = kargs['t_units']

  res={}
  msg=''

  if not isinstance(f,dict) and not f.startswith('http') and not isfile(f):
    msg='file not found %s' % f
    if not quiet: print(msg)
    return res, msg

  # load nc files:
  if not isinstance(f,dict):
    f={'temp':f,'salt':f,'u':f,'v':f,'ssh':f,'misc':f}

  if not f.has_key('xy'):   f['xy']   = f['misc']
  if not f.has_key('z'):    f['z']    = f['misc']
  if not f.has_key('time'): f['time'] = f['misc']

  filesUsed=[]
  ncUsed=[]
  for i in f.keys():
    if not quiet: print('(%s) loading from %s' % (i.ljust(5),f[i]))

    if i=='temp':
      if f[i] in filesUsed: ncTemp=ncUsed[filesUsed.index(f[i])]
      else:
        ncTemp=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncTemp]

    elif i=='salt':
      if f[i] in filesUsed: ncSalt=ncUsed[filesUsed.index(f[i])]
      else:
        ncSalt=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncSalt]

    elif i=='u':
      if f[i] in filesUsed: ncU=ncUsed[filesUsed.index(f[i])]
      else:
        ncU=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncU]

    elif i=='v':
      if f[i] in filesUsed: ncV=ncUsed[filesUsed.index(f[i])]
      else:
        ncV=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncV]

    elif i=='ssh':
      if f[i] in filesUsed: ncSsh=ncUsed[filesUsed.index(f[i])]
      else:
        ncSsh=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncSsh]

    elif i=='xy':
      if f[i] in filesUsed: ncXy=ncUsed[filesUsed.index(f[i])]
      else:
        ncXy=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncXy]

    elif i=='z':
      if f[i] in filesUsed: ncZ=ncUsed[filesUsed.index(f[i])]
      else:
        ncZ=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncZ]

    elif i=='time':
      if f[i] in filesUsed: ncTime=ncUsed[filesUsed.index(f[i])]
      else:
        ncTime=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncTime]

    elif i=='misc':
      if f[i] in filesUsed: ncMisc=ncUsed[filesUsed.index(f[i])]
      else:
        ncMisc=netcdf.ncopen(f[i])
        filesUsed+=[f[i]]
        ncUsed+=[ncMisc]


  # load dims:
  if not quiet: print('  loading dims...')
  dimsXy=netcdf.fdim(ncXy)
  dimsZ =netcdf.fdim(ncZ)

  res['NX']=dimsXy[sett.xdim]
  res['NY']=dimsXy[sett.ydim]
  ###if sett.z_name:
  if sett.zdim:
    res['NZ']=dimsZ[sett.zdim]
  else:
    res['NZ']=1

  # about horizontal inds:
  if inds.has_key(sett.xdim) and len(inds[sett.xdim])==2 and not isinstance(inds[sett.xdim],basestring):
    if not quiet: print('  calc horizontal inds...')
    xlim=inds[sett.xdim]
    ylim=inds[sett.ydim]

    inds.pop(sett.xdim)
    inds.pop(sett.ydim)

    lon=netcdf.use(ncXy,sett.x_name,**inds)
    if np.any(lon>360): lon=np.mod(lon,360.)
    lat=netcdf.use(ncXy,sett.y_name,**inds)
    i0,i1,j0,j1=calc.ij_limits(lon,lat,xlim,ylim,margin=3)
    inds[sett.xdim]='%d:%d' % (i0,i1)
    inds[sett.ydim]='%d:%d' % (j0,j1)


  if not quiet: print('  loading lon, lat, depth...')
  res['lon']  = netcdf.use(ncXy,sett.x_name,**inds)
  if np.any(res['lon']>360): res['lon']=np.mod(res['lon'],360.)
  res['lat']  = netcdf.use(ncXy,sett.y_name,**inds)
  if sett.z_name:
    res['depth'] = -netcdf.use(ncZ,sett.z_name,**inds)
  else: res['depth']=False

  if res['lon'].size!=res['lat'].size:
    res['lon'],res['lat']=np.meshgrid(res['lon'],res['lat'])
    # needed for griddata, later

  # update nx,ny:
  if inds.has_key(sett.xdim):
    res['NY'],res['NX']=res['lon'].shape

  # extra misc vars:
  if len(extra):
    for outKey,fileVar in extra:
      if not quiet: print('  loading extra misc... %s %s' % (outKey,fileVar))
      res[outKey]=netcdf.use(ncMisc,fileVar,**inds)


  # time:
  # file may have one or several times. If several, time dim must be given
  # with kargs inds!
  # but file may also have no time dim or time name !
  if sett.time_name:
    if not quiet: print('  loading time...')
    if t_units:
      times=netcdf.use(ncTime,sett.time_name)
      times=netcdf.num2date(times,t_units)
    else:
      times=netcdf.nctime(ncTime,sett.time_name)

    if inds.has_key(sett.tdim):
      try: tind=dts.parse_date(inds[sett.tdim])
      except: tind=inds[sett.tdim] # is an integer, for instance

      if isinstance(tind,datetime.datetime):
        tind,=np.where(times==tind)
        if tind.size:
          tind=tind[0]
          inds[sett.tdim]=tind # update inds to extract other variables
        else:
          Msg='date not found'
          msg+='\n'+Msg
          return res,msg+' ERROR'

      date=times[tind]
      try:
        len(date)
        ndates=True
      except: ndates=False

      if ndates:
        if not quiet: print('    tind, date= len=%d: %d to %d, %s to %s' % (len(date),tind[0],tind[-1],date[0].isoformat(' '),date[-1].isoformat(' ')))
      else:
        if not quiet: print('    tind, date= %d %s' % (tind,date.isoformat(' ')))

    elif times.size==1:
      date=times[0]
      if not quiet: print('    date= %s' % date.isoformat(' '))
    else: # must provide tind as input!!
      Msg='several dates in file... provice tind!'
      msg+='\n'+Msg
      return res,msg+' ERROR'

    res['date'] = date
  else:
    if not quiet: print('    warning: not using time !!')
    res['date']=0

  empty3d=np.zeros([res['NZ'],res['NY'],res['NX']])
  empty2d=np.zeros([res['NY'],res['NX']])

  if 'temp' in f.keys():
    if not quiet: print('  loading temp...')
    if sett.temp_name in ncTemp.varnames: res['temp'] = netcdf.use(ncTemp,sett.temp_name,**inds)
    else:
      Msg='var %s not found' % 'temp'
      msg+='\n'+Msg
      if not quiet: print(Msg)
      res['temp']=empty3d

  if 'salt' in f.keys():
    if not quiet: print('  loading salt...')
    if sett.salt_name in ncSalt.varnames: res['salt'] = netcdf.use(ncSalt,sett.salt_name,**inds)
    else:
      Msg='var %s not found' % 'salt'
      msg+='\n'+Msg
      if not quiet: print(Msg)
      res['salt']=empty3d

  if 'u' in f.keys():
    if not quiet: print('  loading u...')
    if sett.u_name in ncU.varnames: res['u']    = netcdf.use(ncU,sett.u_name,**inds)
    else:
      Msg='var %s not found' % 'u'
      msg+='\n'+Msg
      if not quiet: print(Msg)
      res['u']=empty3d

  if 'v' in f.keys():
    if not quiet: print('  loading v...')
    if sett.v_name in ncV.varnames: res['v']    = netcdf.use(ncV,sett.v_name,**inds)
    else:
      Msg='var %s not found' % 'v'
      msg+='\n'+Msg
      if not quiet: print(Msg)
      res['v']=empty3d

  if 'ssh' in f.keys():
    if not quiet: print('  loading ssh...')
    if sett.ssh_name in ncSsh.varnames: res['ssh']  = netcdf.use(ncSsh,sett.ssh_name,**inds)
    else:
      Msg='var %s not found' % 'ssh'
      msg+='\n'+Msg
      if not quiet: print(Msg)
      res['ssh']=empty2d

  for nc in ncUsed:
    try:  nc.close()
    except: pass

  return res, msg


def avg_dates(Date,data,data1):
  date  = data['date']
  date1 = data1['date']
  print('averaging to %s (%s, %s)' % (Date.isoformat(' '),date.isoformat(' '),date1.isoformat(' ')))
  a=Date-date
  b=date1-Date

  a=a.days+a.seconds/86400.
  b=b.days+b.seconds/86400.

  print('    ',a,b)

  res=data.copy()
  res['date']=Date
  for i in ('temp','salt','u','v','zeta','ubar','vbar'):
    print('    avg %s' % i)
    res[i] = (data[i]*b +data1[i]*a)/(a+b)

  return res


def data2z(data,**kargs):
  '''
  Interpolates data dict defined at 3D vertical coordinate (like sigma
  levels or s-levels, for instance) to 1D depths (layers of z constant)

  data must contain all the fileds from load_data, z3d and the new
  depth, ie, data['depth'].ndim==1

  **kargs:
  ij : axis for vertical interpolations (*i,j)
  quiet : output messages flag (True by default)
  rep_surf: repeat surface level, ie, add a new surface level with same
    data as the old surface. This will ensure extrapolations in vertical
    slices do not occur horizontally, ie, misture will occur in the vertical
    not along the slice (True by default)
  interp_opts: options for griddata
  '''

  ij=kargs.get('ij','j')
  quiet=kargs.get('quiet',True)
  rep_surf=kargs.get('rep_surf',True)
  interp_opts=kargs.get('interp_opts',{})

  # repeat surface:
  if rep_surf:
    for vname in ['z3d','temp','salt','u','v']:
      if np.ma.isMA(data[vname]):
        data[vname]=np.ma.vstack((data[vname],data[vname][-1][np.newaxis,...]))
      else:
        data[vname]=np.vstack((data[vname],data[vname][-1][np.newaxis,...]))

    data['z3d'][-1]=data['z3d'].max()+1
    data['NZ']=data['NZ']+1


  Z = data['depth'] # new depths 1D
  z = data['z3d'] # 3D depths
  nZ=len(Z)
  nz=data['NZ']
  nx=data['NX']
  ny=data['NY']

  # lon, lat 2d:
  x0=data['lon']
  y0=data['lat']
  if x0.ndim==1:x0,y0 = np.meshgrid(x0,y0)

  # lon,lat 3d:
  x  = np.tile(x0,(nz,1,1))
  y  = np.tile(y0,(nz,1,1))
  xx = np.tile(x0,(nZ,1,1))
  yy = np.tile(y0,(nZ,1,1))

  # new z 3d:
  zz=np.tile(Z[:,np.newaxis,np.newaxis],(ny,nx))

  # new data vars:
  Temp = np.ma.masked_all((nZ,ny,nx),data['temp'].dtype)
  Salt = np.ma.masked_all((nZ,ny,nx),data['salt'].dtype)
  U    = np.ma.masked_all((nZ,ny,nx),data['u'].dtype)
  V    = np.ma.masked_all((nZ,ny,nx),data['v'].dtype)

  if ij=='j':
    for j in range(ny):
      if not quiet and j%10==0: print('z interp (%d z levels) j=%d of %d' % (nZ,j,ny))

      v=data['temp'][:,j,:]
      v=np.ma.masked_where(v==0,v)
      Temp[:,j,:]=calc.griddata(x[:,j,:],z[:,j,:],v,xx[:,j,:],zz[:,j,:],extrap=True,**interp_opts)

      v=data['salt'][:,j,:]
      v=np.ma.masked_where(v==0,v)
      Salt[:,j,:]=calc.griddata(x[:,j,:],z[:,j,:],v,xx[:,j,:],zz[:,j,:],extrap=True,**interp_opts)

      v=data['u'][:,j,:]
      v=np.ma.masked_where(v==0,v)
      U[:,j,:]=calc.griddata(x[:,j,:],z[:,j,:],v,xx[:,j,:],zz[:,j,:],extrap=True,**interp_opts)

      v=data['v'][:,j,:]
      v=np.ma.masked_where(v==0,v)
      V[:,j,:]=calc.griddata(x[:,j,:],z[:,j,:],v,xx[:,j,:],zz[:,j,:],extrap=True,**interp_opts)
  elif ij=='i':
    for i in range(nx):
      if not quiet and i%10==0: print('z interp (%d z levels) i=%d of %d' % (nZ,i,nx))

      v=data['temp'][:,:,i]
      v=np.ma.masked_where(v==0,v)
      Temp[:,:,i]=calc.griddata(y[:,:,i],z[:,:,i],v,yy[:,:,i],zz[:,:,i],extrap=True,**interp_opts)

      v=data['salt'][:,:,i]
      v=np.ma.masked_where(v==0,v)
      Salt[:,:,i]=calc.griddata(y[:,:,i],z[:,:,i],v,yy[:,:,i],zz[:,:,i],extrap=True,**interp_opts)

      v=data['u'][:,:,i]
      v=np.ma.masked_where(v==0,v)
      U[:,:,i]=calc.griddata(y[:,:,i],z[:,:,i],v,yy[:,:,i],zz[:,:,i],extrap=True,**interp_opts)

      v=data['v'][:,:,i]
      v=np.ma.masked_where(v==0,v)
      V[:,:,i]=calc.griddata(y[:,:,i],z[:,:,i],v,yy[:,:,i],zz[:,:,i],extrap=True,**interp_opts)


  res = data.copy()
  res['temp'] = Temp
  res['salt'] = Salt
  res['u']    = U
  res['v']    = V
  res['NZ']   = nZ
  res.pop('z3d')

  return res


def data2roms(data,grd,sparams,**kargs):
  '''
  Interpolates data to roms 3D grid.

  The dict data must contain the prognostic variables temp, salt, u,
  v (3d) and ssh (zeta, 2d), as well as lon, lat (2d), depth (1d) and
  time/date info: date (data date), date0 (reference date) and time
  (difference between date and date0). The input data can be provided
  by load_data.

  Parameters
  ----------
  data : dict with prognostic variables
  grd : ROMS netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels

  **kargs:
  ij : axis for vertical interpolations (*i,j)
  ij_ind : list of i or j index for vertical interpolation, all by
           default (ij_ind=False)
  horizAux : if True, the original data horizontally interpolated is
             returned and can be used for next data2roms call with
             this same karg
  quiet : output messages flag (false by default)
  proj : projection - False, name or basemap proj - lcc by default
         if False, horizontal interpolations will use lonxlat instead of distances
  interp_opts: options for griddata
  rep_surf: repeat surface level (new upper level)
  '''

  ij          = kargs.get('ij','j')
  ij_ind      = kargs.get('ij_ind',False)
  horizAux    = kargs.get('horizAux',False)
  quiet       = kargs.get('quiet',False)
  proj        = kargs.get('proj','lcc') # lonxlat to distance before
                                        # horizontal interpolation
  interp_opts = kargs.get('interp_opts',{})
  rep_surf    = kargs.get('rep_surf',True) # create a surface upper level
                                           # before interpolation

  if not quiet: print('using grid %s' % grd)
  g=roms.Grid(grd)
  xr,yr,hr,mr=g.vars('r')
  xu,yu,hu,mu=g.vars('u')
  xv,yv,hv,mv=g.vars('v')
  ny,nx=hr.shape
  nz=sparams[3]

  if proj:
    print('projecting coordinates...')
    if isinstance(proj,basestring):
       lonc=(xr.max()+xr.min())/2.
       latc=(yr.max()+yr.min())/2.
       from mpl_toolkits.basemap import Basemap
       proj=Basemap(projection=proj,width=1,height=1,resolution=None,
                    lon_0=lonc,lat_0=latc, lat_1=latc)

    xr,yr=proj(xr,yr)
    xu,yu=proj(xu,yu)
    xv,yv=proj(xv,yv)
    dlon,dlat=proj(data['lon'],data['lat'])
    Rdz=1/100. # distance to depth ratio (300km vs 3000m)
    distance=lambda x,y: np.append(0.,np.sqrt(np.diff(x)**2+np.diff(y)**2).cumsum())
  else:
    dlon,dlat=data['lon'],data['lat']
    distance=calc.distance

  # needed for s_levels and for rep_surf!
  sshr=calc.griddata(dlon,dlat,data['ssh'],xr,yr,extrap=True,**interp_opts)

  # repeat surface:
  if rep_surf:
    # copy data cos dont want to change the original dataset:
    import copy
    data=copy.deepcopy(data)
    for vname in ['temp','salt','u','v','depth']:
      if data[vname].ndim==1: # depth !
        if np.ma.isMA(data[vname]): vstack=np.ma.hstack
        else: vstack=np.hstack
      else:
        if np.ma.isMA(data[vname]): vstack=np.ma.vstack
        else: vstack=np.vstack

      if data['depth'][0]>data['depth'][1]: # surf at ind 0
        data[vname]=vstack((data[vname][0][np.newaxis],data[vname]))
        if vname=='depth': data[vname][0]=sshr.max()
      else:
        data[vname]=vstack((data[vname],data[vname][-1][np.newaxis]))
        if vname=='depth': data[vname][-1]=sshr.max()

    data['NZ']=data['NZ']+1

  NX=data['NX']
  NY=data['NY']
  NZ=data['NZ']

  if not quiet: print('calc s levels...')
  Zr = g.s_levels(sparams,sshr,hr,'rr')
  Zu = g.s_levels(sparams,sshr,hr,'ur')
  Zv = g.s_levels(sparams,sshr,hr,'vr')

  # interp horiz:
  retHorizAux=horizAux is True
  if horizAux in (True,False):
    TEMP = np.ma.masked_all((NZ,ny,nx),data['temp'].dtype)
    SALT = np.ma.masked_all((NZ,ny,nx),data['salt'].dtype)
    U    = np.ma.masked_all((NZ,ny,nx),data['u'].dtype)
    V    = np.ma.masked_all((NZ,ny,nx),data['v'].dtype)

    if not quiet: print('horizontal interpolation:')
    for i in range(NZ):
      if not quiet and i%10==0: print('   lev %d of %d' % (i,NZ))
      #import pylab
      #pylab.figure()
      #pylab.pcolormesh(data['lon'],data['lat'],data['temp'][i,...])

      try: TEMP[i,...] = calc.griddata(dlon,dlat,data['temp'][i,...],xr,yr,extrap=True,**interp_opts)
      except: pass

      try: SALT[i,...] = calc.griddata(dlon,dlat,data['salt'][i,...],xr,yr,extrap=True,**interp_opts)
      except: pass

      try: U[i,...] = calc.griddata(dlon,dlat,data['u'][i,...],xr,yr,extrap=True,**interp_opts)
      except: pass

      try: V[i,...] = calc.griddata(dlon,dlat,data['v'][i,...],xr,yr,extrap=True,**interp_opts)
      except: pass

    # rotate U,V:
    if not quiet: print('rotating U,V to grid angle')
    angle=g.use('angle')  # rad
    U,V=calc.rot2d(U,V,angle)
    U=rt.rho2uvp3d(U,'u')
    V=rt.rho2uvp3d(V,'v')

    horizAux={}
    horizAux['TEMP'] = TEMP
    horizAux['SALT'] = SALT
    horizAux['U']    = U
    horizAux['V']    = V

  else:
    TEMP = horizAux['TEMP']
    SALT = horizAux['SALT']
    U    = horizAux['U']
    V    = horizAux['V']

  # interp vert:
  nxu=nx-1
  nyv=ny-1
  #> -----------------------------------------------------------------
  useInd=not ij_ind is False
  if ij_ind is False:
    if   ij=='j': ij_ind=range(ny)
    elif ij=='i': ij_ind=range(nx)
  else:
    try: iter(ij_ind)
    except: ij_ind=[ij_ind]

    if   ij=='j': ny=nyv=len(ij_ind)
    elif ij=='i': nx=nxu=len(ij_ind)
  # -----------------------------------------------------------------<

  Temp = np.zeros((nz,ny ,nx ),data['temp'].dtype)
  Salt = np.zeros((nz,ny ,nx ),data['salt'].dtype)
  Uvel = np.zeros((nz,ny ,nxu),data['u'].dtype)
  Vvel = np.zeros((nz,nyv,nx ),data['v'].dtype)


  jslice=lambda x,ind: x[:,ind,:]
  islice=lambda x,ind: x[:,:,ind]

  ZZr = np.tile(data['depth'],(nx,ny,1)).T
  ZZu = np.tile(data['depth'],(nxu,ny,1)).T
  ZZv = np.tile(data['depth'],(nx,nyv,1)).T

  if not useInd is False: #>------------------------------------------
    if   ij=='j':
      slice=jslice
      sshr=sshr[ij_ind,:]
      hr  =hr[ij_ind,:]
    elif ij=='i':
      slice=islice
      sshr=sshr[:,ij_ind]
      hr  =hr[:,ij_ind]

    Zr,Zu,Zv,TEMP,SALT,U,V=[slice(k,ij_ind) for k in [Zr,Zu,Zv,TEMP,SALT,U,V]]
  # -----------------------------------------------------------------<

  if useInd: # then store distances for a possible bry file
    dtype=Temp.dtype
    distr=np.zeros((nz,ny, nx ),dtype)
    distu=np.zeros((nz,ny, nxu),dtype)
    distv=np.zeros((nz,nyv,nx ),dtype)

  if not quiet: print('vertical interpolation:')
  if ij=='j':
    for j in range(ny):
      if not quiet and (ny<10 or (ny>=10 and j%10==0)): print('  j=%3d of %3d' % (j,ny))
      ind=ij_ind[j]
      dr=np.tile(distance(xr[ind,:],yr[ind,:]),(nz,1))
      du=np.tile(distance(xu[ind,:],yu[ind,:]),(nz,1))
      Dr=np.tile(distance(xr[ind,:],yr[ind,:]),(NZ,1))
      Du=np.tile(distance(xu[ind,:],yu[ind,:]),(NZ,1))

      if useInd:
        distr[:,j,:]=dr;
        distu[:,j,:]=du;

      Temp[:,j,:]   = calc.griddata(Rdz*Dr,ZZr[:,j,:],TEMP[:,j,:],Rdz*dr,Zr[:,j,:],extrap=True,**interp_opts)
      Salt[:,j,:]   = calc.griddata(Rdz*Dr,ZZr[:,j,:],SALT[:,j,:],Rdz*dr,Zr[:,j,:],extrap=True,**interp_opts)
      if 0 and j%10==0:
        print(Dr.shape, ZZr[:,j,:].shape)
        import pylab as pl
        pl.figure(1)
        pl.clf()
        pl.pcolormesh(Dr,ZZr[:,j,:],SALT[:,j,:])
        pl.colorbar()
        clim=pl.gci().get_clim()
      
        pl.figure(2)
        pl.clf()
        pl.pcolormesh(dr,Zr[:,j,:],Salt[:,j,:])
        pl.clim(clim)
        pl.colorbar()
        raw_input()
      
      Uvel[:,j,:]   = calc.griddata(Rdz*Du,ZZu[:,j,:],U[:,j,:],   Rdz*du,Zu[:,j,:],extrap=True,**interp_opts)
      if j<Vvel.shape[1]:
        dv=np.tile(distance(xv[ind,:],yv[ind,:]),(nz,1))
        Dv=np.tile(distance(xv[ind,:],yv[ind,:]),(NZ,1))
        Vvel[:,j,:] = calc.griddata(Rdz*Dv,ZZv[:,j,:],V[:,j,:],   Rdz*dv,Zv[:,j,:],extrap=True,**interp_opts)
        if useInd:
          distv[:,j,:]=dv

      if np.any(np.isnan(Temp[:,j,:])): print('found nan in temp',j)
      if np.any(np.isnan(Salt[:,j,:])): print('found nan in salt',j)
      if np.any(np.isnan(Uvel[:,j,:])): print('found nan in u',j)
      if j<Vvel.shape[1] and np.any(np.isnan(Vvel[:,j,:])): print('found nan in v',j)


  elif ij=='i':
    for i in range(nx):
      if not quiet and (nx<10 or (nx>=10 and i%10==0)): print('  i=%3d of %3d' % (i,nx))
      ind=ij_ind[i]
      dr=np.tile(distance(xr[:,ind],yr[:,ind]),(nz,1))
      dv=np.tile(distance(xv[:,ind],yv[:,ind]),(nz,1))
      Dr=np.tile(distance(xr[:,ind],yr[:,ind]),(NZ,1))
      Dv=np.tile(distance(xv[:,ind],yv[:,ind]),(NZ,1))

      if useInd:
        distr[:,:,i]=dr;
        distv[:,:,i]=dv;

      Temp[:,:,i]   = calc.griddata(Rdz*Dr,ZZr[:,:,i],TEMP[:,:,i],Rdz*dr,Zr[:,:,i],extrap=True,**interp_opts)
      Salt[:,:,i]   = calc.griddata(Rdz*Dr,ZZr[:,:,i],SALT[:,:,i],Rdz*dr,Zr[:,:,i],extrap=True,**interp_opts)
      Vvel[:,:,i]   = calc.griddata(Rdz*Dv,ZZv[:,:,i],V[:,:,i],   Rdz*dv,Zv[:,:,i],extrap=True,**interp_opts)
      if i<Uvel.shape[2]:
        du=np.tile(distance(xu[:,ind],yu[:,ind]),(nz,1))
        Du=np.tile(distance(xu[:,ind],yu[:,ind]),(NZ,1))
        Uvel[:,:,i] = calc.griddata(Rdz*Du,ZZu[:,:,i],U[:,:,i],   Rdz*du,Zu[:,:,i],extrap=True,**interp_opts)
        if useInd:
          distu[:,:,i]=du


  # uv bar:
  if not quiet: print('calc uvbar')
  if useInd is False:
    ubar,vbar=rt.uvbar(Uvel,Vvel,sshr,hr,sparams)
  else: #>------------------------------------------------------------
    sshu=calc.griddata(dlon,dlat,data['ssh'],xu,yu,extrap=True,**interp_opts)
    sshv=calc.griddata(dlon,dlat,data['ssh'],xv,yv,extrap=True,**interp_opts)

    if ij=='j':
      sshu=sshu[ij_ind,:]
      sshv=sshv[ij_ind,:]
      hu  =hu[ij_ind,:]
      hv  =hv[ij_ind,:]
    elif ij=='i':
      sshu=sshu[:,ij_ind]
      sshv=sshv[:,ij_ind]
      hu  =hu[:,ij_ind]
      hv  =hv[:,ij_ind]

    ubar=rt.barotropic(Uvel,sshu,hu,sparams)
    vbar=rt.barotropic(Vvel,sshv,hv,sparams)
  # -----------------------------------------------------------------<


  Vars=cb.odict()
  Vars['temp'] = Temp
  Vars['salt'] = Salt
  Vars['u']    = Uvel
  Vars['v']    = Vvel
  Vars['zeta'] = sshr
  Vars['ubar'] = ubar
  Vars['vbar'] = vbar

  Vars['date']   = data['date']

  if not useInd is False: #>------------------------------------------
    Vars['depth']  = Zr
    Vars['depthu'] = Zu
    Vars['depthv'] = Zv

    Vars['dist']  = distr
    Vars['distu'] = distu
    Vars['distv'] = distv
  # -----------------------------------------------------------------<


  if retHorizAux: return Vars, horizAux
  else: return Vars


def data2romsbry(data,grd,sparams,**kargs):
  '''
  Returns prognostic variables at model boundaries.
  Data may be provided by load_data.

  Parameters
  ----------
  same as data2roms except the karg obc which defined the boundaries to
  be used, default is 'nsew' for north, south, east and west boundaries,
  ie, eta,xi [:,-1],[:,0],[-1,:] and [0,:]

  '''

  obc='nsew'
  if 'obc' in kargs.keys(): obc=kargs['obc']

  kargs['ij_ind']=0,-1
  res={}

  names='temp','salt','u','v','ubar','vbar','zeta','dist','distu',\
        'distv','depth','depthu','depthv'

  # load west and east data:
  if 'w' in obc or 'e' in obc:
    kargs['ij']='i'
    if not 'horizAux' in kargs.keys() or kargs['horizAux'] is False: kargs['horizAux']=True
    if kargs['horizAux'] is True:
      vars,horizAux=data2roms(data,grd,sparams,**kargs)
      kargs['horizAux']=horizAux
    else:
      vars=data2roms(data,grd,sparams,**kargs)

    if 'w' in obc:
      for n in names: res[n+'_west']=vars[n][...,0]

    if 'e' in obc:
      for n in names: res[n+'_east']=vars[n][...,1]


  # load south and north data:
  if 's' in obc or 'n' in obc:
    kargs['ij']='j'
    vars=data2roms(data,grd,sparams,**kargs)
    if 's' in obc:
      for n in names:
        if   vars[n].ndim==3: res[n+'_south']=vars[n][:,0,:]
        elif vars[n].ndim==2: res[n+'_south']=vars[n][0,:]

    if 'n' in obc:
      for n in names:
        if   vars[n].ndim==3: res[n+'_north']=vars[n][:,1,:]
        elif vars[n].ndim==2: res[n+'_north']=vars[n][1,:]

  res['date'] = vars['date']

  return res


def make_ini(ini,data,grid,sparams,quiet=0,**kargs):
  '''
  Create model initial file from data dict as obtained from data2roms

  Parameters
  ----------
  ini : file to create
  data : data as obtained with data2roms
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels
  quiet : output messages flag

  **kargs:
  type : global attribute
  title : global attribute
  date : date as datetime (or something convertible into datetime through
         date_tools.parse_date) or number (0)
  tunits : time units ('days since 1970-01-01')

  '''
  err=False

  q=gennc.GenIni(ini,grid,sparams,**kargs)
  q.create()
  q.fill(data,quiet=quiet)


def make_clm(clm,data,grid,sparams,quiet=0,**kargs):
  '''
  Create model climatology file from data dicts as obtained from data2roms

  Parameters
  ----------
  clm : file to create
  data : data as obtained with data2roms
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels
  quiet : output messages flag

  **kargs:
  type : global attribute
  title : global attribute ('ROMS Climatology file')
  cycle: time variables cycle_length (False, attribute not created)
  tunits : time units ('days since 1970-01-01')
  create: create netcdf file (True)
  '''

  create=True
  if 'create' in kargs.keys(): create=kargs['create']

  q=gennc.GenClm(clm,grid,sparams,**kargs)
  if create: q.create()

  q.fill(data,quiet=quiet)



def make_bry(bry,data,grid,sparams,quiet=0,**kargs):
  '''
  Create model boundary conditions file from source data files or from
  data dicts as obtained from data2romsbry

  Parameters
  ----------
  bry : file to create
  data : data as obtained with data2romsbry
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels
  quiet : output messages flag

  **kargs:
  type : global attribute
  title : global attribute
  cycle: time variables cycle_length (False, attribute not created)
  tunits : time units ('days since 1970-01-01')
  obc : open boundaries, the 4 boundaries are used by default ('nsew',
        for north, south, east and west)
  addxy : include distance x depth variables in the netcdf file (True)
  create: create netcdf file (True)

  '''
  create=True
  if 'create' in kargs.keys(): create=kargs['create']

  q=gennc.GenBry(bry,grid,sparams,**kargs)
  if create: q.create()

  q.fill(data,quiet=quiet)
