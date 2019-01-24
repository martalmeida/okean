'''
NetCDF tools

pythn 2x, 3x
'''

import numpy as np
from os.path import realpath
from functools import reduce
from .nc.pync4 import Pync
from . import calc
from .cookbook import isstr
lck=False

def ncopen(filename,perm='r',interface='auto',**kargs):
  if interface=='auto': return Pync(filename,perm=perm,**kargs)
  else: return Pync(filename,perm=perm,interface=interface,**kargs)


def __open(f,interface='auto',perm='r'):
  if isinstance(f,Pync):
    nc=f
    close=False
  else:
    nc=ncopen(f,interface=interface,perm=perm)
    close=True

  return nc,close


def vatt(fname,vname,attname=False,interface='auto'):
  nc,close=__open(fname,interface)

  if attname: res=nc.vars[vname].atts[attname].value
  else:       res=nc.vars[vname].atts
  if close: nc.close()
  return res


def fatt(fname,attname=False,interface='auto'):
  nc,close=__open(fname,interface)

  if attname: res=nc.atts[attname].value
  else: res=nc.atts
  if close: nc.close()
  return res


def vdim(fname,vname,dimname=False,interface='auto'):
  nc,close=__open(fname,interface)

  if dimname: res=nc.vars[vname].dims[dimname]
  else: res=nc.vars[vname].dims
  if close: nc.close()
  return res


def fdim(fname,dimname=False,interface='auto'):
  nc,close=__open(fname,interface)

  if dimname: res=nc.dims[dimname]
  else: res=nc.dims
  if close:  nc.close()
  return res


def var(fname,vname=False,perm='r',interface='auto'):
  nc,close=__open(fname,interface,perm)

  if vname: res=nc.vars[vname]
  else:     res=nc.vars
  return res##, nc # nc should be closed when not needed !!


def varnames(fname,interface='auto'):
  nc,close=__open(fname,interface)

  res=nc.varnames
  if close: nc.close()
  return res


def num2date(tnum,tunits,calendar='standard'):
  from netcdftime import utime

  #return utime(tunits,calendar).num2date(tnum)
  # for numpy 1.15.3 there are problems with masked arrays without mask!
  # for some reason, removing mask or dividing by 1 solves the problem!!!
  return utime(tunits,calendar).num2date(tnum/1)

def date2num(tdate,tunits,calendar='standard'):
  from netcdftime import utime
  return utime(tunits,calendar).date2num(tdate)


def nctime(filename,varname,interface='auto',**kargs):
  time=use(filename,varname,interface=interface,**kargs)
  if time is None: return
  units=vatt(filename,varname,'units')
  try: cal=vatt(filename,varname,'calendar')
  except: cal='standard'

  # for ROMS output files the calender attribute must be fixed!
  if cal=='gregorian_proleptic': cal='proleptic_gregorian'

  # dates like 2004-01-01T00:00:00Z are not supported by old varsions
  # of netcdftime (older than 0.9.2, don't really know which version).
  # units also cannot have multiple spaces... so:

  units=units.replace('T',' ').replace('Z',' ')
  units=' '.join(units.split())

  return num2date(time,units,cal)


def use(filename,varname,interface='auto',**kargs):
  nc,close=__open(filename,interface)

  if varname not in nc.varnames:
    print(':: variable %s not found'%varname)
    return

  v=nc.vars[varname]
  shape=v.shape()

  if not shape:
    try: return v[:][0] # may be needed for dtype='|S1'
    except: return v[:]
  else:
    size=reduce(lambda i,j: i*j,shape)
    if size==0: return


  d=v.dims
  dimStr=[':' for i in d]
  for k in kargs:
    dval=kargs[k]

    if k.isdigit(): k=list(d.keys())[int(k)] # allow use dim indice
    elif k.startswith('SEARCH'):  # allow *dimname or dimname*
      kk=k[len('SEARCH'):]
      for dn in nc.vars[varname].dimnames:
        if dn.find(kk)>=0:
          k=dn
          break
    elif k.endswith('SEARCH'):
      kk=k[:-len('SEARCH')]
      for dn in nc.vars[varname].dimnames:
        if dn.find(kk)>=0:
          k=dn
          break

    if k in d:
      i=list(d.keys()).index(k)
      if isstr(dval):
        dimStr[i]=dval
      elif isinstance(dval,int):
        if dval<0: dval=shape[i]+dval
        if  dval>shape[i]-1:
          print(':: max allowed '+k+' = '+str(shape[i]-1))
          return
        else:
          dimStr[i]=str(dval)
      elif calc.isiterable(dval):
        exec(k+'_val=dval')
        dimStr[i]=k+'_val'
      else:
        dimStr[i]=str(dval)

  cmd='v['+','.join(dimStr)+']'
  if lck:
    lck.acquire()
    try:
      res=eval(cmd)
    finally:
      lck.release()

  else:
    res=eval(cmd)


  if interface in ('pycdf','scientific'):
    # about missing value:
    miss=False
    if   '_FillValue'    in v.attnames: miss = v.atts['_FillValue']['value']
    elif 'missing_value' in v.attnames: miss = v.atts['missing_value']['value']
    maskMissing=kargs.get('maskMissing',True)
    if not miss is False and maskMissing and v.nctype()!='STRING':
      res=np.ma.masked_where(res==miss,res)

    ## ensure strings have no mask:
    #if v.nctype()=='STRING' and np.ma.isMA(res):
    #  res=np.array(res)


    # about scale and offset:
    if v.nctype()!='STRING':
      scale=1
      offset=0
      if 'scale_factor' in v.attnames: scale=v.atts['scale_factor']['value']
      if 'add_offset'   in v.attnames: offset=v.atts['add_offset']['value']
      if (scale,offset)!=(1,0): res = res*scale + offset

  if close: nc.close()

  if 1 in res.shape and res.ndim>1: res=np.squeeze(res)

  # mask nan
  maskNaN=kargs.get('maskNaN',True)
  if maskNaN and not res.dtype.type==np.string_ and not np.ma.isMA(res) and np.any(np.isnan(res)):
    res=np.ma.masked_where(np.isnan(res),res)

  return res


from .netcdf_new import ncshow as showfile
from .netcdf_new import dim
def showfile_old(filename,**kargs):
  interface = kargs.get('interface','auto')
  lmax      = kargs.get('lmax',False) # max len of attname of varname
  Lmax      = kargs.get('Lmax',70) # max len of file att or lon_name/units

  nc,close=__open(filename,interface)

  print("\n# Contents of the NetCDF file")
  #print('   '+realpath(nc.filename))
  print('   '+nc.filename)

  print("\n:: Global Attributes:")
  at=nc.atts
  if at:
    l1=reduce(max,[len(x) for x in at])

    try:
      l2=reduce(max,[len(str(x)) for x in [y.value for y in at.values()]])
    except:
      l2=reduce(max,[len(unicode(x)) for x in [y.value for y in at.values()]])

    if lmax: l1=min(lmax,l1)
    if Lmax: l2=min(Lmax,l2)

    format='   %-'+str(l1) + 's  %-'+str(l2)+'s'
    for k in at:
       try: v=str(at[k].value)
       except: v=unicode(at[k].value)
       if len(k)>l1: k=k[:l1-1]+'+'
       if len(v)>l2: v=v[:l2-1]+'+'
       print(format % (k,v))

  print("\n:: Dimensions:")
  di = nc.dims
  if di:
    l1=reduce(max,[len(x) for x in di])
    l2=reduce(max,[len(str(x)) for x in di.values()])
    format='   %-'+str(l1) + 's  %'+str(l2)+'d'
    for k in di:
       print(format % (k,di[k]))

  print('\n:: Variables:')
  va=nc.vars
  varnames = va.keys()
  if varnames:
    # find max lens
    # for vname:
    l1=reduce(max,[len(x) for x in varnames])
    # for long_name, units and shape:
    l2=14 # min len for long_name
    l3= 7 # min len for units
    l4= 7 # min len for str(shape)
    for v in varnames:
      vv=va[v]
      at=vv.atts

      if 'long_name' in at:
        longn=at['long_name'].value
        l2=max(l2,len(longn))
      if 'units' in at:
        units=at['units'].value
        l3=max(l3,len(units))

      l4=max(len(str(vv.shape())),l4)

    if Lmax:       
      l2=min(l2,Lmax)
      l3=min(l3,Lmax)
      l4=min(l4,Lmax)

    format='   %-'+str(l1)+'s | %-'+str(l2)+'s | %-'+str(l3)+'s | %-'+str(l4)+'s |'
    format1='   %-'+str(l1)+'s   %-'+str(l2)+'s   %-'+str(l3)+'s   %-'+str(l4)+'s'
    print(format1 % ('','long_name'.center(l2),'units'.center(l3),'shape'.center(l4)))
    for v in varnames:
      vv=va[v]
      at=vv.atts

      if 'long_name' in at: longn=at['long_name'].value
      else: longn=''
      if len(longn)>l2: longn=longn[:l2-1]+'+'

      if 'units' in at: units=at['units'].value
      else: units=''
      if len(units)>l3: units=units[:l3-1]+'+'

      shape=str(vv.shape())
      if len(shape)>l4: shape=shape[:l4-1]+'+'

      print(format % (v,longn,units,shape))

  if close: nc.close()


def showvar(filename,varname,**kargs):
  interface = kargs.get('interface','auto')
  lmax      = kargs.get('lmax',False)
  Lmax      = kargs.get('Lmax',70)

  nc,close=__open(filename,interface)
  vv=nc.vars[varname]
  at=vv.atts
  di=vv.dims

  print("\n:: Dimensions:")
  if di:
    l1=reduce(max,[len(x) for x in di])
    l2=reduce(max,[len(str(x)) for x in di.values()])
    format='   %-'+str(l1) + 's   %'+str(l2)+'d'
    for k in di:
       print(format % (k,di[k]))

  print("\n:: Attributes:")
  if at:
    l1=reduce(max,[len(x) for x in at])
    try:
      l2=reduce(max,[len(str(x)) for x in [y.value for y in at.values()]])
    except:
      l2=reduce(max,[len(unicode(x)) for x in [y.value for y in at.values()]])

    if lmax: l1=min(lmax,l1)
    if Lmax: l2=min(Lmax,l2)

    format='   %-'+str(l1) + 's   %-'+str(l2)+'s'
    for k in at:
       try: v=str(at[k].value)
       except: v=unicode(at[k].value)
       if len(k)>l1: k=k[:l1-1]+'+'
       if len(v)>l2: v=v[:l2-1]+'+'
       print(format % (k,v))

  print("\n:: Shape: %s" % str(vv.shape()))

  if vv.ndim()>0:
    size=reduce(lambda x, y: x*y, vv.shape())
    if size<1e4:
      try: print("\n:: Range: %s" % str(vv.range()))
      except: pass
  else:
    print("\n:: Value: %s" % str(vv[:]))

  if close: nc.close()


def show(filename,varname='',**kargs):
  if not varname: showfile(filename,**kargs)
  else: showvar(filename,varname,**kargs)


def dict2nc(fname,vars,dims={},atts={},vatts={},**kargs):
  '''
  Example:

  import numpy as np
  v1=np.empty((10,20,30),'float32') # netcdf float
  v2=np.empty((10,20),'float64')    # netcdf double

  # variables
  vars  = {'temp':v1,'ssh':v2}

  # dimensions:
  dims  = {'lon':10,'lat':20,'depth':30}

  # file attributes:
  atts  = {'history': 'some info...','time':12345}

  # variables attributes:
  vatts = {'temp':{'long_name':'temperature','units':'deg'}}

  # create the file:
  dict2nc('ncfile.nc',vars,dims,atts,vatts)

  # or simply: (will create dim_0, dim_1, ...)
  dict2nc('ncfile2.nc',vars)

  to overwrite a previous file, use karg perm='t'
  '''

  ncversion=kargs.get('ncversion','NETCDF4_CLASSIC')
  perm=kargs.get('perm','w')

  nc=Pync(fname,perm,version=ncversion)

  # add dims:
  for d in dims:
    nc.add_dim(d,dims[d])

  newDims={}
  # add vars:
  for vname in vars:
    # find var dims:
    try:
      sh=vars[vname].shape
    except:
      vars[vname]=np.asarray(vars[vname])
      sh=vars[vname].shape
    #if not sh: sh=1,
    dimnames=[]
    for s in sh:
      if s in dims.values(): dimnames+=[list(dims.keys())[list(dims.values()).index(s)]]
      else:
        # check inside new dims:
        if  s in newDims.values(): dimnames+=[list(newDims.keys())[list(newDims.values()).index(s)]]
        else:
          # or create a new dim:
          newDimName='dim_'+str(len(newDims.values()))
          newDims[newDimName]=s
          nc.add_dim(newDimName,s)
          dimnames+=[newDimName]
    try:
      fv=vatts[vname]['_FillValue']
    except: fv=None

    try:
      ncv=nc.add_var(vname,vars[vname].dtype,dimnames,ncver=ncversion,fill_value=fv)
    except:
      print('--> error creating variable %s of type %s (nc version = %s)'%(vname,vars[vname].dtype,str(ncversion)))
      return

    # add var atts:
    if vname in vatts:
      at=vatts[vname]
      for d in at:
        if d=='_FillValue': continue
        ncv.add_att(d,at[d])

  # add atts:
  for d in atts:
    nc.add_att(d,atts[d])

  # fill data:
  for vname in vars:
    nc.vars[vname][:]=vars[vname]

  nc.close()
