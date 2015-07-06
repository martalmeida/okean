from nc.pync4 import Pync
import calc
from string import join as sjoin
import numpy as np
from os.path import realpath


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
  return utime(tunits,calendar).num2date(tnum)

def date2num(tdate,tunits,calendar='standard'):
  from netcdftime import utime
  return utime(tunits,calendar).date2num(tdate)


def nctime(filename,varname,interface='auto',**kargs):
  time=use(filename,varname,interface=interface,**kargs)
  units=vatt(filename,varname,'units')

  # dates like 2004-01-01T00:00:00Z not supported by old varsions of
  # netcdftime (older than 0.9.2, dont really know which version).
  # units also cannot have multiple spaces... so:

  units=units.replace('T',' ').replace('Z',' ')
  units=' '.join(units.split())
  return num2date(time,units)


def use(filename,varname,interface='auto',**kargs):
  nc,close=__open(filename,interface)

  if varname not in nc.varnames: return

  v=nc.vars[varname]
  shape=v.shape()

  if not shape:
    try: return v[:][0] # may be needed for dtype='|S1'
    except: return v[:]


  d=v.dims
  dimStr=[':' for i in d.keys()]
  for k in kargs.keys():
    dval=kargs[k]

    if k.isdigit(): k=d.keys()[int(k)] # allow use dim indice
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

    if k in d.keys():
      i=d.index(k)
      if isinstance(dval,basestring):
        dimStr[i]=dval
      elif isinstance(dval,int):
        if dval<0: dval=shape[i]+dval
        if  dval>shape[i]-1:
          print ':: max allowed '+k+' = '+str(shape[i]-1)
          return
        else:
          dimStr[i]=str(dval)
      elif calc.isiterable(dval):
        exec(k+'_val=dval')
        dimStr[i]=k+'_val'
      else:
        dimStr[i]=str(dval)

  cmd='res=v['+sjoin(dimStr,',')+']'
  exec(cmd)

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


def show(filename,**kargs):
  lmax      = False
  Lmax      = 100
  interface = 'auto'

  if 'lmax'      in kargs.keys(): lmax      = kargs['lmax']
  if 'Lmax'      in kargs.keys(): Lmax      = kargs['Lmax']
  if 'interface' in kargs.keys(): interface = kargs['interface']

  nc=ncopen(filename,interface=interface)

  print "\n# Contents of the NetCDF file"
  print '   '+realpath(filename)

  print "\n:: Global Attributes:"
  at=nc.atts
  if at:
    l1=reduce(max,[len(x) for x in at.keys()])
    try:
      l2=reduce(max,[len(str(x)) for x in [y.value for y in at.values()]])
    except:
      l2=reduce(max,[len(unicode(x)) for x in [y.value for y in at.values()]])

    if lmax: l1=min(lmax,l1)
    if Lmax: l2=min(Lmax,l2)

    format='   %-'+str(l1) + 's   %-'+str(l2)+'s'
    for k in at.keys():
       try: v=str(at[k].value)
       except: v=unicode(at[k].value)
       if len(k)>l1: k=k[:l1-1]+'+'
       if len(v)>l2: v=v[:l2-1]+'+'
       print format % (k,v)

  print "\n:: Dimensions:"
  di = nc.dims
  if di:
    l1=reduce(max,[len(x) for x in di.keys()])
    l2=reduce(max,[len(str(x)) for x in di.values()])
    format='   %-'+str(l1) + 's   %'+str(l2)+'d'
    for k in di.keys():
       print format % (k,di[k])

  print "\n:: Variables:"
  va=nc.vars
  varnames = va.keys()
  if varnames:
    l1=reduce(max,[len(x) for x in varnames])
    l2=40
    if Lmax: l2=min(Lmax,l2)
    format='   %-'+str(l1)+'s  %-'+str(l2)+'s  %-15.15s  %-20s'
    for v in varnames:
      vv=va[v]
      at=vv.atts

      if 'long_name' in at.keys(): longn=at['long_name'].value
      else: longn=''
      if len(longn)>l2: longn=longn[:l2-1]+'+'

      if 'units' in at.keys(): units=at['units'].value
      else: units=''

      shape=vv.shape()
      print format % (v,longn,units,str(shape))

  nc.close()


def showvar(filename,varname,**kargs):
  lmax      = False
  Lmax      = 100
  interface = 'auto'

  if 'lmax'      in kargs.keys(): lmax      = kargs['lmax']
  if 'Lmax'      in kargs.keys(): Lmax      = kargs['Lmax']
  if 'interface' in kargs.keys(): interface = kargs['interface']

  nc=ncopen(filename,interface=interface)
  vv=nc.vars[varname]
  at=vv.atts
  di=vv.dims

  print "\n:: Dimensions:"
  if di:
    l1=reduce(max,[len(x) for x in di.keys()])
    l2=reduce(max,[len(str(x)) for x in di.values()])
    format='   %-'+str(l1) + 's   %'+str(l2)+'d'
    for k in di.keys():
       print format % (k,di[k])

  print "\n:: Attributes:"
  if at:
    l1=reduce(max,[len(x) for x in at.keys()])
    try:
      l2=reduce(max,[len(str(x)) for x in [y.value for y in at.values()]])
    except:
      l2=reduce(max,[len(unicode(x)) for x in [y.value for y in at.values()]])

    if lmax: l1=min(lmax,l1)
    if Lmax: l2=min(Lmax,l2)

    format='   %-'+str(l1) + 's   %-'+str(l2)+'s'
    for k in at.keys():
       try: v=str(at[k].value)
       except: v=unicode(at[k].value)
       if len(k)>l1: k=k[:l1-1]+'+'
       if len(v)>l2: v=v[:l2-1]+'+'
       print format % (k,v)

  print "\n:: Shape: %s" % str(vv.shape())

  try: print "\n:: Range: %s" % str(vv.range())
  except: pass

  nc.close()


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

  perm='w'
  ncversion=3
  if 'perm'      in kargs.keys(): perm      = kargs['perm']
  if 'ncversion' in kargs.keys(): ncversion = kargs['ncversion']

  nc=Pync(fname,perm,version=ncversion)

  # add dims:
  for d in dims.keys():
    nc.add_dim(d,dims[d])

  newDims={}
  # add vars:
  for vname in vars.keys():
    # find var dims:
    sh=vars[vname].shape
    dimnames=[]
    for s in sh:
      if s in dims.values(): dimnames+=[dims.keys()[dims.values().index(s)]]
      else:
        # check inside new dims:
        if  s in newDims.values(): dimnames+=[newDims.keys()[newDims.values().index(s)]]
        else:
          # or create a new dim:
          newDimName='dim_'+str(len(newDims.values()))
          newDims[newDimName]=s
          nc.add_dim(newDimName,s)
          dimnames+=[newDimName]
    try:
      fv=vatts[vname]['_FillValue']
    except: fv=None

    #print 'adding ',vname, vars[vname].dtype
    ncv=nc.add_var(vname,vars[vname].dtype,dimnames,ncver=ncversion,fill_value=fv)

    # add var atts:
    if vname in vatts.keys():
      at=vatts[vname]
      for d in at.keys():
        if d=='_FillValue': continue
        ncv.add_att(d,at[d])

  # add atts:
  for d in atts.keys():
    nc.add_att(d,atts[d])

  # fill data:
  for vname in vars.keys():
    nc.vars[vname][:]=vars[vname]

  nc.close()
