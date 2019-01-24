import netCDF4
import numpy as np

isstr=lambda s: isinstance(s,(''.__class__,u''.__class__))


def ncopen(f,mode='r'):
  if isinstance(f,netCDF4.Dataset) or isinstance(f,netCDF4.MFDataset):
    return f

  if not isstr(f) or any([i in f for i in '*?']):
    nc=netCDF4.MFDataset(f)
  else:
    nc=netCDF4.Dataset(f,mode)

  return nc

def __open(f,mode='r'):
  if isinstance(f,netCDF4.Dataset):
    nc=f
    close=False
  else:
    nc=ncopen(f,mode)
    close=True

  return nc,close


def dim(f,dimname=False):
  nc,close=__open(f)

  if dimname: res=nc.dimensions[dimname].size
  else: res=dict([(i,nc.dimensions[i].size) for i in nc.dimensions])
  if close:  nc.close()
  return res


def ncshow(f,**kargs):
  from functools import reduce

  lmax      = kargs.get('lmax',False) # max len of attname of varname
  Lmax      = kargs.get('Lmax',70) # max len of file att or lon_name/units

  nc=ncopen(f)

  print('\n# Contents of the NetCDF file')
  try:
    print('   '+f)
  except: 
    print('   '+str(f.__class__))

  print('\n:: Global Attributes:')
  atn=list(nc.ncattrs())
  atv=[getattr(nc,i) for i in atn]
  if atn:
    l1=reduce(max,[len(x) for x in atn])
    try:
      l2=reduce(max,[len(str(x)) for x in atv])
    except:
      l2=reduce(max,[len(unicode(x)) for x in atv])

    if lmax: l1=min(lmax,l1)
    if Lmax: l2=min(Lmax,l2)

    format='   %-'+str(l1) + 's  %-'+str(l2)+'s'
    for i,k in enumerate(atn):
       try: v=str(atv[i])
       except: v=unicode(at[i])
       if len(k)>l1: k=k[:l1-1]+'+'
       if len(v)>l2: v=v[:l2-1]+'+'
       print(format % (k,v))

  print('\n:: Dimensions:')
  din = list(nc.dimensions)
  div=[]
  unlim=[]
  for i in din:
    if hasattr(nc.dimensions[i],'size'):
      div+=[nc.dimensions[i].size]
    else: div+=[nc.dimensions[i].dimtotlen]
    if nc.dimensions[i].isunlimited(): unlim+=[1]
    else: unlim+=[0]

  if din:
    l1=reduce(max,[len(x) for x in din])
    l2=reduce(max,[len(str(x)) for x in div])
    format='   %-'+str(l1) + 's  %'+str(l2)+'d'
    for i,k in enumerate(din):
       if unlim[i]:
         print(format % (k,div[i]) + ' (unlimited)')
       else:
         print(format % (k,div[i]))

  print('\n:: Variables:')
  varnames = list(nc.variables)
  if varnames:
    # find max len
    # for vname:
    l1=reduce(max,[len(x) for x in varnames])
    # for long_name, units and shape:
    l2=14 # min len for long_name
    l3= 7 # min len for units
    l4= 7 # min len for str(shape)
    anyUnits=False
    for v in varnames:
      atn=list(nc.variables[v].ncattrs())
      atv=[getattr(nc.variables[v],i) for i in atn]

      if 'long_name' in atn:
        longn=atv[atn.index('long_name')]
        l2=max(l2,len(longn))
      if 'units' in atn:
        units=atv[atn.index('units')]
        l3=max(l3,len(units))
        anyUnits=True

      l4=max(len(str(nc.variables[v].shape)),l4)

    if Lmax:
      l2=min(l2,Lmax)
      l3=min(l3,Lmax)
      l4=min(l4,Lmax)

    if anyUnits:
      format ='   %-'+str(l1)+'s | %-'+str(l2)+'s | %-'+str(l3)+'s | %-'+str(l4)+'s |'
      format1='   %-'+str(l1)+'s   %-'+str(l2)+'s   %-'+str(l3)+'s   %-'+str(l4)+'s'
      print(format1 % ('','long_name'.center(l2),'units'.center(l3),'shape'.center(l4)))
    else:
      format ='   %-'+str(l1)+'s | %-'+str(l2)+'s | %-'+str(l4)+'s |'
      format1='   %-'+str(l1)+'s   %-'+str(l2)+'s   %-'+str(l4)+'s'
      print(format1 % ('','long_name'.center(l2),'shape'.center(l4)))

    for v in varnames:
      atn=list(nc.variables[v].ncattrs())
      atv=[getattr(nc.variables[v],i) for i in atn]

      if 'long_name' in atn: longn=atv[atn.index('long_name')]
      else: longn=''
      if len(longn)>l2: longn=longn[:l2-1]+'+'

      if 'units' in atn: units=atv[atn.index('units')]
      else: units=''
      if len(units)>l3: units=units[:l3-1]+'+'

      shape=str(nc.variables[v].shape)
      if len(shape)>l4: shape=shape[:l4-1]+'+'

      if anyUnits:
        print(format % (v,longn,units,shape))
      else:
        print(format % (v,longn,shape))

  nc.close()
