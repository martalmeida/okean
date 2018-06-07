'''
grib utils
Tools to extract info and data form grib files

Martinho MA, UFBA, Salvador, Brasil, 2009
TAMU, Texas, USA, 2011
'''

import numpy as np
try:
  from grib2 import Grib2Decode
  isPygrib=False
except:
  import pygrib
  Grib2Decode=pygrib.open
  isPygrib=True

from okean import calc


def show_all(fname):
  '''
  Lists variables inside grib file
  '''

  G=Grib2Decode(fname)
  tmp=[]
  for g in G:
    if isPygrib:
      abb    = g['shortName']
      name   = g['parameterName']
      level  = ''
      vlevel = str(g['level'])
      units  = g['parameterUnits']
    else:
      abb    = g.parameter_abbrev
      name   = g.parameter
      level  = g.vertical_level_descriptor
      vlevel = g.vertical_level
      units  = g.parameter_units

#    print '%10s %30s %40s %10s %20s' % (abb.ljust(10),name.ljust(30)[:30],level.ljust(40)[:40],vlevel.ljust(10)[:10],units.ljust(20)[:20])
    print('%-10s %-20s %-20s %-10s %-10s' % (abb,name[:20],level[:20],vlevel[:10],units[:10]))


def show(fname):
  # organize variables with same short name:
  G=Grib2Decode(fname)
  from collections import OrderedDict
  vars=OrderedDict()
  varsn={}
  for g in G:
    abb    = g['shortName']
    if abb in vars.keys():
      varsn[abb]+=1
    else:
      vars[abb]=g['parameterName'],str(g['level']),g['parameterUnits']
      varsn[abb]=1


  print('%-10s %-25s %5s   %-10s' % ('','name','n lev','units'))
  for abb in vars.keys():
    name,vlevel,units=vars[abb]
    vlevel=str(varsn[abb])
    print('%-10s %-25s %5s   %-10s' % (abb,name[:25],vlevel[:10],units[:10]))



def findvar(fname,name,*args):
  '''
  Locates variable inside grib file
  Inputs:
    name, part of variable short name (parameter abbv)
    args, any string inside variable string representation

  Example:
    findvar('file.grib2','temperatue','10 m')
  '''

  G=Grib2Decode(fname)
  out=[]
  for g in G:
    if isPygrib:
      #if g['parameterName'].lower().find(name.lower())>=0:
      if g['shortName'].lower()==name:
        if all([str(g).lower().find(i.lower())>=0 for i in args]): out+=[g]
    else:
      #if g.parameter.lower().find(name.lower())>=0:
      if g.parameter_abbrev.lower()==name:
        if all([str(g).lower().find(i.lower())>=0 for i in args]): out+=[g]

  return out


def cross_lon0(x,y,v):
  '''
  cross longitude 0
  used by getvar
  '''
  #i=np.where(x[0]==180.)[0]
  i=np.where(x[0]>180.)[0]
  L=len(i)
#  print x[0].min(),x[0].max()
  i=i[0]
  X=x.copy()
  Y=y.copy()
  V=v.copy()
  #X[:,i-1:]=x[:,:i+1]
  #X[:,:i-1]=x[:,i+1:]
  #Y[:,i-1:]=y[:,:i+1]
  #Y[:,:i-1]=y[:,i+1:]
  #V[:,i-1:]=v[:,:i+1]
  #V[:,:i-1]=v[:,i+1:]
#  print X.shape, L
  X[:,:L]=x[:,i:]
  X[:,L:]=x[:,:i]
  Y[:,:L]=y[:,i:]
  Y[:,L:]=y[:,:i]
  V[:,:L]=v[:,i:]
  V[:,L:]=v[:,:i]
  X=np.where(X>180.,X-360.,X)
  return X,Y,V


def extract_region(lon,lat,data,lons,lats):
  '''
  extract data inside region lonsxlats (xlim x ylim)
  used by getvar
  '''
  if not lons: lons=lon.min(),lon.max()
  if not lats: lats=lat.min(),lat.max()
  i1,i2,j1,j2=calc.ij_limits(lon,lat,lons,lats,margin=1)
  return lon[j1:j2,i1:i2],lat[j1:j2,i1:i2],data[j1:j2,i1:i2]


def getvar(*args,**kargs):
  '''
  Returns data from grib file

  Inputs:
    Name, grib filename
    kargs:
      lons, x limits
      lats, y limits
      tags: string or list of string used by findvar
      neglon, if true west is negative (True)

  Example:
      getvar('file.grib2','temperature',tags='2 m')
      getvar('file.grib2','temperature',tags='2 m',lons=(-60,-30),lats=(-50,0))
      or:
      var=findvar('file.grib2','temperature','2 m')
      getvar(var)
  '''


  lons   = kargs.get('lons',False)
  lats   = kargs.get('lats',False)
  neglon = kargs.get('neglon',True)
  tags   = kargs.get('tags',[])
  quiet  = kargs.get('quiet',False)

  if len(args)==2:
    fname,name=args
#    if not isinstance(fname,basestring): fname=fname['name'] # what for is this?
    if isinstance(tags,str): tags=[tags]
    var=findvar(fname,name,*tags)
  elif len(args)==1:
    var=args[0]

  try: len(var)
  except: var=[var] # an iterable is expected, like the output of findvar

  if len(var)>1:
    if not quiet: print('more than one var found !! found', len(var))
    return False,False,False
  elif len(var)==0:
    if not quiet: print("no variable found")
    return False,False,False
  else: var=var[0]

  if isPygrib:
    lat,lon=var.latlons()
    data=var.values
  else:
    lat,lon=var.grid()
    data=var.data()

  if neglon:
    #lon=np.where(lon>180.,lon-360.,lon)
    lon,lat,data=cross_lon0(lon,lat,data)

  if lons or lats:
    lon,lat,data=extract_region(lon,lat,data,lons,lats)

  return lon,lat,data


