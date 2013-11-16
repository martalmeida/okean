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


def show(fname):
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

    print '%10s %30s %40s %10s %20s' % (abb.ljust(10),name.ljust(30)[:30],level.ljust(40)[:40],vlevel.ljust(10)[:10],units.ljust(20)[:20])


def findvar(fname,name,*args):
  '''
  Locates variable inside grib file
  Inputs:
    name, part of variable parameter
    args, any string inside variable string representation

  Example:
    findvar('file.grib2','temperatue','10 m')
  '''

  G=Grib2Decode(fname)
  out=[]
  for g in G:
    if isPygrib:
      if g['parameterName'].lower().find(name.lower())>=0:
        if all([str(g).lower().find(i.lower())>=0 for i in args]): out+=[g]
    else:
      if g.parameter.lower().find(name.lower())>=0:
        if all([str(g).lower().find(i.lower())>=0 for i in args]): out+=[g]

  return out


def cross_lon0(x,y,v):
  '''
  cross longitude 0
  used by getvar
  '''
  i=np.where(x[0]==180.)[0]
  X=x.copy()
  Y=y.copy()
  V=v.copy()
  X[:,i-1:]=x[:,:i+1]
  X[:,:i-1]=x[:,i+1:]
  Y[:,i-1:]=y[:,:i+1]
  Y[:,:i-1]=y[:,i+1:]
  V[:,i-1:]=v[:,:i+1]
  V[:,:i-1]=v[:,i+1:]
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


def getvar(fname,name,quiet=1,**kargs):
  '''
  Returns data from GFS grib file

  Inputs:
    Name, grib filename
    quiet, print info flag
    kargs:
      lons, x limits
      lats, y limits
      tags: string or list of string used by findvar
      neglon, if true west is negative (True)

  Example:
      gfs_getvar('file.grib2','temperature',tags='2 m')
      gfs_getvar('file.grib2','temperature',tags='2 m',lons=(-60,-30),lats=(-50,0))
  '''

  if not isinstance(fname,basestring): fname=fname['name']
  lons=False
  lats=False
  neglon=True
  tags=[]
  for k in kargs.keys():
    if k=='lons':   lons   = kargs[k]
    if k=='lats':   lats   = kargs[k]
    if k=='neglon': neglon = kargs[k]
    if k=='tags':   tags   = kargs[k]

  if isinstance(tags,str): tags=[tags]
  var=findvar(fname,name,*tags)

  if len(var)==1:
    var=var[0]

    if isPygrib:
      lat,lon=var.latlons()
      data=var.values
    else:
      lat,lon=var.grid()
      data=var.data()

    if neglon:
      lon=np.where(lon>180.,lon-360.,lon)
      lon,lat,data=cross_lon0(lon,lat,data)

    if lons or lats:
      lon,lat,data=extract_region(lon,lat,data,lons,lats)

    return lon,lat,data

  elif len(var)>1:
    if not quiet: print 'more than one var found !! found', len(var)
    return False,False,False

  else:
    if not quiet: print "not found"
    return False,False,False
