import datetime
import numpy as np
from okean import netcdf

def get_ij_inds(grd,**kargs):
  f=kargs.get('url','http://tds.hycom.org/thredds/dodsC/glb_analysis')
  vlon=kargs.get('vlon','Longitude')
  vlat=kargs.get('vlat','Latitude')
  lon=kargs.get('lon',False)
  lat=kargs.get('lat',False)
  lon_add=kargs.get('lon_add',-360)
  fsave=kargs.get('fsave','ijinds.pickle')

  if lon is False:
    lon=netcdf.use(f,vlon)
    if np.any(lon>360): lon=np.mod(lon,360)
    lon+=lon_add

  if lat is False:
    lat=netcdf.use(f,vlat)

  rlon=netcdf.use(grd,'lon_rho')
  rlat=netcdf.use(grd,'lat_rho')
  xlim=rlon.min(),rlon.max()
  ylim=rlat.min(),rlat.max()
  from okean import calc
  i1,i2,j1,j2=calc.ij_limits(lon,lat,xlim,ylim,margin=1)

  i1=i1-2
  i2=i2+2
  j1=j1-2
  j2=j2+2

  if fsave:
    np.asarray([ i1,i2,j1,j2]).dump(fsave)
    print 'saved %s'%fsave

  return np.asarray([ i1,i2,j1,j2])



def src(date,agg=False):
  if agg: return src_agg(date)
  else: return src_files(date)

def src_files(date,agg=False):
  if date>=datetime.datetime(2011,1,3):
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_90.9/%d'%date.year
  elif date>=datetime.datetime(2009,5,7):
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_90.8/%d'%date.year
  elif date>=datetime.datetime(2008,9,18):
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_90.6/%d'%date.year
  elif date>=datetime.datetime(2007,4,27):
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_90.3/%d'%date.year
  elif date>=datetime.datetime(2007,1,1):
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_90.2/%d'%date.year
  else:
    url0='http://dap.hycom.org:8080/opendap/nph-dods/datasets/GLBa0.08/expt_60.5/%d'%date.year

  res={}
  yd=(date-datetime.datetime(date.year,1,1)).days+1
  for v in ['temp','salt','u','v','ssh']:
    if   v=='temp': fname=url0+'/temp/archv.%d_%d_00_3zt.nc'%(date.year,yd)
    elif v=='salt': fname=url0+'/salt/archv.%d_%d_00_3zs.nc'%(date.year,yd)
    elif v=='u':    fname=url0+'/uvel/archv.%d_%d_00_3zu.nc'%(date.year,yd)
    elif v=='v':    fname=url0+'/vvel/archv.%d_%d_00_3zv.nc'%(date.year,yd)
    elif v=='ssh':  fname=url0+'/2d/archv.%d_%d_00_2d.nc'%(date.year,yd)

    res[v]=fname

  return res

def src_agg(date):
  '''may be very very slow !!'''

  if date>=datetime.datetime(2011,1,3):
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9'
  elif date>=datetime.datetime(2009,5,7):
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8/%d'%date.year
  elif date>=datetime.datetime(2008,9,18):
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6/%d'%date.year
  elif date>=datetime.datetime(2007,4,27):
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.3/%d'%date.year
  elif date>=datetime.datetime(2007,1,1):
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.2/%d'%date.year
  else:
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5/%d'%date.year


  res={}
  for v in ['temp','salt','u','v','ssh']:
    res[v]=url

  return res
