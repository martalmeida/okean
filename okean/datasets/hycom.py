import datetime

def get_ij_inds(grd,**kargs):
  f=kargs.get('url','http://tds.hycom.org/thredds/dodsC/glb_analysis')
  vlon=kargs.get('vlon','Longitude')
  vlat=kargs.get('vlat','Latitude')
  lon=kargs.get('lon',False)
  lat=kargs.get('lat',False)
  lon_add=kargs.get('lon_add',-360)

  if lon is False:
    lon=netcdf.use(f,vlon)
    if lon_add: lon=np.mod(lon,360)+lon_add

  if lat is False:
    lat=netcdf.use(f,vlat)

  g=roms.Grid(grd)
  xlim=g.lon.min(),g.lon.max()
  ylim=g.lat.min(),g.lat.max()
  from okean import calc
  i1,i2,j1,j2=calc.ij_limits(lon,lat,xlim,ylim,margin=1)

  i1=i1-2
  i2=i2+2
  j1=j1-2
  j2=j2+2

  np.array([ i1,i2,j1,j2]).dump(fsave)
  print '%d %d %d %d  saved in %s'%(i1,i2,j1,j2,fsave)



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

''' OLD HYCOM.py. Probably not needed anymore!!! REMOVE some day. (jan2014)
import numpy as np
import os
from okean import dateu, netcdf


def hycom_vname(roms_vname):
  if   roms_vname=='temp': return 'temperature'
  elif roms_vname=='salt': return 'salinity'
  elif roms_vname=='zeta': return 'ssh'
  else: return roms_vname # u, v


def gen_url(date,varname,dataset='GOM'):

  names={}
  names['temperature']='temp'
  names['salinity']='salt'
  names['u']='uvel'
  names['v']='vvel'

  if varname in names.keys(): varname=names[varname]

  date=dateu.parse_date(date)
  y=date.year
  if dataset=='GOM':
    date=date.strftime('%Y%m%d')
    if date<='20120406':
      return 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1'
    else:
      return 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_31.0'


    url0='http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1/'
    date=date.strftime('%Y%m%d')

    if date>='20120101':
      url='http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2012'
    elif date>='20110101':
      return 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2011'

  elif dataset=='GLOBAL':
    return 'http://tds.hycom.org/thredds/dodsC/glb_analysis'


def get_lonlat(date,vname,dataset='GOM'):
  url=gen_url(date,vname,dataset)
  lon=netcdf.use(url,'Longitude')
  lat=netcdf.use(url,'Latitude')

  if dataset=='GLOBAL': lon=np.mod(lon,360)-360

  return lon,lat


def get_ij(lon,lat,lons,lats):
  if 0:
    # accurate but slow!! better for non regular grids..
    from pyp.tools import num_tools as nt
    i1,i2,j1,j2=nt.ij_limits(lon,lat,lons,lats)
  else:
    if lon.ndim==1:
      i1=np.where(lon<lons[0])[0][-1]
      i2=np.where(lon>lons[1])[0][0]

      j1=np.where(lat<lats[0])[0][-1]
      j2=np.where(lat>lats[1])[0][0]

    elif lon.ndim==2:
      d=(lon-lons[0])**2+(lat-lats[0])**2
      mm=d==d.min()
      try:
        j1,i1=np.where(d==d.min())
      except:
        j1,i1=np.ma.where(d==d.min())

      i1=i1-1
      j1=j1-1

      d=(lon-lons[1])**2+(lat-lats[1])**2
      try:
        j2,i2=np.where(d==d.min())
      except:
        j2,i2=np.ma.where(d==d.min())

      i2=i2+1
      j2=j2+1

  return i1,i2,j1,j2
'''
