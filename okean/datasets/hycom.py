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
