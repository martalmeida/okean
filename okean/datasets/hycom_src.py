import datetime


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
