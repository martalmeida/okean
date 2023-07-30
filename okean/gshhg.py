'''
data extraction from the Global Self-consistent Hierarchical
High-resolution Geography, GSHHG

https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/

Uses the shapefiles from gshhg-shp-<VERSION>.zip
Tested with gshhg-shp-2.3.7.zip

About:
  Currently Cartopy can use several Natural Earth features and coastlines from GSHHS
  (https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.feature.GSHHSFeature.html)

  However, Natural Earth resolution is very poor in many regions, where GSHHG seems much better.
  GSHHG dataset inclues GSHHS (coastlines) and WDBII (national boundaries and rivers).

  The purpose of this model is to extract national boundaries and rivers (as well as coastlines)
  from GSHHG to be used with Cartopy.

  GSHHS and WDBII are availabe in Basemap and the binary datasets are dowloaded during installation,
  except for higest resolution. This module, however, uses the shapefiles which must be downloaded
  and the data path provided to the gshhg class.
'''

import shapefile
import numpy as np
import os

def clip(lims,lon,lat):
  '''
  Clip vectorial data inside to rectangle
  lims: x0,x1,y0,y1
  '''

  X0,X1,Y0,Y1=lims
  _ll_lon=X0
  ll_lat=Y0
  _ur_lon=X1
  ur_lat=Y1

  # Clip data to window
  lon = np.clip(lon, _ll_lon, _ur_lon)
  lat = np.clip(lat, ll_lat, ur_lat)

  # shrink corners to a single point
  b = np.array([0.5, 0.5])
  idx = ( ((lon==_ll_lon)&(lat==ll_lat))  # South-west corner
        | ((lon==_ur_lon)&(lat==ll_lat))  # South-east corner
        | ((lon==_ll_lon)&(lat==ur_lat))  # North-west corner
        | ((lon==_ur_lon)&(lat==ur_lat)) )# North-east corner
  idx=np.convolve(idx, b, mode='same')
  lon=lon[idx<1]; lat=lat[idx<1]

  # Only keep first and last points in series along each boundary
  b = np.array([0.25, 0.5, 0.25])
  if lon.size>=b.size:
    idx = (lon==_ll_lon)
    idx = np.convolve(idx, b, mode='same')
    lon=lon[idx<1]; lat=lat[idx<1]

  if lon.size>=b.size:
    idx = (lon==_ur_lon)
    idx = np.convolve(idx, b, mode='same')
    lon=lon[idx<1]; lat=lat[idx<1]

  if lon.size>=b.size:
    idx = (lat==ll_lat)
    idx = np.convolve(idx, b, mode='same')
    lon=lon[idx<1]; lat=lat[idx<1]

  if lon.size>=b.size:
    idx = (lat==ur_lat)
    idx = np.convolve(idx, b, mode='same')
    lon=lon[idx<1]; lat=lat[idx<1]

  return lon,lat


def _get_common(f,lims,min_area=0,clip_data=0):
  '''
  Extract data from GSHHG shapefile f
  lims: x0,x1,y0,y1
  min_area: smallest area of polygons if data to return has area
  '''

  X0,X1,Y0,Y1=lims
  a=shapefile.Reader(f)
  n=len(a.shapes())
  shapes=a.shapes()
  n=len(shapes)

  records = a.records()
  out=[]
  #a.fields
  for i in range(n):
    #s=a.shape(i)
    # or
    s=shapes[i]
    x0,y0,x1,y1=s.bbox
    if x0>X1 or x1<X0 or y0>Y1 or y1<Y0: continue

    x,y=np.asarray(s.points).T
    #sType=s.shapeType # POLYGON = 5
    sTypeN=s.shapeTypeName # POLYGON

    r=records[i]
    if 'area' in r.as_dict() and min_area:
      if not r.area<min_area:
        out+=[(x,y)]
    else:
      # lines could have parts. Not the gshhg-shp-2.3.7.zip dataset, but just in case:
      if len(s.parts)>1:
        for j in range(len(s.parts)):
          j0=s.parts[j]
          try: j1=s.parts[j+1]
          except: j1=-1
          out+=[(x[j0:j1],y[j0:j1])]
      else:
        out+=[(x,y)]

    #level=r.level

  if clip_data:
    out=[clip(lims,*i) for i in out]

  return out


class gshhg():
  def __init__(self,lims,data_path):
    '''
    lims: X0,X1,Y0,Y1
    data_path: path to the folders  GSHHS_shp and WDBII_shp

    All data sets come in 5 different resolutions:
      f : Full resolution.  These contain the maximum resolution
          of this data and has not been decimated.
      h : High resolution.  The Douglas-Peucker line reduction was
          used to reduce data size by ~80% relative to full.
      i : Intermediate resolution.  The Douglas-Peucker line reduction was
          used to reduce data size by ~80% relative to high.
      l : Low resolution.  The Douglas-Peucker line reduction was
          used to reduce data size by ~80% relative to intermediate.
      c : Crude resolution.  The Douglas-Peucker line reduction was
          used to reduce data size by ~80% relative to low.

    Usage example:
      import pylab as pl
      lims=-57,-30,-8,10
      X0,X1,Y0,Y1=lims
      pl.plot([X0,X1,X1,X0,X0],[Y0,Y0,Y1,Y1,Y0])

      a=gshhg(lims,'./gshhs')
      out=a.get_coast(res='i',lev=1,min_area=0,clip_data=1)
      [pl.plot(*i,'k',lw=2) for i in out]

      out=a.get_boundary(res='i',lev=1,clip_data=0)
      [pl.plot(*i,'g',lw=2) for i in out]

      out=a.get_river(res='i',lev=1,clip_data=1)
      [pl.plot(*i,'b',lw=2) for i in out]

    '''

    self.lims=lims
    self.data_path=data_path

  def _file(self,data_type,res,lev):
    if data_type == 'coast':
      f=os.path.join(self.data_path,'GSHHS_shp','%s'%res,'GSHHS_%s_L%d'%(res,lev))
    elif data_type=='boundary':
      f=os.path.join(self.data_path,'WDBII_shp','%s'%res,'WDBII_border_%s_L%d'%(res,lev))
    elif data_type=='river':
      f=os.path.join(self.data_path,'WDBII_shp','%s'%res,'WDBII_river_%s_L%02d'%(res,lev))

    return f

  def get_coast(self,res,lev,min_area=0,clip_data=True):
    '''
    Coastline, islands and lakes

    res: resolution (c,l,i,h,f)
    lev: data level:

      Level 1: Continental land masses and ocean islands, except Antarctica.
      Level 2: Lakes
      Level 3: Islands in lakes
      Level 4: Ponds in islands within lakes
      Level 5: Antarctica based on ice front boundary.
      Level 6: Antarctica based on grounding line boundary.

    min_area: smallest area of the polygons to return
    '''
    f=self._file('coast',res,lev)
    return _get_common(f,self.lims,min_area,clip_data)

  def get_boundary(self,res,lev,clip_data=True):
    '''
    National and international borders

    res: resolution (c,l,i,h,f)
    lev: data level:

      Level 1: National boundaries.
      Level 2: Internal (state) boundaries for the 8 largest countries only.
      Level 3: Maritime boundaries.
    '''
    f=self._file('boundary',res,lev)
    return _get_common(f,self.lims,clip_data=clip_data)

  def get_river(self,res,lev,clip_data=True):
    '''
    Rivers

    lev: data level:

      Level  1: Double-lined rivers (river-lakes).
      Level  2: Permanent major rivers.
      Level  3: Additional major rivers.
      Level  4: Additional rivers.
      Level  5: Minor rivers.
      Level  6: Intermittent rivers - major.
      Level  7: Intermittent rivers - additional.
      Level  8: Intermittent rivers - minor.
      Level  9: Major canals.
      Level 10: Minor canals.
      Level 11: Irrigation canals.
    '''
    f=self._file('river',res,lev)
    return _get_common(f,self.lims,clip_data=clip_data)
