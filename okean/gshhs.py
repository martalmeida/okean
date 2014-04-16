"""
Module for reading gshhs binary coastline formats
Martinho MA
original from Robert Hetland, Texas A&M

data available at http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/

The env variable GSHHS_DATA must be set (folder where files like
gshhs_c.b, etc, should be placed)
"""

import struct
from numpy import *
import numpy as np
import pylab as pl
import os
from matplotlib.ticker import FuncFormatter

def deg2str(deg, dir='E', fmt="%3.1f"):
    min = 60 * (deg - floor(deg))
    deg = floor(deg)
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        if dir=='E':
            dir='W'
        if dir=='N':
            dir='S'
    return (u"%d\N{DEGREE SIGN}" + fmt + "' %s") % (abs(deg), abs(min), dir)


class gshhs(object):
    """
    Read binary coastline in various resolutions from gshhs database.
   
    c = gshhs([12, 15.5, 53.25, 55], resolution='f',
              area_thresh=0., max_level=1)
   
    INPUTS:
        window = [lon_min, lon_max, lat_min, lat_max]
                 The window of the coastline extracted, same format as
                 a matplotlib axis.
       
        resolution = ['i' | 'c', 'l', 'h', 'f']
                     Resolution of coastline:
                     'c'ourse, 'l'ow, 'i'ntermidiate, 'h'igh and 'f'ine
       
        clip = [True | False]
               Flag to clip the coastline to the specified window.  This will
               cause the coastline segments to be shorter, as points outside the
               view are not include.
       
        max_level = [inf | integer]
                    The maximum gshhs level to include.  Level 1 is land only,
                    level 2 includes lakes, level 3 islands in lakes, etc.
       
        area_thresh = 0.
                      The minimum area a segment must have to be included.  Islands
                      or lakes smaller than area_thresh are excluded from the final
                      coastline polynomials.
   
    ATTRIBUTES:
       
        c.lon, c.lat are list of arrays containing the polynomials defining the land.
        c.level is the level of the polygon (1=land, 2=lake, etc.)
        c.greenwich is True if the polygon crosses greenwich, False otherwise
   
   
    METHODS:
       
        c.plot(ax=None, water_color=(0.7, 0.8, 1.0),
                 land_color=3*(0.5, ), edge_color = 3*(0.2,))
       
        plot polygons in a new figure (or in axis=ax if ax is not None).  Other
        inputs determine the colors of land, water, and the coastlines.
    """

    def _get_header(self):
        nb_i = struct.calcsize('i') # 4 bytes in 'int'
        nb_h = struct.calcsize('h') # 2 bytes in 'short int'
        if   self.version=='1.2':   hdr_len = nb_i*8 + nb_h*2
        elif self.version=='2.1.1': hdr_len = nb_i*11
        bhdr = self.f.read(hdr_len)
        if not bhdr:
            return None
        hdr = {}; n=0
        if self.version=='1.2':
          # integer header information (big-endien)
          hdr['id']     = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['n']      = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['level']  = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['west']   = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['east']   = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['south']  = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['north']  = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['area']   = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          # short integer header information
          hdr['greenwich'] = bool(struct.unpack('>h',bhdr[n:n+nb_h])[0]); n += nb_h
          hdr['source']   = struct.unpack('>h',bhdr[n:n+nb_h])[0];
        elif self.version=='2.1.1':
          # id:   Unique polygon id number, starting at 0
          # n:    Number of points in this polygon
          # flag: level + version << 8 + greenwich << 16 + source << 24 + river << 25
          # west, east, south, north: min/max extent in micro-degrees
          # area:      Area of polygon in 1/10 km^2
          # area_full: Area of original full-resolution polygon in 1/10 km^2
          # container: Id of container polygon that encloses this polygon (-1 if none)
          # ancestor:  Id of ancestor polygon in the full resolution set that was the source of this polygon (-1 if none)

          hdr['id']        = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['n']         = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['flag']      = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['west']      = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['east']      = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['south']     = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['north']     = struct.unpack('>i',bhdr[n:n+nb_i])[0]*1.0e-6; n += nb_i
          hdr['area']      = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['area_full'] = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['container'] = struct.unpack('>i',bhdr[n:n+nb_i])[0];        n += nb_i
          hdr['ancestor']  = struct.unpack('>i',bhdr[n:n+nb_i])[0]

          flag=hdr['flag']
          hdr['level']     = flag & 255        # Values: 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
          hdr['version']   = (flag >> 8) & 255 # Values: Should be 7 for GSHHS release 7
          hdr['greenwich'] = (flag >> 16) & 1  # Values: Greenwich is 1 if Greenwich is crossed
          hdr['source']    = (flag >> 24) & 1  # Values: 0 = CIA WDBII, 1 = WVS
          hdr['river']     = (flag >> 25) & 1  # Values: 0 = not set, 1 = river-lake and level = 2

        return hdr

    def _get_poly(self, hdr):
        n = hdr['n']
        n_i4 = struct.calcsize('i')
        bdata = self.f.read(2*n_i4*n)
        # unpack big-endien micro degrees
        data = array(struct.unpack('>'+str(2*n)+'i', bdata)[:], dtype='d')
        lon = array(+data[::2]*1.0e-6)
        lat = array(+data[1::2]*1.0e-6)
       
        if hdr['greenwich']:
            lon = (180.0/pi)*unwrap(lon*pi/180.0)
            if lon.max() > 360.0: lon -= 360.0
       
        if self.clip:
            # Clip data to window
            lon = clip(lon, self._ll_lon, self._ur_lon)
            lat = clip(lat, self.ll_lat, self.ur_lat)
           
            # shrink corners to a single point
            b = array([0.5, 0.5])
            idx = ( ((lon==self._ll_lon)&(lat==self.ll_lat))  # South-west corner
                  | ((lon==self._ur_lon)&(lat==self.ll_lat))  # South-east corner
                  | ((lon==self._ll_lon)&(lat==self.ur_lat))  # North-west corner
                  | ((lon==self._ur_lon)&(lat==self.ur_lat)) )# North-east corner
            idx=convolve(idx, b, mode='same')
            lon=lon[idx<1]; lat=lat[idx<1]
           
            # Only keep first and last points in series along each boundary
            b = array([0.25, 0.5, 0.25])
            if lon.size>=b.size:
              idx = (lon==self._ll_lon)
              idx = convolve(idx, b, mode='same')
              lon=lon[idx<1]; lat=lat[idx<1]
           
            if lon.size>=b.size:
              idx = (lon==self._ur_lon)
              idx = convolve(idx, b, mode='same')
              lon=lon[idx<1]; lat=lat[idx<1]
           
            if lon.size>=b.size:
              idx = (lat==self.ll_lat)
              idx = convolve(idx, b, mode='same')
              lon=lon[idx<1]; lat=lat[idx<1]
           
            if lon.size>=b.size:
              idx = (lat==self.ur_lat)
              idx = convolve(idx, b, mode='same')
              lon=lon[idx<1]; lat=lat[idx<1]
       
        return (lon,lat)
           
    def __init__(self, window, resolution='i',
                 clip=True, max_level=inf, area_thresh=0.,version='2.1.1'):

        self.version=version

        self.window = window
        self.ll_lon = self.window[0]
        self.ur_lon = self.window[1]
        self.ll_lat = self.window[2]
        self.ur_lat = self.window[3]
       
        self.clip = clip
        self.area_thresh = area_thresh
        self.max_level = max_level
       
        # if eastern limit is negative, wrap to (0, 360) range
        # western limit may stay negative if eastern limit is positive.
        if self.ur_lon < 0.0:
            # if ur_lon is negative, return lat/lon in negative deg. east.
            corr = -360.0
            self._ur_lon = self.ur_lon+360.0
            self._ll_lon = self.ll_lon+360.0
        else:
            corr = 0.0
            self._ur_lon = self.ur_lon
            self._ll_lon = self.ll_lon

        filename = 'gshhs_%s.b'%resolution
        gshhs_file=os.sep.join((os.environ['GSHHS_DATA'],filename))

        self.aspect = 1.0/float(cos(0.5*(self.ll_lat+self.ur_lat)*pi/180.))

        self.f = open(gshhs_file, 'rb')
        self.lon=[]; self.lat=[];
        self.area=[]; self.level=[]; self.greenwich=[]
        hdr = gshhs._get_header(self)
        while hdr:
            if (hdr['east']>self._ll_lon and hdr['west']<self._ur_lon and
                hdr['north']>self.ll_lat and hdr['south']<self.ur_lat and
                hdr['area']>self.area_thresh and hdr['level']<=self.max_level):
                self.area.append(hdr['area'])
                self.level.append(hdr['level'])
                self.greenwich.append(hdr['greenwich'])
                lon,lat = self._get_poly(hdr)
                lon += corr
                self.lon.append(lon)
                self.lat.append(lat)
            else:
                self.f.seek(2*4*hdr['n'],1)
            hdr = gshhs._get_header(self)
        self.f.close()


    def plot(self, ax=None, water_color=(0.7, 0.8, 1.0),
             land_color=3*(0.5, ), edge_color = 3*(0.2,)):
        pl.ioff()
        if not ax:
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_axis_bgcolor(water_color)
        for lon, lat, level in zip(self.lon, self.lat, self.level):
            if (level == 1):
                ax.fill(lon, lat, ec=edge_color, fc=land_color)
            else:
                ax.fill(lon, lat, ec=edge_color, fc=water_color)
        pl.axis(self.window)
        ax.set_aspect(self.aspect, adjustable='box', anchor='C')
        def xformat(x, pos=None): return deg2str(x, 'E', fmt="%2.0f")
        xformatter = FuncFormatter(xformat)
        ax.xaxis.set_major_formatter(xformatter)
        def yformat(y, pos=None): return deg2str(y, 'N', fmt="%2.0f")
        yformatter = FuncFormatter(yformat)
        ax.yaxis.set_major_formatter(yformatter)
        ax.fmt_xdata = lambda x: deg2str(x, 'E', fmt="%5.3f")
        ax.fmt_ydata = lambda y: deg2str(y, 'N', fmt="%5.3f")
   
    def get_lcc_proj(self, **kwargs):
        """return Basemap lcc projection based on window"""
        pass
   
    def get_merc_proj(self, **kwargs):
        """return Basemap mercator projection based on window"""
        pass


def add_coastline(ax,res='i',version='2.1.1',clip=True):
  '''Add land and lakes to axes'''

  g=gshhs(ax.axis(), resolution=res,area_thresh=0., max_level=inf,
          clip=clip,version=version)
  edge_color  = '#8f0d0d'
  land_color  = '#dcc683'
  water_color = '#90c7d5'
  for lon, lat, level in zip(g.lon, g.lat, g.level):
    if (level == 1):
      ax.fill(lon, lat, ec=edge_color, fc=land_color)
    else:
      ax.fill(lon, lat, ec=edge_color, fc=water_color)


def get_coastline0(xlim,ylim,res='h',version='2.1.1',clip=True):
  '''Load coastline from gshhs data folder, set in environment variable
  GSHHS_DATA. Supported data versions: 1.2 and 2.1.1

  There may be some problems with the data extraction... If so, use
  get_coastline (uses Basemap tools to extract data)
  '''

  lims=xlim[0],xlim[-1],ylim[0],ylim[-1]
  g=gshhs(lims, resolution=res,area_thresh=0., max_level=inf,clip=clip,
          version=version)
  res=np.zeros((0,2),'f')
  for lon, lat, level in zip(g.lon, g.lat, g.level):
    for i in range(len(lon)): res=np.vstack((res,(lon[i],lat[i])))
    res=np.vstack((res,(999.,999.)))

  return np.ma.masked_where(res==999,res).T


def get_coastline(xlim,ylim,res='h'):
  '''Load coastline from Basemap GSHHS data'''

  from mpl_toolkits.basemap import Basemap
  m=Basemap(projection='cyl',llcrnrlat=ylim[0],urcrnrlat=ylim[1],
            llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution=res)

  out=np.zeros((0,2),'f')
  for seg in m.coastsegs:
    for i in range(len(seg)): out=np.vstack((out,(seg[i][0],seg[i][1])))
    out=np.vstack((out,(999.,999.)))

  return np.ma.masked_where(out==999,out).T


def gen_bna0(fname,xylim,res='i'):
  '''Create Boundary File (BNA).
  This file is ready to be used by the NOAA oil spill model GNOME.

  This tool can deal with outdated gshhs datasets. Use gen_bna instead
  (based on Basemap).
  '''

  g=gshhs(xylim, resolution=res,area_thresh=0., max_level=2,clip=True)

  s=''
  n=0
  for lon, lat, level in zip(g.lon, g.lat, g.level):
    n+=1
    if (level == 1): # land
      sign=1
    else: # water
      sign=-1

    s+='"%s","1",%s\n' % (n,sign*len(lon))
    for i in range(len(lon)):
      s+='%8.4f, %8.4f\n' % (lon[i],lat[i])


  open(fname,'w').write(s)


def gen_bna(fname,xylim,res='i'):
  '''Create Boundary File (BNA).
  This file is ready to be used by the NOAA oil spill model GNOME.
  '''

  xlim,ylim=xylim[:2],xylim[2:]
  from mpl_toolkits.basemap import Basemap
  m=Basemap(projection='cyl',llcrnrlat=ylim[0],urcrnrlat=ylim[1],
            llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution=res)

  s=''
  n=0
  for i in range(len(m.coastpolygons)):
    n+=1
    level=m.coastpolygontypes[i]
    lon=m.coastpolygons[i][0]
    lat=m.coastpolygons[i][1]
    if (level == 1): # land
      sign=1
    else: # water
      sign=-1

    s+='"%s","1",%s\n' % (n,sign*len(lon))
    for i in range(len(lon)):
      s+='%8.4f, %8.4f\n' % (lon[i],lat[i])


  open(fname,'w').write(s)



if __name__ == '__main__':
    from pylab import *
    import time
    windows = ([12, 15.5, 53.25, 55],
              [-10, -5, 50, 60],
              # [-10, 10, 40, 60],
              [-95., -87., 28., 31.0])
    for window in windows:
        print 'window = ', window
        t1 = time.clock()
        m = gshhs(window=window, resolution='f', area_thresh=0., max_level=1)
        t2 = time.clock()
        print 'gshhs time: ' + str(t2-t1) + ' seconds'
        m.plot()
        t3 = time.clock()
        print 'plot time: ' + str(t3-t2) + ' seconds'
        t4 = time.clock()
        print 'show time: ' + str(t4-t3) + ' seconds'
        draw()
        print 'Total time: ' + str(t4-t1) + ' seconds'
    show()

