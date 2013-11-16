from glob import glob
import os
from string import join as sjoin
import numpy as np
from okean import calc

def parse_locations(files):
  '''Reads locations file(s)
     locations file first line must be:
     #<label>, ex: #CITIES, #RIVERS
     Then each line should include <name>, <country>, <lon> <lat>
   '''

  if isinstance(files,basestring): files=[files]  

  locs=[]
  for c in files:
    f=open(c).readlines()
    for i in f:
      if i.strip().startswith('#'): continue
      tmp0=i.strip().split(',')
      name=tmp0[0].strip()

      tmp=tmp0[1].strip().split()
      country=tmp[1]

      for j in range(len(tmp)):
        if tmp[j].isdigit(): break

      country=sjoin(tmp[:j])

      lat=float(tmp[j])+float(tmp[j+1])/60.
      slat=tmp[j+2]
      if slat=='S': lat=-lat
      lon=float(tmp[j+3])+float(tmp[j+4])/60.
      slon=tmp[j+5]
      if slon=='W': lon=-lon

      locs+=[{'name': name,'country':country,'lon':lon,'lat':lat}]

  return locs 


class Locations:
  '''Loads locations from txt files
     Files must start with #<label> (ex #CITIES). Then each line
     should include <name>, <country>, <lon> <lat>, ex:
     Couto de Esteves, Portugal 40 45.46 N  8 18.44 W
  '''

  def __init__(self,what='',filespath='auto'):
    self.filespath=filespath
    self.what=what

    if self.filespath=='auto':
      here=os.path.dirname(os.path.abspath(__file__))
      self.filespath=os.path.join(here,'data')

    self.load()

  def load(self):
    files=glob(os.path.join(self.filespath,'*.txt'))

    # read 1st line to see if it matches self.what
    if self.what:
      Files=[]
      for f in files:
        l=open(f).readline()
        if l.strip()[1:].lower().startswith(self.what.lower()): Files+=[f]

      files=Files

    # parse files:
    self.locs=parse_locations(files)

    lon=[]
    lat=[]
    name=[]
    country=[]
    for c in self.locs:
      lon+=[c['lon']]
      lat+=[c['lat']]
      name+=[c['name']]
      country+=[c['country']]

    self.lon     = np.asarray(lon)
    self.lat     = np.asarray(lat)
    self.name    = np.asarray(name)
    self.country = np.asarray(country)

  def inside(self,x,y):
    '''Returns data inside the polygon x,y
      If len(x)==2, the polygon is the rectangle
      with limits x[0],x[1] and y[0],y[1]
      Return lon, lat, name and country

      -180<=x<=180; -90<=y<=90
    ''' 
    if len(x)==2: # rectangle
      x=x[0],x[1],x[1],x[0]
      y=y[0],y[0],y[1],y[1]

    inp=calc.inpolygon(self.lon,self.lat,x,y)
    return self.lon[inp],self.lat[inp],self.name[inp],self.country[inp]
