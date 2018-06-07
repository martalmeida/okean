'''
SMOS sea surface salinity data access from the Centre Aval de
Traitement des Donnees SMOS (CATDS), French ground segment for
the SMOS Level 3 and 4 data.

Login and password required. Check the webpage for more info:
http://www.catds.fr/Products/Products-access (support@catds.fr)

ESA's Soil Moisture Ocean Salinity (SMOS) Earth Explorer mission is
a radio telescope in orbit, but pointing back to Earth not space.
It's Microwave Imaging Radiometer using Aperture Synthesis (MIRAS)
radiometer picks up faint microwave emissions from Earth's surface to
map levels of land soil moisture and ocean salinity.

https://earth.esa.int/web/guest/missions/esa-operational-eo-missions/smos


M Marta-Almeida, feb 2014
'''


import ftplib
import os
import glob
from okean import netcdf, calc
import numpy as np
import datetime
try:
  from collections import OrderedDict
except:
  from okean .cookbook import odict as OrderedDict

import netrc

try:
  login,account,passw=netrc.netrc().authenticators('smos_ifremer')
except:
  msg='''login and pass missing in file $HOME/.netrc
Add to that file the lines:
machine smos_ifremer
login LOGIN
password PASSWORD\n
To get access to the data, check the site:
http://www.catds.fr/Products/Products-acces'''
  print(msg)

url='eftp.ifremer.fr'
path2='salinity/sss_smos_l#LEVEL#_V02/#YEAR#/10_Day_composite/Quarter_degree'
path1='salinity/sss_smos_l#LEVEL#_v01/#YEAR#/10_Day_Composite/Quarter_degree'
#                                -^-^-             -^-                        # differences in paths

class SMOSdownload:
  def __init__(self,dest=''):
    self.dest=dest
    self.f=False

  def connect(self):
    self.f=ftplib.FTP(url)
    self.f.login(login,passw)

  def goto(self,year,level=3,version=2):
    if not self.f: self.connect()

    if version==1: path=path1
    elif version==2: path=path2

    p=path.replace('#LEVEL#','%d'%level).replace('#YEAR#','%d'%year)
    if version==1: p=p.replace('composite','Composite') #  sick !!

    print('entering folder %s'%p)
    res=''
    try:
      self.f.cwd(p)
    except: res='cannot access %s'%p
    return res

  def list(self,year,level=3,version=2):
    msg=self.goto(year,level,version)
    if msg:
      print(msh)
      return

    res=self.f.nlst()
    for r in res: print(r)

  def destination_folder(self,year,level=3,version=2,gen=True):
    dest=os.path.join(self.dest,'%d'%year,'level_%d_version_%d'%(level,version))
    if gen and not os.path.isdir(dest):
        print('creating folder %s'%dest)
        os.makedirs(dest)

    return dest

  def download(self,year,level=3,version=2):
    msg=self.goto(year,level,version)
    if msg:
      print(msg)
      return

    files=self.f.nlst('*.nc')
    for r in files:
      destname=os.path.join(self.destination_folder(year,level,version),r)

      if os.path.isfile(destname):
        print('file %s already exists'%destname)
      else:
        print('downloading %s'%r)
        dest=open(destname,'wb')
        self.f.retrbinary('RETR '+r,dest.write)


class SMOS:
  def __init__(self,baseFolder='.'):
    self.path=baseFolder

  def extract(self,year,level=3,version=2,lims=False):
    a=SMOSdownload(self.path)
    p=a.destination_folder(year,level,version,gen=False)
    files=glob.glob(os.path.join(p,'*.nc'))
    files.sort()

    res=OrderedDict()
    for f in files:
      print(' -- extracting from %s'%f)
      if f==files[0]:
        lon=netcdf.use(f,'longitude')
        lat=netcdf.use(f,'latitude')
        if lims:
          i0,i1,j0,j1=calc.ij_limits(lon,lat,lims[:2],lims[2:],margin=2)
          ii='%d:%d'%(i0,i1+1)
          jj='%d:%d'%(j0,j1+1)
          lon=netcdf.use(f,'longitude',longitude=ii)
          lat=netcdf.use(f,'latitude',latitude=jj)
        else: ii,jj=':',':'

      date0=netcdf.nctime(f,'date_start')[0]
      date1=netcdf.nctime(f,'date_stop')[0]
      date=netcdf.nctime(f,'time')[0]
      sss=netcdf.use(f,'sss',longitude=ii,latitude=jj)
      res[date]=date0,date1,sss

    return lon,lat,res

if __name__=='__main__':
  import sys
  a=SMOSdownload()

  try:
    year=int(sys.argv[1])
    level=int(sys.argv[2])
    version=int(sys.argv[3])
  except:
    print('Usage: python smos.py YEAR LEVEL VERSION')
    print('Example: python smos.py 2012 3 2')
    sys.exit()

  a.download(year,level,version)
