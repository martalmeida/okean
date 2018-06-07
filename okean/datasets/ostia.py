'''
OSTIA SST from myocean

Login and password required. Check the webpage for more info:
http://www.myocean.eu

The Operational Sea Surface Temperature and Sea Ice Analysis (OSTIA) system is run by the UK Met Office. Both a high resolution (1/20 - approx. 6 km) daily analysis of sea surface temperature (SST) and a reduced resolution (1/4 - approx. 28 km) daily analysis of sea surface temperature (SST) are produced for the global ocean and some lakes. SST anomalies are calculated on the reduced resolution analysis from the daily Pathfinder climatology. 

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
  login,account,passw=netrc.netrc().authenticators('ostia_myocean')
except:
  msg='''login and pass missing in file $HOME/.netrc
Add to that file the lines:
machine ostia_myocean
login LOGIN
password PASSWORD\n
To get access to the data, check the site:
http://www.myocean.eu'''
  print(msg)

url='data.ncof.co.uk'
path_nrt='SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001/#YEAR#/sst' # from 2007
path_rep='SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/#YEAR#/sst'

class OSTIAdownload:
  def __init__(self,dest=''):
    self.dest=dest
    self.f=False

  def connect(self):
    self.f=ftplib.FTP(url)
    self.f.login(login,passw)

  def goto(self,year,type):#level=3,version=2):
    '''type: rep or nrt'''
    if not self.f: self.connect()

    if type=='rep': path=path_rep
    elif type=='nrt': path=path_nrt

    p=path.replace('#YEAR#','%d'%year)
############    if version==1: p=p.replace('composite','Composite') #  sick !!

    print('entering folder %s'%p)
    res=''
    try:
      self.f.cwd(p)
    except: res='cannot access %s'%p
    return res

  def list(self,year,type):#######level=3,version=2):
    msg=self.goto(year,type)###########level,version)
    if msg:
      print(msh)
      return

    res=self.f.nlst()
    for r in res: print(r)

  def destination_folder(self,year,type,gen=True):##level=3,version=2,gen=True):
    dest=os.path.join(self.dest,'%d'%year,'%s'%type)#####level_%d_version_%d'%(level,version))
    if gen and not os.path.isdir(dest):
        print('creating folder %s'%dest)
        os.makedirs(dest)

    return dest

  def download(self,year,type):##########level=3,version=2):
    msg=self.goto(year,type)#########level,version)
    if msg:
      print(msg)
      return

    files=self.f.nlst('*.nc')
    if len(files)==0: files=self.f.nlst('*.bz2')
    for r in files:
      destname=os.path.join(self.destination_folder(year,type),r)

      if os.path.isfile(destname):
        print('file %s already exists'%destname)
      else:
        print('downloading %s'%r)
        dest=open(destname,'wb')
        self.f.retrbinary('RETR '+r,dest.write)


class OSTIA:
  def __init__(self,baseFolder='.'):
    self.path=baseFolder

  def extract(self,year,type,lims=False,date=False):
    a=OSTIAdownload(self.path)
    p=a.destination_folder(year,type)#########level,version,gen=False)
    files=glob.glob(os.path.join(p,'*.nc'))
    files.sort()

    res=OrderedDict()
    c=-1
    for f in files:
      c+=1
      print(' -- extracting from %s'%f)
      if c==0:
        lon=netcdf.use(f,'lon')
        lat=netcdf.use(f,'lat')
        if lims:
          i0,i1,j0,j1=calc.ij_limits(lon,lat,lims[:2],lims[2:],margin=2)
          ii='%d:%d'%(i0,i1+1)
          jj='%d:%d'%(j0,j1+1)
          lon=netcdf.use(f,'lon',lon=ii)
          lat=netcdf.use(f,'lat',lat=jj)
        else: ii,jj=':',':'


      date=netcdf.nctime(f,'time')[0]
      u=netcdf.use(f,'analysed_sst',lon=ii,lat=jj)
      if c==0:
        sst=np.ma.zeros((len(files),)+u.shape,u.dtype)
        time=np.zeros(len(files),datetime.datetime)

      sst[c]=u
      time[c]=date

###      date0=netcdf.nctime(f,'date_start')[0]
###      date1=netcdf.nctime(f,'date_stop')[0]
#      res[date]=sst
#####  return lon,lat,res

    return time,lon,lat,sst


if __name__=='__main__':
  import sys
  a=OSTIAdownload()

  try:
    year=int(sys.argv[1])
    type=sys.argv[2]
  except:
    print('Usage: python ostia.py YEAR TYPE')
    print('TYPE: nrt or rep')
    print('Example: python ostia.py 2013 nrt')
    sys.exit()

  a.download(year,type)#########level,version)
