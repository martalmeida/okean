'''
Real-time, global, sea surface temperature (RTG_SST_HR) analysis
http://polar.ncep.noaa.gov/sst/rtg_high_res/

A daily, high-resolution, real-time, global, sea surface temperature
(RTG_SST) analysis has been developed at the National Centers for
Environmental Prediction/Marine Modeling and Analysis Branch
(NCEP / MMAB). The analysis was implemented in the NCEP parallel
production suite 16 August 2005. It became fully operational on
September 27, 2005.

M Marta-Almeida, feb 2014
'''


import ftplib
import os
import glob
from okean import netcdf, calc,gribu
import numpy as np
import datetime
try:
  from collections import OrderedDict
except:
  from okean .cookbook import odict as OrderedDict

#import netrc
#
#try:
#  login,account,passw=netrc.netrc().authenticators('ostia_myocean')
#except:
#  msg='''login and pass missing in file $HOME/.netrc
#Add to that file the lines:
#machine ostia_myocean
#login LOGIN
#password PASSWORD\n
#To get access to the data, check the site:
#http://www.myocean.eu'''
#  print msg
#
#url='data.ncof.co.uk'
#path_nrt='SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001/#YEAR#/sst' # from 2007
#path_rep='SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/#YEAR#/sst'
#
#ftp://polar.ncep.noaa.gov/pub/history/sst/ophi/rtg_sst_grb_hr_0.083.200510.gz
###url='ftp://polar.ncep.noaa.gov'
url='polar.ncep.noaa.gov'
path='pub/history/sst/ophi'
path_mask='pub/history/sst/ophi/lsmask'

class RTG_SSTdownload:
  def __init__(self,dest=''):
    self.dest=dest
    self.f=False

  def connect(self):
    self.f=ftplib.FTP(url)
    self.f.login()##login,passw)

  def goto(self,type=False):#####,year,type):#level=3,version=2):
    if not self.f: self.connect()

    if type=='mask': p=path_mask
    else: p=path
##    if type=='rep': path=path_rep
##    elif type=='nrt': path=path_nrt
##
##    p=path.replace('#YEAR#','%d'%year)
##############    if version==1: p=p.replace('composite','Composite') #  sick !!
##
    print('entering folder %s'%p)
    res=''
    try:
      self.f.cwd(p)
    except: res='cannot access %s'%p
    return res

  def list(self,year):####,type):#######level=3,version=2):
    msg=self.goto()####year,type)###########level,version)
    if msg:
      print(msh)
      return

    res=self.f.nlst('rtg_sst_*.%d*'%year)
    for r in res: print(r)

  def destination_folder(self,year,gen=True):#,type,gen=True):##level=3,version=2,gen=True):
    if year=='mask':
      dest=os.path.join(self.dest,'lsmask')
    else:
      dest=os.path.join(self.dest,'%d'%year)

    if gen and not os.path.isdir(dest):
        print('creating folder %s'%dest)
        os.makedirs(dest)

    return dest

  def download(self,year):####,type):##########level=3,version=2):
    msg=self.goto()####year,type)#########level,version)
    if msg:
      print(msg)
      return

    files=self.f.nlst('rtg_sst_*.%d*'%year)
####    if len(files)==0: files=self.f.nlst('*.bz2')
    for r in files:
      destname=os.path.join(self.destination_folder(year),r)

      if os.path.isfile(destname):
        print('file %s already exists'%destname)
      else:
        print('downloading %s'%r)
        dest=open(destname,'wb')
        self.f.retrbinary('RETR '+r,dest.write)

  def download_mask(self):
    msg=self.goto('mask')
    if msg:
      print(msg)
      return    

    files=['ls.dat']
    for r in files:
      destname=os.path.join(self.destination_folder('mask'),r)

      if os.path.isfile(destname):
        print('file %s already exists'%destname)
      else:
        print('downloading %s'%r)
        dest=open(destname,'wb')
        self.f.retrbinary('RETR '+r,dest.write)


class RTG_SST:
  def __init__(self,baseFolder='.'):
    self.path=baseFolder

  def load_mask(self,shape=False):
    fmask='mask_rtg.npy'
    if os.path.isfile(fmask):
      print('loading mask from %s'%fmask)
      return np.load(fmask)


    f=os.path.join(self.path,'lsmask','ls.dat')
    print('reading mask file %s'%f)
    w=open(f).read()
    w=w.replace('\n','')
    mask=np.asarray([int(i) for i in w])
    if not shape is False:
      mask.shape=shape

    mask=np.flipud(mask)

    mask.dump(fmask)  
    print('  done')
    return mask  



  def extract(self,year,lims=False,date=False):
    a=RTG_SSTdownload(self.path)
    p=a.destination_folder(year)
    if date:
      files=glob.glob(os.path.join(p,'*.%s'%date.strftime('%Y%m%d')))
    else:
      files=glob.glob(os.path.join(p,'*'))

    files.sort()

    # find ntimes:
    nt=0
    for f in files:
      nt=nt+len(gribu.findvar(f,'temp'))

    c=-1
    for f in files:
      print(' -- extracting from %s'%f)
      q=gribu.findvar(f,'temp')
      for V in q: # times per file
        c+=1

        if c==0:
          lat,lon=V.latlons()
          mask=self.load_mask(lon.shape)
          if lims:
            if lims[1]<0 and lims[1]<0:
              lon=lon-360 # if not both lon lims <0,, a few more lines are needed !
            print('calc ij inds...')
            ijname='ijinds.npy'
            if os.path.isfile(ijname):
              i0,i1,j0,j1=np.load(ijname)
            else:
              i0,i1,j0,j1=calc.ij_limits(lon,lat,lims[:2],lims[2:],margin=2)
              np.asarray([i0,i1,j0,j1]).dump(ijname)

            print('done')
            lon=lon[j0:j1,i0:i1]
            lat=lat[j0:j1,i0:i1]
            mask=mask[j0:j1,i0:i1]

          else: i0=False

          sst=np.ma.zeros((nt,)+lon.shape,lon.dtype)
          time=np.zeros(nt,datetime.datetime)

        if not i0 is False:
          s=V.values[j0:j1,i0:i1]
        else: s=V.values

        s=np.ma.masked_where(mask==3,s)

        sst[c]=s
        time[c]=V.analDate
        print('=> done for %s'%time[c].isoformat(' '))

    return time,lon,lat,sst


if __name__=='__main__':
  import sys
  a=RTG_SSTdownload()

  try:
    year=sys.argv[1]
    if year!='mask': year=int(year)
  except:
    print('Usage: python rtg_sst.py YEAR')
    print('or, to get land sea mask: python rtg_sst.py mask')
    print('Example: python rtg_sst.py 2013')
    sys.exit()

  if year=='mask':
    a.download_mask()
  else:
    a.download(year)
