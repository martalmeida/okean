'''
Tools to download GFS and store data from NOMADS server and to extract
the fields required for ocean model atmospheric bulk forcing.
Can retrieve files and data for/from analysis or forecast simulations.

Created by Martinho MA at UFBA, Salvador, Brasil, 2009
Improved in Feb 2011, TAMU, Texas, USA
Updated Aug 2013 to new OKEAN features, Guayaquil, Ecuador
'''

import os
import datetime
import numpy as np
from okean import dateu, air_sea, gribu
try:
  from collections import OrderedDict as odict
except:
  from okean.cookbook import odict

CNVGRIB='cnvgrib'


class Data:
  '''
  Very simple class for x,y,v data.
  Includes the attributes x,y and data.
  By default, also includes units and info.
  '''

  def __init__(self,x,y,data,units,info=''): 
    self.x     = x
    self.y     = y
    self.data  = data
    self.units = units
    self.info  = info


class GFSDownload:
  '''
  Support for GFS data download

  Inputs:
    basefolder, musc contain:
      - conf file tags.info
      - scripts folder with get_inv.pl and get_grib.pl, needed for downloading
  '''

  def __init__(self,basefolder,**kargs):
    self.basefolder    = basefolder
    self.datafolder    = os.path.join(self.basefolder,'gribs')
    self.logFolder     = os.path.join(self.basefolder,'log')
    self.scriptsFolder = os.path.join(self.basefolder,'scripts')

    self.attmax        = 10
    self.ngetBefore    = 7

    # vars for ROMS-AGRIF
    self.egrep='| egrep "(PRATE|ALBDO|RH:2 m|TMP:2 m|TMP:s|PRES:s|DSWRF:s|DLWRF:s|UGRD:10 m |VGRD:10 m )" |'
    # other vars for ROMS:
    # Total Cloud Cover            TCDC:e
    # Upward Short-Wave Rad. Flux  USWRF:s
    # Upward Long-Wave Rad. Flux   ULWRF:s
    # Total Precipitation          APCP:s
    # Ground Heat Flux             GFLUX:s
    # Sensible Heat Net Flux       SHTFL:s
    # Latent Heat Net Flux         LHTFL:s
#    self.egrep='| egrep "(PRATE|ALBDO|RH:2 m|TMP:2 m above|TMP:s|PRES:s|DSWRF:s|DLWRF:s|UGRD:10 m a|VGRD:10 m a|TCDC:e|USWRF:s|ULWRF:s|APCP:s|GFLUX:s|SHTFL:s|LHTFL:s )" |'

    self.egrep='| egrep "(PRATE:surface.*ave|ALBDO|RH:2 m above ground|:TMP:2 m above|:TMP:surface|:PRES:surface|DSWRF:s|DLWRF:s|UGRD:10 m a|VGRD:10 m a|TCDC:e|USWRF:s|ULWRF:s|APCP:s|GFLUX:s|SHTFL:s|LHTFL:s)" |'

    # get info from tags info file:
    self.get_tags()

    # info available from get_tags:
    if 'url'         in kargs.keys(): self.url           = kargs['url']
    if 'nameTagDest' in kargs.keys(): self.nameTagDest   = kargs['nameTagDest']
    if 'nameTagSrc'  in kargs.keys(): self.nameTagSrc    = kargs['nameTagSrc']
    if 'gribVersion' in kargs.keys(): self.gribVersion   = kargs['gribVersion']
    if 'invExt'      in kargs.keys(): self.invExt        = kargs['invExt']
    if 'nforec'      in kargs.keys(): self.nforec        = kargs['nforec']
    if 'dt_start'    in kargs.keys(): self.dt_startc     = kargs['dt_start']
    if 'dt_sim'      in kargs.keys(): self.dt_sim        = kargs['dt_sim']
    # other arguments:
    if 'logpath'     in kargs.keys(): self.logFolder     = kargs['logpath']
    if 'scriptspath' in kargs.keys(): self.scriptsFolder = kargs['logpath']
    if 'attmax'      in kargs.keys(): self.attmax        = kargs['attmax']
    if 'egrep'       in kargs.keys(): self.egrep         = kargs['egrep']
    if 'nbefore'     in kargs.keys(): self.ngetBefore    = kargs['nbefore']


  def get_tags(self,quiet=1):
    '''
    parse GFS download configuration file ("tags.info")
    '''
    f=os.path.join(self.basefolder,'tags.info')

    tags  ='url', 'nameTagDest', 'nameTagSrc', 'gribVersion', 'invExt', 'dt_start', 'dt_sim', 'nforec'
    types = str , str          , str         , int          , str     , int       , int     , int

    if os.path.isfile(f):
      f=open(f).readlines()
      for i in range(len(f)-1):
        for j in range(len(tags)):
          if f[i].find(tags[j])>=0 :
            attname  = tags[j]
            attvalue = types[j](f[i+1].strip())
            setattr(self,attname,attvalue)
            if not quiet: print(attname,' : ',attvalue)


  def nameof(self,type,date,hour_start,hour_sim):
    '''
    Name of source and destination files
    '''
    if isinstance(date,datetime.datetime): date=date.strftime('%Y%m%d')

    if type in ('src','srcfname','srcpath'):
      ptmp=self.url
      ntmp=self.nameTagSrc
    elif type in ('dest','destfname','destpath'):
      ptmp=os.path.join(self.datafolder,'#DATE#')
      ntmp=self.nameTagDest

    ptmp=ptmp.replace('#DATE#',date)
    ptmp=ptmp.replace('#HOUR_START#','%02d' % hour_start)

    if date>='20150115': # actually, hour sim is not in the path, only in filname.... so, to this in ntmp!
      ptmp=ptmp.replace('#HOUR_SIM#','%03d' % hour_sim)
    else:
      ptmp=ptmp.replace('#HOUR_SIM#','%02d' % hour_sim)

    ptmp=ptmp.replace('#VERSION#','%d' % self.gribVersion)

    ntmp=ntmp.replace('#DATE#',date)
    ntmp=ntmp.replace('#HOUR_START#','%02d' % hour_start)
    if date>='20150115':
      ntmp=ntmp.replace('#HOUR_SIM#','%03d' % hour_sim)
    else:
      ntmp=ntmp.replace('#HOUR_SIM#','%02d' % hour_sim)
    ntmp=ntmp.replace('#VERSION#','%d' % self.gribVersion)

    if   type in ('srcfname','destfname'): name = ntmp
    elif type in ('srcpath','destpath'):   name = ptmp
    elif type in ('src','dest'): name = os.path.join(ptmp,ntmp)

    return {'name':name,'hour_start':hour_start,'hour_sim':hour_sim,'date':date,'type':type}


  def src2dest(self,src):
    '''
    Name of destination file from source file
    '''
    type       = src['type']
    if   'name' in type: type='destname'
    elif 'path' in type: type='destpath'
    else: type='dest'

    date       = src['date']
    name       = src['name']
    hour_start = src['hour_start']
    hour_sim   = src['hour_sim']

    return self.nameof(type,date,hour_start,hour_sim)


  def prev_option(self,src):
    '''
    Previous source file corresponding to the same time.
    Example, if the runs start each 6h and give outputs each 3h, the
    previous option of the 3 hours forecast of the run starting at 12h
    (time=15h) is the 9h forecast of the run starting at 6h.
    '''
    date       = src['date']
    hour_start = src['hour_start']
    hour_sim   = src['hour_sim']

    if hour_start==0:
      date       = dateu.next_date(date,-1)
      hour_start = 24-self.dt_start
      hour_sim   = hour_sim+self.dt_start
    else:
      hour_start = hour_start-self.dt_start
      hour_sim   = hour_sim+self.dt_start

    return self.nameof('src',date,hour_start,hour_sim)


  def daily_files(self,date,FA='af'):
    '''
    Files to download each day (date) for analysis (FA='a') and
    forecast (FA='f'), and corresponding destination files.
    '''
    aux=range(self.dt_sim,self.dt_start+self.dt_sim,self.dt_sim) # 3,6 usually
    H=odict()
    if FA=='af' or FA=='f':
      H[0]  = range(self.dt_sim,self.nforec*24+self.dt_sim,self.dt_sim)
    elif FA=='a':
      H[0]  = aux

    if FA=='a' or FA=='af':
      for i in range(self.dt_start,24,self.dt_start):
        H[i]  = aux

    srcFnames=[]
    destFnames=[]
    Paths=[]
    for h in H.keys(): # hour start
       for hh in H[h]: # hour sim
        sfname = self.nameof('src', date,h,hh)
        dfname = self.nameof('dest',date,h,hh)

        srcFnames+=[sfname]
        destFnames+=[dfname]

    return srcFnames,destFnames


  def __download_fast_once(self,source,checkinv,log,del1,quiet):
    '''
    Donloads the self.egrep variables from source file
    Used by download_fast
    '''
    err=True
    src=source['name']

    # download scripts:
    get_inv  = os.path.join(self.scriptsFolder,'get_inv.pl')
    get_grib = os.path.join(self.scriptsFolder,'get_grib.pl')

    # inventory filename:
    source_inv=src+'.'+self.invExt

    # detinations files:
    dest=self.src2dest(source)['name']
    dest1=dest[:-1]+'%d' % self.gribVersion
    dest2=dest[:-1]+'2'

    if not quiet: print('getting ',dest2)
    if not os.path.isfile(dest2):
      if not quiet: print('Downloading '+src)
      cmdGet  = ' '.join((get_inv,source_inv,self.egrep,get_grib,src,dest1))
      mayDownload=1
      if checkinv:
        mayDownload=0
        cmdCheck=' '.join((get_inv,source_inv))
        res=os.system(cmdCheck+'>/dev/null 2>&1')
        if res is 0: mayDownload=1

      if mayDownload:
        # create destination folder:
        p=os.path.dirname(dest2)
        if not os.path.isdir(p): os.makedirs(p)

        # downloading:
        open(log,'a').write(':: Downloading '+src+'\n')
        err=os.system(cmdGet  + '>> '+log+' 2>&1')

        # convert to grib 2 and remove grib1 if required:
        if self.gribVersion==1:
          # converting:
          cmdConv = ' '.join((CNVGRIB,'-g12',dest1,dest2))
          open(log,'a').write(':: Converting to '+dest2+'\n')
          os.system(cmdConv + '>> '+log+' 2>&1')

          # removing grib1:
          if del1:
            open(log,'a').write(':: Removing to '+dest1+'\n')
            try:
              os.remove(dest1)
            except:
              open(log,'a').write(':: Problems removing '+dest1+'\n')

      else:
        if not quiet: print('  -- cannot download')

    else:
      if not quiet: print('  -- already exists')
      err=False

    return err,dest2


  def download_fast(self,date,FA='af',del1=True,checkinv=False,quiet=True,prevopt=True):
    '''
    Downloads all the self.egrep variables for date. By default both
    analysis and forecast data is downloaded.
    Used by download_current and download_range (use then instead)

    Inputs:
      date, FA ('f'orecast of 'a'nalysis or 'af')
      del1, if the original version is grib1, the conversion to grib2
        is done and the version 1 files are removed if del1 is True
      checkinv will check if inv file exists... increases speed for
        current day downloads!!
      prevopt, check for previous options file is the "best one" is note
        present
    '''

    # daily files:
    targets,destinations=self.daily_files(date,FA)

    # start download log:
    log=os.path.join(self.logFolder,'download.log')
    if not os.path.isdir(self.logFolder): os.makedirs(self.logFolder)
    open(log,'a').write('::::'+dateu.currday(local=True).strftime("%b %d %Y %H:%M:%S")+'\n')


    for i in range(len(targets)):
      err,dest2best=self.__download_fast_once(targets[i],checkinv,log,del1,quiet)

      # if download failes, check previous files, nAttenpts, if best file is not present:
      nAttempts=0
      Target=targets[i]
      while err and nAttempts<=self.attmax and prevopt:
        nAttempts+=1
        prevtarget=self.prev_option(Target)
        Target=prevtarget
        err,dest2=self.__download_fast_once(prevtarget,False,log,del1,quiet)
        if not err:
          # link file:
          open(log,'a').write(':: Linking '+dest2+'\n')
          os.symlink(os.path.realpath(dest2),os.path.realpath(dest2best))
          if not quiet: print('linking ',dest2,dest2best)


  def download_current(self,date=False,del1=True, quiet=True):
    '''
    Download files for today analysis and forecast.
    If date is not provided (the current day is used), also download
    analysis data of last self.ngetBefore days.

    Inputs:
      date, default is the current day
      del1, if the original version is grib1, the conversion to grib2
        is done and the version 1 files are removed if del1 is True
      quiet, print info flag
    '''

    if not date: date=dateu.currday()
    else: self.ngetBefore=1 # just get the selected day!

    for i in range(self.ngetBefore):
      day=dateu.next_date(date,-i,samefmt=False)
      if not quiet: print('Downloading GFS files for date '+day.strftime('%Y%m%d'))
      if i==0:
        # download files for analysis and forecast:
        self.download_fast(day,FA='a',del1=del1,checkinv=True,quiet=quiet,prevopt=False)
        self.download_fast(day,FA='f',del1=del1,checkinv=True,quiet=quiet,prevopt=False)
      else:
        # download any missing analysis file from prev days:
        self.download_fast(day,FA='a',del1=del1,checkinv=False,quiet=quiet,prevopt=True)


  def download_range(self,date0,date1,quiet=True):
    '''
    Download analysis data in the interval [date0... date1[
    '''
    dates=dateu.drange(date0,date1,inclast=False)
    for date in dates:
      self.download_fast(date,FA='a',del1=True,checkinv=True,quiet=quiet,prevopt=True)


  def FTP_download_progress(self,sofar,size,totalsize):
    '''used by FTP_daily_download: DEPRECATED'''
    p=sofar*size/float(totalsize) * 100
    if int(ceil(mod(p,10)))==10:
      if int(ceil(p))>self.tmp:
        print('%3d %% (%d secs)' % (int(ceil(p)), time.time()-self.tmp3))
        self.tmp=int(ceil(p))
        self.tmp3=time.time()

    if int(p)>self.tmp2:
      sys.stdout.write('.')
      sys.stdout.flush()

    self.tmp2=int(p)


  def FTP_daily_download(self,date,FA='a'):
    '''DEPRECATED'''

    import urllib2
    import urllib

    self.tmp=0
    self.tmp2=0
    self.tmp3=time.time()

    H  = 0,6,12,18
    H={}
    if FA=='f':
      H[0]  = range(0,self.nforec*24+self.dt,self.dt)
    else:
      H[0]  = 3,6

    H[6]  = 3,6
    H[12] = 3,6
    H[18] = 3,6

    path = os.path.join(self.url,'gfs'+date)
    try:
      url=urllib2.urlopen(path)
      url=' '.join(url.readlines())
      mayGo=True
    except:
      mayGo=False

    if mayGo:
      Fnames=[]
      for h in H.keys():
         for hh in H[h]:
          fname    = 'gfs.t%02dz.pgrbf%02d' % (h,hh)
          Fnames+=[fname]


      for F in Fnames:
        destination=os.path.join(self.datafolder,date,F)
        target=os.path.join(path,F)
        if not os.path.isdir(os.path.join(self.datafolder,date)): os.makedirs(os.path.join(self.datafolder,date))
        if F in url and not os.path.isfile(destination):
         print('Getting ',target,destination)
         urllib.urlretrieve(target,destination,self.FTP_download_progress)
         self.tmp=0
         self.tmp2=0
         self.tmp3=time.time()



def get_date(fname):
  '''
  Reads date from saved fname, which is in a format like
  etc/20110101/gfs_20110101_18_06.grib2
  or
  etc/gfs_20110101_18_06.grib2
  '''
  if isinstance(fname,basestring):
    # expected date to be in filename or in path of file
    # expected hour_start and hour sim to be in filename!
    # hour_start must have 2 digist and hour_sim musta have at least 2 digits!

    date=False
    hour_start=False
    hour_sim=False
    name=os.path.basename(fname)
    path=os.path.dirname(fname)

    # search  date in filename
    for i in range(len(name)-8+1):
      if name[i].isdigit():
        try:
          int(name[i:i+8])
          date=name[i:i+8]
          break
        except: pass

    # search for hour start:
    if not date is False: i0=i+8
    else: i0=0
    for i in range(i0,len(name)-2+1):
      if name[i].isdigit():
        try:
          hour_start=int(name[i:i+2])
          break
        except: pass


    # search for hour sim:
    if not hour_start is False:
      for i in range(i+2,len(name)-2+1):
        if name[i].isdigit():
          try:
            hour_sim=int(name[i:i+3])
            break
          except:
            hour_sim=int(name[i:i+2])
            break


    # search  date in path
    if date is False:
      for i in range(len(path)-8+1):
        if path[i].isdigit():
          try:
            int(path[i:i+8])
            date=path[i:i+8]
            break
          except: pass


  h=hour_start+hour_sim
  return dateu.parse_date(date)+datetime.timedelta(hours=h)



class GFSData:
  '''
  GFS data extraction
  '''
  def __init__(self,basefolder):
    self.basefolder=basefolder

  def __files(self,date0,date1=False,FA='a',nforec='auto'):
    '''
    Used by files_analysis and files_forecast
    '''

    if FA=='f': date1=False
    if nforec=='auto': args={}
    else: args={'nforec':nforec}

    a=GFSDownload(basefolder=self.basefolder,**args)
    if date1 is False: dates=[date0]
    else: dates=dateu.drange(date0,date1)

    files=[]
    time=[]
    isbest=[]

    # first file, 00h data (last of previous day)
    datePrev=dateu.next_date(date0,-1)
    file0=a.daily_files(datePrev,FA='a')[1][-1]['name']

    for d in dates:
      Src,Dest=a.daily_files(d,FA=FA)
      for dest in Dest:
        files+=[dest['name']]

    if files: files=[file0]+files

    for f in files:
      time+=[get_date(f)]
      if os.path.isfile(f): isbest+=[not os.path.islink(f)]
      else: isbest+=[None]

    return files,time,isbest


  def files_analysis(self,date0,date1=False):
    '''
    Returns all the required filenames for analysis data
    Also returns the datetime of each file and isbest flag.
    isbest is True if the file is not a symbolic link, which
    means is a file from a previous run (having however the correct time)

    if date1 is not False, all dates in range [date0 ... date1[ will be
    considered (date1==False means date1=date0+1day)

    dates can be strings (yyyymmdd) or datetime objects
    '''
    return self.__files(date0,date1,FA='a')


  def files_forecast(self,date,nforec='auto'):
    '''
    Returns all the required filenames for analysis data
    Also returns the datetime of each file and isbest flag,
    (see files_analysis for more info).

    if npred is auto, use all files (as set in config file tags.info)
    date as string or datetime

    date can be string (yyyymmdd) or datetime object
    '''
    return self.__files(date,FA='f',nforec=nforec)


  def __data(self,date0,date1=False,FA='a',nforec='auto',xlim=False,ylim=False,quiet=True):
    '''
    Used by data_analysis and data_forecast
    '''

    if FA=='a':   files,time,isbest=self.files_analysis(date0,date1)
    elif FA=='f': files,time,isbest=self.files_forecast(date0,nforec)

    res=odict()
    missing=odict()

    for i in range(len(files)):
      if not quiet: print('|-> getting from file %s' % files[i])
      if isbest[i]!=None:
        data=gfs_file_data(files[i],xlim,ylim,quiet=quiet)
        data['INFO_isbest'] = isbest[i]
        data['INFO_file']   = files[i]
        res[time[i]]=data
        if not quiet: print('     ** isbest =',isbest[i],'**')
      else:
        if not quiet: print('     file is missing !!')
        missing[time[i]]=files[i]

    return res,missing


  def data_analysis(self,date0,date1=False,xlim=False,ylim=False,quiet=True):
    '''
    Returns atm data for analysis forcing in the dates range [date0... date1[
    If date1 is False, only date0 is considered (one day).

    xlim and ylim are data x,y limits. If not given, all data in file
    is given.
    '''
    return self.__data(date0,date1,'a',None,xlim,ylim,quiet)


  def data_forecast(self,date,nforec='auto',xlim=False,ylim=False,quiet=True):
    '''
    Returns atm data for forecast forcing for day date.

    xlim and ylim are data x,y limits. If not given, all data in file
    is given.
    '''
    return self.__data(date,False,'f',nforec,xlim,ylim,quiet)


  def is_ready(self,date,FA,nforec='auto'):
    '''
    Check is required data files with data for analysis (FA='a') or
    forecast (FA='f'), for day date, are present.
    nforec is the number of forecast days, if auto, the default from
    config file (tags.info) is used.
    '''

    if FA=='a':
      files,time,isbest=self.files_analysis(date)
    elif FA=='f':
      files,time,isbest=self.files_forecast(date,nforec)

    return not any([i is None for i in isbest])


def is_ready(basefolder,date,FA,nforec='auto'):
  '''
  Check is required data files with data for analysis or forecast are
  all present

  Inputs
    basefolder, GFS data, scripts and config base folder
    date, datetime or string (yyyymmdd)
    FA, 'a' or 'f' for analysis and forecast
    nforec, number of forecast days, default auto from the config file
  '''

  return GFSData(basefolder).is_ready(date,FA,nforec)


def gfs_file_data(fname,xlim=False,ylim=False,quiet=False):
  '''
  Returns bulk data from one GFS file
  '''
  # the way to ectract differes if using pygrib or grib2 ! so check first:
  try:
    import pygrib
    isPygrib=True
  except: isPygrib=False


  out={}

  # T air 2m [K->C]
  if not quiet: print(' --> T air')
  if isPygrib:
    #x,y,tair=gribu.getvar(fname,'temperature',tags=(':2 metre',),lons=xlim,lats=ylim)
    #newest gribu:
    x,y,tair=gribu.getvar(fname,'2t',lons=xlim,lats=ylim)
  else:
    x,y,tair=gribu.getvar(fname,'temperature',tags=(':2 m','TMP'),lons=xlim,lats=ylim)
  tair=tair-273.15
  out['tair']=Data(x,y,tair,'C')

  # R humidity 2m [%-->0--1]
  if not quiet: print(' --> R humidity')
  if 0:
    # kg/kg
    x,y,rhum=gribu.getvar(fname,'humidity',tags=(':2 m','kg'),lons=xlim,lats=ylim)
    rhum=rhum/air_sea.qsat(tair)
    rhum=np.where(rhum>1.0,1.0,rhum)
  else:
    # %
    #x,y,rhum=gribu.getvar(fname,'humidity',tags=('2 m','%'),lons=xlim,lats=ylim)
    x,y,rhum=gribu.getvar(fname,'2r',lons=xlim,lats=ylim)
    rhum=rhum/100.

  out['rhum']=Data(x,y,rhum,'0--1')

  # surface pressure [Pa]
  if not quiet: print(' --> Surface pressure')
  #x,y,pres=gribu.getvar(fname,'pressure',tags='surface',lons=xlim,lats=ylim)
  x,y,pres=gribu.getvar(fname,'sp',lons=xlim,lats=ylim)
  out['pres']=Data(x,y,pres,'Pa')

  # P rate [kg m-2 s-1 -> cm/d]
  if not quiet: print(' --> P rate')
  x,y,prate=gribu.getvar(fname,'prate',tags='avg',lons=xlim,lats=ylim)
  # Conversion kg m^-2 s^-1  to cm/day
  prate=prate*86400*100/1000.
  prate=np.where(abs(prate)<1.e-4,0,prate)
  out['prate']=Data(x,y,prate,'cm/d')

  # Net shortwave flux  [W/m^2]
  if not quiet: print(' --> Net shortwave flux')
  if not quiet: print('       SW down')
  if isPygrib:
    #x,y,sw_down = gribu.getvar(fname,'',tags='Downward short-wave radiation flux',lons=xlim,lats=ylim)
    x,y,sw_down = gribu.getvar(fname,'dswrf',lons=xlim,lats=ylim)
  else:
    x,y,sw_down = gribu.getvar(fname,'downward short-wave',lons=xlim,lats=ylim)

  if not quiet: print('       SW up')
  x,y,sw_up   = gribu.getvar(fname,'uswrf',tags='surface',lons=xlim,lats=ylim)
  if sw_up is False:
    if not quiet: print('       SW up not found: using albedo')
    #x,y,albedo  = gribu.getvar(fname,'albedo',lons=xlim,lats=ylim)
    x,y,albedo  = gribu.getvar(fname,'al',lons=xlim,lats=ylim)
    albedo=albedo/100.
    sw_net=sw_down*(1-albedo)
  else:
    sw_net=sw_down-sw_up

  sw_net=np.where(sw_net<1.e-10,0,sw_net)
  out['radsw']=Data(x,y,sw_net,'W m-2',info='positive downward')

  # Net longwave flux  [W/m^2]
  if not quiet: print(' --> Net longwave flux')
  if not quiet: print('       LW down')
  if isPygrib:
    #x,y,lw_down = gribu.getvar(fname,'',tags='Downward long-wave radiation flux',lons=xlim,lats=ylim)
    x,y,lw_down = gribu.getvar(fname,'dlwrf',lons=xlim,lats=ylim)
  else:
    x,y,lw_down = gribu.getvar(fname,'downward long-wave',lons=xlim,lats=ylim)

  if not quiet: print('       LW up')
  x,y,lw_up   = gribu.getvar(fname,'ulwrf',tags='surface',lons=xlim,lats=ylim)
  if lw_up is False:
    if not quiet: print('       LW up not found: using sst')
    if isPygrib:
      #x,y,sst=gribu.getvar(fname,'temperature',tags='surface',lons=xlim,lats=ylim) # K
      x,y,sst=gribu.getvar(fname,'t',lons=xlim,lats=ylim) # K
    else:
      x,y,sst=gribu.getvar(fname,'temperature',tags='water surface',lons=xlim,lats=ylim) # K

    lw_net=air_sea.lwhf(sst,lw_down)
    lw_up=lw_down-lw_net
  else:
    lw_net=lw_down-lw_up

  # ROMS convention: positive upward
  # GFS convention: positive downward --> * (-1)
  lw_net=np.where(np.abs(lw_net)<1.e-10,0,lw_net)
  out['radlw']=Data(x,y,-lw_net,'W m-2',info='positive upward')

  # downward lw:
  out['dlwrf']=Data(x,y,-lw_down,'W m-2',info='negative... downward')


  # U and V wind speed 10m
  if not quiet: print(' --> U and V wind')
  #x,y,uwnd  = gribu.getvar(fname,'u',tags=':10 m',lons=xlim,lats=ylim)
  #x,y,vwnd  = gribu.getvar(fname,'v',tags=':10 m',lons=xlim,lats=ylim)
  x,y,uwnd  = gribu.getvar(fname,'10u',lons=xlim,lats=ylim)
  x,y,vwnd  = gribu.getvar(fname,'10v',lons=xlim,lats=ylim)

  if not quiet: print(' --> calc wind speed and stress')
  speed = np.sqrt(uwnd**2+vwnd**2)
  taux,tauy=air_sea.wind_stress(uwnd,vwnd)

  out['wspd']=Data(x,y,speed,'m s-1')
  out['uwnd']=Data(x,y,uwnd,'m s-1')
  out['vwnd']=Data(x,y,vwnd,'m s-1')
  out['sustr']=Data(x,y,taux,'Pa')
  out['svstr']=Data(x,y,tauy,'Pa')


  # Cloud cover [0--100 --> 0--1]:
  if not quiet: print(' --> Cloud cover')
  x,y,clouds  = gribu.getvar(fname,'tcc',tags='atmosphere:level 0',lons=xlim,lats=ylim)
  if clouds is False:
    if not quiet: print('CALC clouds from LW,TAIR,TSEA and RH')
    # first get sst (maybe already done to calc lw_up)
    try: sst
    except:
      if not quiet: print('  get TSEA')
      if isPygrib:
        #x,y,sst=gribu.getvar(fname,'temperature',tags='surface',lons=xlim,lats=ylim) # K
        x,y,sst=gribu.getvar(fname,'t',lons=xlim,lats=ylim) # K
      else:
        x,y,sst=gribu.getvar(fname,'temperature',tags='water surface',lons=xlim,lats=ylim) # K

    clouds=air_sea.cloud_fraction(lw_net,sst-273.15,tair,rhum,'net')
  else: clouds=clouds/100.

  out['cloud']=Data(x,y,clouds,'fraction (0--1)')

  return out
