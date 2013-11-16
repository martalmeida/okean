import ftplib
import os
import pygrib
import datetime
import numpy as np
import pylab as pl
from okean import calc, dateu, roms, netcdf
import os


Server='polar.ncep.noaa.gov'
Path='pub/history/waves'


def set_files(year,month):
  vars='hs','tp','dp'
  res={}
  for v in vars: res[v]='multi_1.glo_30m.%s.%d%02d.grb2'%(v,year,month)
  return res


class WW3Download:
  def __init__(self,destination='WW3_data',**kargs):
    self.destination=destination

    self.server=Server
    self.path=Path

    if 'server' in kargs.keys(): self.server=kargs['server']
    if 'path' in kargs.keys(): self.path=kargs['path']

  def connect(self):
    print 'connecting ...'
    self.ftp=ftplib.FTP(self.server)
    self.ftp.login()
    self.ftp.cwd(self.path)
    self.contents=self.ftp.nlst()
    print 'done'

  def download(self,year,month,overwrite=False):
    files=set_files(year,month).values()
    self.connect()  

    pdest=os.path.join(self.destination,str(year))
    if not os.path.isdir(pdest):
      print 'creating folder %s'%pdest
      os.makedirs(pdest)

    for fname in files:
      dest=os.path.join(pdest,fname)
      if os.path.isfile(dest) and not overwrite:
         print ' -- %s already exists'%dest
         continue
        
      print 'downloading %s'%dest
      f=open(dest,'w')
      self.ftp.retrbinary('RETR %s' % fname, f.write)
      f.close()


  def files(self,year,month):
    f=set_files(year,month).values()
    pdest=os.path.join(self.destination,str(year))
    return [os.path.join(pdest,i) for i in f]

  def files2(self,date0,date1):
    ym=dateu.mrange(date0.year,date0.month,date1.year,date1.month)
    res={}
    for d in ym:
      pdest=os.path.join(self.destination,str(d.year))
      tmp=set_files(d.year,d.month)
      for k in tmp.keys():
        f=os.path.join(pdest,tmp[k])
        if res.has_key(k): res[k]+=[f]
        else: res[k]=[f]

    return res
      

class WW3Data:
  def __init__(self,datafolder):
    self.datafolder=datafolder

  def data(self,date0,date1,xlim=False,ylim=False,quiet=False):

    Time=False
    Res={}

    files=WW3Download().files2(date0,date1)
    for k in files.keys():
      Res[k]=False
      for f in files[k]:

        if not os.path.isfile(f):
          print 'ERROR: missing data file %s'%f
          return

        lon,lat,time,val=ww3_file_data(f,xlim,ylim,quiet)

        cond=(time>=date0)&(time<=date1)
        time=time[cond]
        val=val[cond,...]
        if time.size:
          if Res[k] is False:
            Res[k]=val
            Time=time
          else:
            cond=time>Time[-1]
            if cond[cond].size:
              Res[k]=np.ma.vstack((Res[k],val[cond,...]))
              Time=np.hstack((Time,time[cond]))
        
    return lon,lat,Time,Res    


def ww3_file_data(fname,xlim=False,ylim=False,quiet=False):

  # fast load aux tmp file:
  xylab=''
  if xlim: xylab+='_%.2f_%.2f'%xlim
  if ylim: xylab+='_%.2f_%.2f'%ylim
  faux=fname+xylab+'.npz'
  if os.path.isfile(faux):
    if not quiet: print 'loading from %s'%faux
    a=np.load(faux)
    val=a['val']
    mval=a['mval']
    return a['lon'],a['lat'],a['time'],np.ma.masked_where(mval,val)


  f=pygrib.open(fname)

  nt=f.messages
  time=np.zeros(nt,datetime.datetime)
  
  for i in range(nt):
    o=f.message(i+1)
    v=o.values
    if not quiet and i%10==0: print '%03d of %03d %s'%(i,nt,o.validDate.isoformat(' '))
    if i==0:
      lat,lon=o.latlons()
      lon[lon>180]=lon[lon>180]-360

      if not xlim: lons=lon.min(),lon.max()
      if not ylim: lats=lat.min(),lat.max()
      i1,i2,j1,j2=calc.ij_limits(lon,lat,xlim,ylim,margin=1)
      lon=lon[j1:j2,i1:i2]
      lat=lat[j1:j2,i1:i2]
      ny,nx=v[j1:j2,i1:i2].shape
      val=np.ma.zeros((nt,ny,nx),'f')

    time[i]=o.validDate
    val[i,...]=v[j1:j2,i1:i2]

  # save aux tmp file for fast loading:
  if not quiet: print 'saving aux file %s'%faux
  np.savez(faux,lon=lon,lat=lat,time=time,val=val.data,mval=val.mask)

  return lon,lat,time,val


def ww3_plot(fname,xlim=False,ylim=False,res='c'):
  '''
  ex.:
    f='WW3_data/2012/multi_1.glo_30m.dp.201201.grb2'
    swanu.ww3_plot(f,xlim=[-100,-80],ylim=[20,32],res='i')
  '''

  from mpl_toolkits.basemap import Basemap
  import pylab as pl

  f=pygrib.open(fname)
  lat,lon=f.message(1).latlons()
  lon[lon>180]=lon[lon>180]-360.

  if not xlim: xlim=-180,180
  if not ylim: ylim=-90,90

  m = Basemap(projection='cyl',llcrnrlat=ylim[0],urcrnrlat=ylim[1],
              llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution=res)

  m.drawcoastlines()
  m.fillcontinents(color='#cccccc',lake_color='w')
  # draw parallels and meridians.
  m.drawparallels(np.arange(-90.,91.,30.))
  m.drawmeridians(np.arange(-180.,181.,60.))

  x,y=m(lon,lat)
  pl.plot(x,y,'b.',ms=.5)
  pl.show()


def in_comp_grid(XY,IJ,add=(0.001,0.001),plot=False):
  # got to be sure the point is inside the comp grid avoiding the message:
  # "** Error            : Boundary point outside comp. grid"

  i0=IJ[:,0]
  j0=IJ[:,1]
  i1=IJ[:,2]
  j1=IJ[:,3]

  x0=XY[:,0]
  y0=XY[:,1]
  x1=XY[:,2]
  y1=XY[:,3]

  M=IJ[:,::2].max() # nx
  N=IJ[:,1::2].max()

  if plot:
    pl.figure()
    pl.plot(x0,y0,'go')
    pl.plot(x1,y1,'go')

  addx,addy=add
  for i in range(len(i0)):
    if i0[i]==0 and i1[i]==0: # west
      x0[i]=x0[i]+addx
      x1[i]=x1[i]+addx
      if   j0[i]==0: y0[i]=y0[i]+addy
      elif j1[i]==N: y1[i]=y1[i]-addy

    elif i1[i]==M and i1[i]==M: # east
      x0[i]=x0[i]-addx
      x1[i]=x1[i]-addx
      if   j0[i]==0: y0[i]=y0[i]+addy
      elif j1[i]==N: y1[i]=y1[i]-addy

    if j0[i]==0 and j1[i]==0: # south
      y0[i]=y0[i]+addy
      y1[i]=y1[i]+addy
      if   i0[i]==0: x0[i]=x0[i]+addx
      elif i1[i]==M: x1[i]=x1[i]-addx

    if j0[i]==N and j1[i]==N: # north
      y0[i]=y0[i]-addy
      y1[i]=y1[i]-addy
      if   i0[i]==0: x0[i]=x0[i]+addx
      elif i1[i]==M: x1[i]=x1[i]-addx

  if plot:
    pl.plot(x0,y0,'k.')
    pl.plot(x1,y1,'kx')

  return XY


def ww3_specpoints(romsgrd,nseg=(5,5),addxy=(.001,.001)):
  gx=netcdf.use(romsgrd,'lon_rho')
  gy=netcdf.use(romsgrd,'lat_rho')
  mask=netcdf.use(romsgrd,'mask_rho')
  eta,xi=gx.shape

  if 0: # use range (step=nseg)
    spec_res=nseg
    dxi,deta=spec_res

    ix=range(0,xi,dxi)+[xi-1]
    iy=range(0,eta,deta)+[eta-1]
  else: # use linspace (nseg=n segments)
    nNS,nEW=nseg
    ix=np.round(np.linspace(0,xi-1,nNS)).astype('i')
    iy=np.round(np.linspace(0,eta-1,nEW)).astype('i')


  ix,iy=np.meshgrid(ix,iy)
  ix=calc.var_border(ix)
  iy=calc.var_border(iy)

  # unsorted unique:
  _,i=np.unique(ix+iy*1j,return_index=True)
  i=np.sort(i)

  ix=ix[i]
  iy=iy[i]

  # add 1st:
  ix=np.append(ix,ix[0])
  iy=np.append(iy,iy[0])

  # create segments:
  segx=[]
  segy=[]
  for i in range(len(ix)-1):
    I=ix[i],ix[i+1]
    J=iy[i],iy[i+1]
    i0,i1=np.min(I),np.max(I)
    j0,j1=np.min(J),np.max(J)

    if i0==i1: mseg=mask[j0:j1,i0]
    else:      mseg=mask[j0,i0:i1]


    if 1:
      # not use fully masked segments:
      if np.all(mseg==0):
        print 'masked segment %d %d %d %d'%(i0,j0,i1,j1)
        continue
    else:
      # not use if segment starts with mask:
      if mseg.size and mseg[0]==0:
        print 'masked 1st point of segment %d %d %d %d'%(i0,j0,i1,j1)
        continue

   
    segx+=[[i0,i1]]
    segy+=[[j0,j1]]

  # XY and IJ:
  XY=[]
  IJ=[]
  for i in range(len(segx)):
    i0,i1=segx[i][0],segx[i][1]
    j0,j1=segy[i][0],segy[i][1]

    IJ+=[[i0,j0,i1,j1]]
    XY+=[[gx[j0,i0], gy[j0,i0], gx[j1,i1], gy[j1,i1]]]

  IJ=np.asarray(IJ)
  XY=np.asarray(XY)

  # got to be sure the point is inside the comp grid avoiding the message:
  # "** Error            : Boundary point outside comp. grid"
  if not addxy is False:
    XY=in_comp_grid(XY,IJ,addxy)

  return XY,IJ


def plot_segments(grd,nseg=(10,10),addxy=(.001,.001)):
  g=roms.Grid(grd)
  xb,yb=g.border()

  pl.figure()
  pl.plot(xb,yb,'b.') 

  xy,ij=ww3_specpoints(grd,nseg,addxy)

  x=xy[:,0]
  y=xy[:,1]

  xx=xy[:,2]
  yy=xy[:,3]

  pl.plot(x,y,'o')
  pl.plot(xx,yy,'.')

  for i in range(len(x)):
    #print [x[i],xx[i]],[y[i],yy[i]]
    pl.plot([x[i],xx[i]],[y[i],yy[i]],'r')
    raw_input()



def ww3gb_2TPAR(datafolder,date0,date1,romsgrid,nseg=(10,10),addxy=(0.001,0.001),dspr=20,path='',name='TPAR'):
  '''dspr, directional spreading; see https://www.myroms.org/forum/viewtopic.php?f=14&t=2335
  '''

  g=roms.Grid(romsgrid)

  # loading data:
  xlim=g.lon.min()-2,g.lon.max()+2
  ylim=g.lat.min()-2,g.lat.max()+2

  lon,lat,Time,Res=WW3Data(datafolder).data(date0,date1,xlim,ylim)

  # calc spec points:
  xy,ij=ww3_specpoints(romsgrid,nseg,addxy)

  # interpolate data at initial point of each segment:
  x=xy[:,0]
  y=xy[:,1]

  
  nt=len(Time)
  nv=len(Res)
  nx=xy.shape[0]

  TPAR=np.zeros((nt,nv,nx),'f')
  vnames='hs','tp','dp'
  # Significant_height_of_combined_wind_waves_and_swell
  # Primary_wave_mean_period
  # Primary_wave_direction

  for it in range(nt):
    if it%10==0: print 'time %03d of %03d : %s'%(it,nt,Time[it].isoformat(' '))
    for iv in range(nv):
      vname=vnames[iv]
      #print '  -- %s'%vname
      v=Res[vname][it,...]

      # 1st extrap masked data:
      v=calc.mask_extrap(lon,lat,v)

      TPAR[it,iv,:]=calc.griddata(lon,lat,v,x,y)

  # create TPAR files:
  for i in range(nx):
    fname=name+'_seg%03d.txt'%i
    j=open(fname,'w')
    j.write('TPAR\n')
    for it in range(nt):

      # time as YYYYMMDD.HHMMSS 
      tiso=Time[it].strftime('%Y%m%d.%H%M%S')

      j.write('%s %6.2f %6.2f  %6.2f  %6.2f\n'%(tiso,TPAR[it,0,i],TPAR[it,1,i],TPAR[it,2,i],dspr))

    j.close()
    print 'created tpar file %s'%fname

  x1=xy[:,2]
  y1=xy[:,3]

  # add to swan INPUT:
  # BOUND SHAPESPEC JONSWAP PEAK DSPR DEGREES
  # BOUNDSPEC SEGMENT XY -54.3906 33.1171 -55.5271 35.8094 VARIABLE FILE 0 './forcings/TPAR2.txt'
  print '\n'
  print 'BOUND SHAPESPEC JONSWAP PEAK DSPR DEGREES'
  for i in range(nx):
    fname=os.path.join(path,name+'_seg%03d.txt'%i)
    print 'BOUNDSPEC SEGMENT XY %8.4f %8.4f %8.4f %8.4f VARIABLE FILE 0 \'%s\''%(x[i],y[i],x1[i],y1[i],fname)

  print '\n'


  i0=ij[:,0]
  j0=ij[:,1]
  i1=ij[:,2]
  j1=ij[:,3]
  for i in range(nx):
    fname=os.path.join(path,name+'_seg%03d.txt'%i)
    print 'BOUNDSPEC SEGMENT IJ %3d %3d %3d %3d VARIABLE FILE 0 \'%s\''%(i0[i],j0[i],i1[i],j1[i],fname)

  print '\n'


def roms2swan_grid(grd,exc=9999.,label='swan_bathy'):
  '''exc (EXCeption) = mask points
  '''

  x=netcdf.use(grd,'lon_rho')
  y=netcdf.use(grd,'lat_rho')
  m=netcdf.use(grd,'mask_rho')
  h=netcdf.use(grd,'h')
  h[m==0]=exc
  ny,nx=h.shape

  fbot=label+'.bot'
  fgrd=label+'.grd'

  # depths file:
  fb=open(fbot,'w')
  for i in range(nx):
    for j in range(ny):
      fb.write('   ')
      fb.write('%12.8f'%h[j,i])

    fb.write('\n')

  fb.close()
  print ' -- created swan depths file %s'%fbot

  # coords file:
  #v=np.hstack((x.T.flatten(),y.T.flatten()))
  v=np.hstack((x.flatten(),y.flatten()))
  np.savetxt(fgrd,v,fmt='%12.6f')
  print ' -- created swan computational grid file %s\n'%fgrd

  # add to swan INPUT:
  print '\n'
  print 'CGRID CURVILINEAR %d %d EXC 9.999000e+003 9.999000e+003 CIRCLE 36 0.04 1.0 24'%(nx-1,ny-1)
  print 'READGRID COORDINATES 1 \'%s\' 4 0 0 FREE '%(fgrd)

  print 'INPGRID BOTTOM CURVILINEAR 0 0 %d %d EXC 9.999000e+003'%(nx-1,ny-1)
  print 'READINP BOTTOM 1 \'%s\' 4 0 FREE '%(fbot)
  print '\n'


def roms2swan_wind(frc,date0,date1,fname='swan_wind.dat',**kargs):
  tname='wind_time'
  uname='Uwind'
  vname='Vwind'
  grd=False # needed if wind is 1d
  dt=1 # hours
  path=''
  if 'tname' in kargs.keys(): tname=kargs['tname']
  if 'uname' in kargs.keys(): uname=kargs['uname']
  if 'vname' in kargs.keys(): vname=kargs['vname']
  if 'grd'   in kargs.keys(): grd  =kargs['grd']
  if 'dt'    in kargs.keys(): dt   =kargs['dt']
  if 'path'  in kargs.keys(): path =kargs['path']
 
  print 'wind: loading time ...' 
  
  time=netcdf.nctime(frc,tname)
  #time=np.load('tfile')
  #cond=(time>=date0)&(time<=date1)
  cond=(time>=date0)&(time<=date1+datetime.timedelta(days=1)) # add one day at the end, just to avoid the "repeating last"
  time=time[cond]
  d=np.diff(pl.date2num(time))
  print 'current max and min dt = %6.2f %6.2f hrs = %6.2f %6.2f mins'%(d.max()*24, d.min()*24, d.max()*24*60, d.min()*24*60)
#  time=time[::dt]
#  d=np.diff(pl.date2num(time))
#  print ' final  max and min dt = %6.2f %6.2f hrs = %6.2f %6.2f mins'%(d.max()*24, d.min()*24, d.max()*24*60, d.min()*24*60)

  print 'wind: loading u ...' 
  u,nc=netcdf.var(frc,uname)
  print 'wind: loading v ...' 
  v,nc=netcdf.var(nc,uname)
#  u=u[cond,...][::dt,...]
#  v=v[cond,...][::dt,...]
  u=u[cond,...]
  v=v[cond,...]
  nc.close()


  if u.ndim==1:
    if not grd:
      print 'must provide grd karg!'
      return

    nt=u.size
    eta=netcdf.fdim(grd)['eta_rho']
    xi=netcdf.fdim(grd)['xi_rho']
  else:
    nt,eta,xi=u.shape

# array may be too big, so do this later (for each it)
#
#    u=np.tile(u[:,np.newaxis,np.newaxis],(1,eta,xi))
#    v=np.tile(v[:,np.newaxis,np.newaxis],(1,eta,xi))



  i=open(fname,'w')

  times=[]
  time0=time[0]-datetime.timedelta(hours=dt)
  ITs=[]
  for it in range(nt):
    time0=time0+datetime.timedelta(hours=dt)
    if time0>date1: break

    if time0>time[-1]:
      print 'Warning : repeating last ...', it

    times+=[time0]
    d=np.abs(time-time0)
    it=np.where(d==d.min())[0][0]
    ITs+=[it]

    if it%100==0: print 'saving u %s %s'%(fname,time[it].isoformat(' '))
    if u[it,...].ndim==0:
      U=np.tile(u[it,...],(eta,xi)).flatten()
    else:
      U=u[it,...].flatten()

    [i.write('%8.4f\n'%uu) for uu in U]

  for it in ITs:
    if it%100==0: print 'saving v %s %s'%(fname,time[it].isoformat(' '))
    if v[it,...].ndim==0:
      V=np.tile(v[it,...],(eta,xi)).flatten()
    else:
      V=v[it,...].flatten()

    [i.write('%8.4f\n'%vv) for vv in V]



  times=np.asarray(times)
  t0iso=times[0].strftime('%Y%m%d.%H%M%S')
  t1iso=times[-1].strftime('%Y%m%d.%H%M%S')
  dt=times[1]-times[0]
  dth=dt.days*24. + dt.seconds/60.**2

  print ' -- created swan wind file %s\n'%fname

  # add to swan INPUT:
  print '\n'
  print 'INPGRID WIND CURVILINEAR 0 0 %d %d EXC 9.999000e+003 &'%(xi-1,eta-1)
  print '       NONSTATIONARY %s %.2f HR %s'%(t0iso,dth,t1iso)
  print 'READINP WIND 1 \'%s\' 4 0 FREE '%(os.path.join(path,fname))
  print '\n'

