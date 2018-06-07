'''
Tools to deal with Telemac output files

Requires parserSELAFIN

mma, IEO-A Coruna, oct2016
'''

import pylab as pl
import numpy as np
import datetime
import parserSELAFIN


class Telemac():
  def __init__(self,fname,**kargs):
    '''
    fname: 3d or 2d file, or both (3d,2d)
      If both are provided, 2d instace will be inside self.d2

    kargs:
      zone, zone number. If set, lonxlat are calculated using utm
    '''

    zone=kargs.get('zone',None)

    if isinstance(fname,basestring):
      q=parserSELAFIN.SELAFIN(fname)
      d2=None
    else:
      q=parserSELAFIN.SELAFIN(fname[0])
      d2=Telemac(fname[1],**kargs)


    self.x=q.MESHX
    self.y=q.MESHY
    self.N=q.NPLAN

    self.q=q
    if self.N>1:
      self.d2=d2

    # load time:
    t0=datetime.datetime(q.DATETIME[0],q.DATETIME[1],q.DATETIME[2],q.DATETIME[3])
    t=q.tags['times']
    time=np.zeros(t.size,datetime.datetime)
    for i in range(time.size):
      time[i]=t0+datetime.timedelta(hours=t[i]/3600.)

    self.time=time

    # try to load bottom:
    if self.N==1:
      self.h=self.load('bottom',0,quiet=1)

    if not zone is None:
      self.convert_to_lonlat(zone)

  def convert_to_lonlat(self,zone_number,northern=False):
    import utm
    lon=np.zeros(self.x.size,dtype=self.x.dtype)
    lat=np.zeros(self.x.size,dtype=self.x.dtype)
    for i in range(self.x.size):
      lat[i],lon[i]=utm.to_latlon(self.x[i], self.y[i], zone_number=zone_number, northern=northern)

    self.lon=lon
    self.lat=lat

  def centres(self):
    '''
    Needed for area weighted averages
    '''
    from okean import calc
    elements = []
    for e in self.q.IKLE2:
      element = []
      for n in e:
        element.append((self.x[n],self.y[n]))

      elements.append(element)

    elements=np.asarray(elements)
    X=elements[:,:,0]
    Y=elements[:,:,1]
    n=len(X)
    Cx=np.zeros(n,'f')
    Cy=np.zeros(n,'f')
    A=np.zeros(n,'f')
    inp=np.zeros(n,'bool')
    for i in range(n):
      A[i]=calc.poly_area(X[i],Y[i])
      Cx[i],Cy[i]=calc.poly_centroid(X[i],Y[i])
      # check if centre is inside polygon:
      inp[i]=calc.inpolygon(Cx[i],Cy[i],X[i],Y[i])

    self.xc=Cx
    self.yc=Cy
    if any(~inp): A=np.ma.masked_where(~inp,A)
    self.area=A
    self.inp=inp


  def load(self,vname,it,k=-1,quiet=0):
    # find var ind:
    I=-1
    for i,v in enumerate(self.q.VARNAMES):
      if v.lower().find(vname)>=0:
        if not quiet: print('loading %s at %s'%(v,self.time[it]))
        I=i
    if I>=0:
      u=self.q.getVariablesAt(it,[I])
      u.shape=self.N,u.size/self.N
      return u[k]

  def load_all_times(self,vname,k):
    nt=self.time.size
    for it in range(nt):
      tmp=self.load(vname,it,k)
      if it==0:
        u=np.zeros((nt,tmp.size),dtype=tmp.dtype)

      u[it]=tmp

    return u

  def plt_domain(self,ax=None,**kargs):
    alpha=kargs.get('alpha',1)
    ec=kargs.get('ec','face')
    fc=kargs.get('fc','b')
    C=kargs.get('C',None)
    clim=kargs.get('clim',None)

    if ax is None:
      newAx=True
      ax=pl.gca()
    else: newAx=False

    elements = []
    for e in self.q.IKLE2:
      element = []
      for n in e:
        element.append((self.x[n],self.y[n]))

      elements.append(element)

    # Collections
    colection = pl.matplotlib.collections.PolyCollection(
              elements,facecolor=fc,edgecolor=ec,alpha=alpha)

    if not C is None: colection.set_array(C)
    if not clim is None: colection.set_clim(clim)

    ax.add_collection(colection)
    if newAx:
      ax.axis('auto')
      if not C is None: pl.colorbar(colection)

    #ax.plot(self.x,self.y,'r.')#,alpha=.2)
    return colection
