'''
Circular basemap polar stereographic projections
Run test() to check

Todo: clip rivers, states and countries.
These can be easily clipped using the .clipping_path

mma may2016

Usage:
  a=SPstere()
  a.fillcontinents()
  a.drawparallels()
  a.drawmeridians()

Requirements:
- matplotlib, basemap
- okean, needed only to calculate the circles if parallels not provided
   in drawparallels
'''

from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as pl

class _pstere(Basemap):
  def __init__(self,pole,resolution='c',boundinglat=0,lon_0=0,**kargs):

    if pole=='north': prj='npstere'
    else: prj='spstere'

    Basemap.__init__(self,resolution=resolution,
                     boundinglat=boundinglat,
                     lon_0=lon_0,projection=prj,
                     **kargs)

    self._tasks_done=False


  def _tasks(self,ax=None):
    if not self._tasks_done:
      if not ax: ax=pl.gca()
      ax.axis('image')
      ax.set_frame_on(0)
      ax.set_xticks([])
      ax.set_yticks([])
      if self.projection[0]=='s':
        ax.invert_xaxis()
        ax.invert_yaxis()

      self._tasks_done=True
    else: pass


  def fillcontinents(self, color='#9887ff', lake_color='b', ax=None, zorder=None, alpha=None,
                     ocean_color='w',ocean_alpha=1,edge_color='k',coast_color='none'):

    if not ax: ax=pl.gca()
    if self.projection[0]=='s': Lat0=-90
    else: Lat0=90

    # boundary:
    xc,yc=self([0,0],[Lat0,self.boundinglat])
    xc0,yc0=xc[0],yc[0]
    r=np.sqrt((xc[1]-xc[0])**2+(yc[1]-yc[0])**2)
    circle = pl.Circle((xc0,yc0),r,facecolor=ocean_color, edgecolor=edge_color,clip_on=0,
             alpha=ocean_alpha,zorder=zorder)
    ax.add_patch(circle)

    for i in range(len(self.coastpolygons)):
      xy=zip(*self.coastpolygons[i])
      if self.coastpolygontypes[i]==1: cor=color
      else: cor=lake_color
      poly=pl.Polygon(xy,edgecolor=coast_color,facecolor=cor,alpha=alpha,zorder=zorder)
      axpoly=ax.add_patch(poly)
      axpoly.set_clip_path(circle)

    self._tasks(ax)

    self.clipping_path=circle


  def drawparallels(self,parallels='auto',color='k',linewidth=1., zorder=None,\
                      dashes=[1,1],ax=None,**kargs):

    if not ax: ax=pl.gca()
    if self.projection[0]=='s': Lat0,Lat1=-90,self.boundinglat
    else: Lat0,Lat1=self.boundinglat,90

    if parallels=='auto':
      try:
        from okean import ticks
        tks=ticks.tight(Lat0,Lat1,4)
      except ImportError:
        print 'Warning: okean needed to calculate parallels (circles)'
        print 'You may get okean from https://github.com/martalmeida/okean'
        print 'or specify them: ex: drawparallels([-80,-60])'
        print 'Will now use 2 linearly spaced ticks...'
        tks=np.linspace(Lat0,Lat1,4)[1:-1]
        print 'circles: ',tks

    else: tks=np.asarray(parallels)
    self.circles=tks

    if self.projection[0]=='s':
      self.central_lat=tks.min()
    else:
      self.central_lat=tks.max()

    if not self.boundinglat in tks: tks=np.hstack((tks,self.boundinglat))
    t=np.linspace(0,360,100)

    for a in tks:
      xx,yy=self(t,a+0*t)
      if a==self.boundinglat:
        ax.plot(xx,yy,color=color,clip_on=0,lw=linewidth,zorder=zorder)
      else:
        ax.plot(xx,yy,color=color,dashes=dashes,lw=linewidth,zorder=zorder)

    self._tasks(ax)


  def drawmeridians(self,meridians=range(0,360,45),color='k',linewidth=1., zorder=None,\
                      dashes=[1,1],fmt='%g',ax=None,**kargs):

    fontsize=kargs.get('fontsize',8)
    labelpad=kargs.get('labelpad',0.04)

    central_lat=kargs.get('central_lat',None)
    if not hasattr(self,'central_lat') and not central_lat:
      print 'Error: central_lat not set. Provide it in kargs of draw parallels first'
      return

    if not ax: ax=pl.gca()
    for a in meridians:
      xx,yy=self([a,a],[self.central_lat,self.boundinglat])
      ax.plot(xx,yy,color=color,dashes=dashes,zorder=zorder)

      if self.projection[0]=='s':
        xx,yy=self(a,self.boundinglat+np.abs(self.boundinglat)*labelpad)
      else:
        xx,yy=self(a,self.boundinglat-np.abs(self.boundinglat)*labelpad)

      A=a-self.projparams['lon_0']
      A=np.mod(A,360)

      if self.projection[0]=='s':
        if A>=90 and A<=270:
          ang=180-A
        else:
          ang=-A
      else:
        if A>=90 and A<=270:
          ang=-180+A
        else:
          ang=A

      if fmt:
        if a in [0,180]:
         lab=('$\mathrm{'+fmt+'^{\circ}}$')%a
        else:
          if a>180:
           lab=('$\mathrm{'+fmt+'^{\circ}W}$')%(360-a)
          else:
           lab=('$\mathrm{'+fmt+'^{\circ}E}$')%a

        ax.text(xx,yy,lab,rotation=ang,va='center',ha='center',fontsize=fontsize)


    self._tasks(ax)


class SPstere(_pstere):
  def __init__(self,resolution='c',boundinglat=-40,lon_0=135,**kargs):
    _pstere.__init__(self,'south',resolution,boundinglat,lon_0,**kargs)


class NPstere(_pstere):
  def __init__(self,resolution='c',boundinglat=40,lon_0=0,**kargs):
    _pstere.__init__(self,'north',resolution,boundinglat,lon_0,**kargs)


def test():
  '''Compare north and south polar stereographic projections with basemap'''

  pl.figure()
  ax=[pl.subplot(2,2,i) for i in range(1,5)]

  a=SPstere()
  a.fillcontinents(ax=ax[0])
  a.drawparallels(ax=ax[0])
  a.drawmeridians(ax=ax[0],labelpad=0.1)

  a=NPstere()
  a.fillcontinents(ax=ax[1])
  a.drawparallels(ax=ax[1])
  a.drawmeridians(ax=ax[1],labelpad=0.1)

  sp=Basemap(resolution='c',boundinglat=-40,lon_0=135,projection='spstere')

  sp.fillcontinents(ax=ax[2],color='#9887ff',lake_color='b')
  sp.drawparallels((-80,-60),ax=ax[2])
  sp.drawmeridians(range(0,360,45),ax=ax[2],labels=(1,1,1,1),fontsize=8)

  npl=Basemap(resolution='c',boundinglat=40,lon_0=0,projection='npstere')

  npl.fillcontinents(ax=ax[3],color='#9887ff',lake_color='b')
  npl.drawparallels((60,80),ax=ax[3])
  npl.drawmeridians(range(0,360,45),ax=ax[3],labels=(1,1,1,1),fontsize=8)
