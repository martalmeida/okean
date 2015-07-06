import pylab as pl
import numpy as np
from okean import dateu as dts, ticks, cache
import matplotlib
from string import join as sjoin
import ConfigParser
from mpl_toolkits.basemap import Basemap
import os
import datetime
try:
  from collections import OrderedDict
except:
  from ordereddict import OrderedDict


class visCfg:

  def __init__(self):
    default={}
    # figure and axes:
    default['figure.figsize']=8,6
    default['figure.edgecolor']=None
    default['figure.facecolor']=None
    default['axes.position']=.1,.1,.75,.8
    default['axes.axis']=False # xlim ylim
    default['axes.labels']={'fontsize':10}
    default['colorbar.ax_position']=.875,.1,.03,.8
    default['colorbar.options']={}
    default['colorbar.bg_position']=False #.86,.08,.08,.84
    default['colorbar.bg']={'edgecolor':'none','facecolor':'white','alpha':0.5,'linewidth':0}

    # map projection:
    default['proj.name']='merc'
    default['proj.resolution']='i'

    default['proj.coast_add']=True
    default['proj.coast']={'color': '#999999', 'linewidth': 0.5}

    default['proj.continents_add'] = True
    default['proj.continents']     = {'color': '#ebebeb'}

    default['proj.countries_add'] = True
    default['proj.countries']     = {'color': '#ebebeb'}

    default['proj.meridians_add'] = True
    default['proj.meridians_vals'] = 'auto'
    default['proj.meridians']     = {'labels':[0,0,0,1],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5}

    default['proj.parallels_add'] = True
    default['proj.parallels_vals'] = 'auto'
    default['proj.parallels']     = {'labels':[1,0,0,0],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5}

    default['proj.rivers_add'] = True
    default['proj.rivers']     = {'color': 'blue', 'linewidth': 0.5}

    default['proj.states_add'] = False
    default['proj.states']     = {'color': 'green', 'linewidth': 0.5}

    # about scalar field:
    default['field.plot']  = 'pcolormesh'
    default['field.clim']  = False
    default['field.cvals'] = False
    default['field.cmap']  = False
    default['field.linecolors']  = False

    # about vector field:
    default['vfield.options']={'units':'x','scale':None}
    default['vfield.key_options']={'coordinates':'figure'}
    default['vfield.key_XYU']=0.05,0.05,'auto'
    default['vfield.key_label']='auto'
    default['vfield.dij0']  = 0,0
    default['vfield.dij']  = 1,1

    # about 1d plots:
    default['1d.options']  = {'lw':1,'facecolor':'#cccccc','edgecolor':'none'}
    default['1d.plot']  = 'plot'
    default['1d.y0']    = 'min'

    self.default=default
    self.config=default.copy()


  def set_param(self,**kargs):
    quiet=kargs.pop('quiet',0)
    for k in kargs.keys():
      if not self.config.has_key(k) and not quiet: print 'new key %s' % k
      self.config[k]=kargs[k]

  def __cp2config(self):
    for k in self.cp.sections():
      for o in self.cp.options(k):
        val=self.cp.get(k,o)
        try: val=eval(val)
        except: pass
        self.config['%s.%s'%(k,o)]=val


  def __config2cp(self):
    config=ConfigParser.RawConfigParser(dict_type=OrderedDict)


    keys=self.config.keys()
    keys.sort()
    print keys
    for k in keys:
      if k.find('.')>0:
         tmp=k.split('.')
         section=tmp[0]
         name=sjoin(tmp[1:],'')
      else:
         section='unknown'
         name=k

      if not config.has_section(section):
        config.add_section(section)
        print 'adding sec ',section
      config.set(section,name,self.config[k])

    self.cp=config

  def write_config(self,fname):
    self.__config2cp()
    self.cp.write(open(fname,'w'))

  def read_config(self,fname):
    if not hasattr(self,'cp'): self.__config2cp()
    self.cp.read(fname)
    self.__cp2config()


class Vis(visCfg):

  def init_figure(self,openfig=True):
    self.handles={}
    for k in ['mappable','quiver','quiver_key','1d']: self.handles[k]=[]

    if not openfig:
      # check if figure exist:
      try:
        exists=self.fig.number in pl.get_fignums()
      except: exists=False
      openfig=not exists

     # figure:
    if openfig:
      self.fig=pl.figure(figsize=self.config['figure.figsize'],
                         facecolor=self.config['figure.facecolor'],
                         edgecolor=self.config['figure.edgecolor'])

    # main ax:
    self.ax=self.fig.add_axes(self.config['axes.position'])

    # colorbar:
    if  self.config['colorbar.bg_position']:
      self.cbbg=self.fig.add_axes(self.config['colorbar.bg_position'])
      self.cbbg_rec=matplotlib.patches.Rectangle((0,0),1,1,**self.config['colorbar.bg'])
      self.cbbg.add_artist(self.cbbg_rec)
      self.cbbg.set(frame_on=0,xticks=[],yticks=[],visible=False)
    else: self.cbbg=False

    if self.config['colorbar.ax_position']:
      self.cbax=self.fig.add_axes(self.config['colorbar.ax_position'])
      self.cbax.set(visible=False)
    else: self.cbax=False


  def clear_figure(self):
    self.fig.clf()
    self.init_figure(openfig=False)

  def map_info_get(self):
    return dict(limits=self.config['axes.axis'],proj=self.config['proj.name'],res=self.config['proj.resolution'])


  def map_info_copy(self,prevObj):
    for k in ['axes.axis','proj.name','proj.resolution']:
      self.config[k]=prevObj.config[k]

    self.map_info_current=self.map_info_get()


  def init_projection(self):
    self.map_info_current=self.map_info_get()

    cacheLab=str(self.map_info_current)
    cch=cache.Cache()
    if cch.is_stored(cacheLab,'localfile'):
      m=cch.load(cacheLab,'localfile')
      print 'LOADING PROJ from cache !'
    else:
      print 'CREATING PROJ'
      if not self.config['axes.axis']:
        xlim=-100,20
        ylim=-30,60
      else: xlim,ylim=self.config['axes.axis'][:2],self.config['axes.axis'][2:]

      m = Basemap(projection=self.config['proj.name'], lat_ts=0.0,
                  resolution=self.config['proj.resolution'],
                  urcrnrlon=xlim[1], urcrnrlat=ylim[1],
                  llcrnrlon=xlim[0], llcrnrlat=ylim[0],
                  lon_0=0.5*(xlim[0]+xlim[1]),
                  lat_0=0.5*(ylim[0]+ylim[1]))

      cch.store(cacheLab,m,'localfile')

    self.map=m


  def draw_projection(self):
    if not self.config['proj.name']:
      self.map=False
      return

    if self.map_info_get()!=self.map_info_current:
      self.init_projection()

    if self.config['proj.coast_add']:
      self.map.drawcoastlines(ax=self.ax,**self.config['proj.coast'])

    if self.config['proj.continents_add']:
      self.map.fillcontinents(ax=self.ax,**self.config['proj.continents'])

    if self.config['proj.countries_add']:
      self.map.drawcountries(ax=self.ax,**self.config['proj.countries'])

    if self.config['proj.states_add']:
      self.map.drawstates(ax=self.ax,**self.config['proj.states'])

    print 'ADDING RIVERS !!! ',self.config['proj.rivers_add']
    if self.config['proj.rivers_add']:
      self.map.drawrivers(ax=self.ax,**self.config['proj.rivers'])

    if self.config['proj.meridians_add']:
      if self.config['proj.meridians_vals']=='auto':
        xlim=self.map.llcrnrlon,self.map.urcrnrlon
        meridians=ticks.loose_label(xlim[0],xlim[1])
      else: meridians=self.config['proj.meridians_vals']
      self.map.drawmeridians(meridians,ax=self.ax,**self.config['proj.meridians'])

    if self.config['proj.parallels_add']:
      if self.config['proj.parallels_vals']=='auto':
        ylim=self.map.llcrnrlat,self.map.urcrnrlat
        parallels=ticks.loose_label(ylim[0],ylim[1])
      else: parllels=self.config['proj.parallels_vals']
      self.map.drawparallels(parallels,ax=self.ax,**self.config['proj.parallels'])


  def _convCoord(self,x,y):
    if x.size!=y.size: x,y=np.meshgrid(x,y)
    if hasattr(self,'map') and self.map: return self.map(x,y)
    else: return x,y


  def draw_scalar_field(self,x,y,v):

    args={}
    if self.config['field.cmap']:
      if isinstance(self.config['field.cmap'],basestring): cmap=eval('pl.cm.'+self.config['field.cmap'])
      else: cmap=self.config['field.cmap']
      args['cmap']=cmap

    x,y=self._convCoord(x,y)
    meth=eval('self.ax.'+self.config['field.plot'])

    if self.config['field.clim']:
      vmin,vmax=self.config['field.clim']
      args['vmin']=vmin
      args['vmax']=vmax

    if self.config['field.linecolors']:
      args['colors']=self.config['field.linecolors']

    if self.config['field.plot'] in ('pcolor','pcolormesh'):
      self.handles['mappable']+=[meth(x,y,v,**args)]
    elif self.config['field.plot'] in ('contour','contourf'):
       if self.config['field.plot'].startswith('contour'):
         # contour(f) does not like bad values in coordinates, even when masked
         # so:
         if np.ma.isMA(x): x[x.mask]=x.mean()
         if np.ma.isMA(y): y[y.mask]=y.mean()

       if self.config['field.cvals'] is False: vals=20
       else: vals=self.config['field.cvals']
       self.handles['mappable']+=[meth(x,y,v,vals,**args)]


  def draw_vector_field(self,x,y,u,v):
    di,dj=self.config['vfield.dij']

    x,y=self._convCoord(x,y)

    m=np.zeros(x.shape,'bool')
    dij=self.config['vfield.dij']
    dij0=self.config['vfield.dij0']
    m[dij0[0]::dij[0],dij0[1]::dij[1]]=True

    self.handles['quiver']+=[self.ax.quiver(x[m],y[m],u[m],v[m],**self.config['vfield.options'])]

    X,Y,U=self.config['vfield.key_XYU']
    if U:
      if U=='auto': U=2*ticks.nicenum(np.sqrt(u**2+v**2).mean(),1)
      if self.config['vfield.key_label']=='auto': Label=str(U)
      else: Label=self.config['vfield.key_label']

      self.handles['quiver_key']+=[self.ax.quiverkey(self.handles['quiver'][-1],X,Y,U,Label,
                                   **self.config['vfield.key_options'])]

  def draw_colorbar(self):
    if self.cbax and len(self.handles['mappable']):
      self.cb=pl.colorbar(self.handles['mappable'][-1],cax=self.cbax,**self.config['colorbar.options'])

    if  self.cbbg: self.cbbg.set(visible=True)
    if  self.cbax: self.cbax.set(visible=True)


  def draw_1d(self,x=None,y=None):
    args=self.config['1d.options']
    if x is None:
      self.handles['1d']+=self.ax.plot(y,**args)
    else:
      plt=self.config['1d.plot']

      if  plt=='fill':
        self.handles['1d']+=[self.ax.fill(x,y,**args)]
      elif  plt=='fill_between':

        y0=self.config['1d.y0']
        if y0=='min': y0=self.ax.get_ylim()[0]
        elif y0=='max': y0=self.ax.get_ylim()[1]

        self.handles['1d']+=[self.ax.fill_between(x,y,y0,**args)]

   # TODO: deal with x as datetimes... choose best xticks

  def draw_label(self,labType,labStr):
    lab=self.ax.set(**{labType:labStr})[0]
    args=self.config['axes.labels']
    lab.set(**args)

  def draw_fill(self,x,y): pass

  def draw_between(self,x,y,y0): pass


class Data(Vis):
  def __init__(self):
    self.x=None
    self.y=None
    self.d=None
    self.z=None
    self.t=None
    self.tnum=None # used with time_series (2d times)
    self.v=None
    self.msg=''
#############################    self.plot=''
    self.extra=[]
    self.alias=''

    self.info={}
    for k in ('x','y','d','z','t','tnum','v'):
      self.info[k]=dict(name='Unk',units='Unk')

    self.info['v']['slice']='Unk'

    Vis.__init__(self)

  def show(self):
    for i in ['x','y','d','z','t','tnum']:
      a=getattr(self,i)
      try: print '%5s shape=%12s size=%d'%(i,str(a.shape),a.size)
      except: print '%5s %s'%(i,str(a))

    if not self.msg: msg='EMPTY'
    else: msg=self.msg
    print ' == msg == %s'%msg

  def plot(self,coords='auto',proj='auto',labels=True,extras=True,isExtra=0):
    x,y,xname,yname=None,None,None,None
    if coords=='auto':
      if self.v.ndim==2:
        if   not self.x    is None and self.x.ndim==2:    x,xname=self.x,'x'
        elif not self.d    is None and self.d.ndim==2:    x,xname=self.d,'d'
        elif not self.tnum is None and self.tnum.ndim==2: x,xname=self.tnum,'tnum'
        elif not self.y    is None and self.y.ndim==2:    x,xname=self.y,'y'

        if   not self.y is None and self.y.ndim==2:  y,yname=self.y,'y'
        elif not self.z is None and self.z.ndim==2:  y,yname=self.z,'z'

      elif self.v.ndim==1:
        if   not self.d is None and self.d.ndim==1: x,xname=self.d,'d'
        elif not self.t is None and self.t.ndim==1: x,xname=self.t,'t'
        elif not self.x is None and self.x.ndim==1: x,xname=self.x,'x'
        elif not self.y is None and self.y.ndim==1: x,xname=self.y,'y'
        elif not self.z is None and self.z.ndim==1: x,xname=self.z,'z'


    if self.v.ndim==2 and proj=='auto' and xname=='x' and yname=='y' and\
       np.all(x>=-360) and np.all(x<=360) and np.all(y>=-90) and np.all(y<=90): proj=True
    else: proj=False

#    # about x distance: --> better do this during slice
#    xUnitsAdd=''
#    if xname=='d' and (x.max()-x.min())>3e3:
#      x=x/1.e3
#      xUnitsAdd=r'$\mathrm{\times10^3}$ '


    if hasattr(self,'fig') and pl.fignum_exists(self.fig.number): pass
    else: self.init_figure()

#    try: self.fig
#    except:
#     print 'starting new fig !!'
#     self.init_figure()

###    useProj=0
    if proj and not isExtra:
      if not hasattr(self,'map') or self.map_info_get()!=self.map_info_current:
        # make a new projection:
        if not self.config['axes.axis']: self.config['axes.axis']=x.min(),x.max(),y.min(),y.max()
        print 'NEW PRJ !!'
        self.init_projection()
###        useProj=1

    if self.v.ndim==2:
      if isinstance(self.v,tuple):
        self.draw_vector_field(x,y,self.v[0],self.v[1])
      else:
        self.draw_scalar_field(x,y,self.v)
        if not isExtra: self.draw_colorbar()

      if labels:
        tit  = '%s (%s)'%(self.info['v']['name'],self.info['v']['units'])
        if not self.t is None: tit+=' $\\bullet$ %s'%self.t.strftime('%Y-%m-%d %H:%M')
        tit+=' $\\bullet$ %s'%self.info['v']['slice']

        self.draw_label('title',tit)
        if not proj:###not useProj: 
          #xlab = '%s (%s%s)'%(self.info[xname]['name'],xUnitsAdd,self.info[xname]['units'])
          xlab = '%s (%s)'%(self.info[xname]['name'],self.info[xname]['units'])
          ylab = '%s (%s)'%(self.info[yname]['name'],self.info[yname]['units'])
          self.draw_label('xlabel',xlab)
          self.draw_label('ylabel',ylab)

      if proj and not isExtra:
        print 'WILL DRAW PROJECTION', self.info['v']
        self.draw_projection()

    else:
      self.draw_1d(x,self.v)
      if labels:
        xlab = '%s (%s)'%(self.info[xname]['name'],xUnitsAdd,self.info[xname]['units'])
        tit  = '%s (%s)'%(self.info['v']['name'],self.info['v']['units'])
        self.draw_label('xlabel',xlab)
        self.draw_label('title',tit)


    if self.extra and extras:
      for e in self.extra:
        e.fig=self.fig
        e.ax=self.ax
        e.handles=self.handles
        if hasattr(self,'map'):
           e.map=self.map
           e.map_info_copy(self)

        print 'will plot extras...'
        e.plot(labels=False,isExtra=1)
        print e.config['1d.plot']
