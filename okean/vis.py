import pylab as pl
import numpy as np
from okean import dateu as dts, ticks, cache, calc
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

if 1:
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
    default['proj.name']='auto'#'merc'
    default['proj.resolution']='auto'#'i'
    default['proj.options']={}

    default['proj.coast_add']=True
    default['proj.coast']={'color': '#999999', 'linewidth': 0.5}

    default['proj.continents_add'] = True
    default['proj.continents']     = {'color': '#f2f2f2'}#ebebeb'}

    default['proj.oceans_add'] = False
    default['proj.oceans'] = {'color': 'aqua'}


    default['proj.countries_add'] = True
    default['proj.countries']     = {'color': '#d4d4d4'}

    default['proj.meridians_add'] = True
    default['proj.meridians_vals'] = 'auto'
    default['proj.meridians']     = {'labels':[0,0,0,1],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5}

    default['proj.parallels_add'] = True
    default['proj.parallels_vals'] = 'auto'
    default['proj.parallels']     = {'labels':[1,0,0,0],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5}

    default['proj.rivers_add'] = False
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
    default['vfield.C']  = None

    # about 1d plots:
    default['d1_fill.options']  = {'lw':1,'facecolor':'#cccccc','edgecolor':'none'}
    default['d1_line.options']  = {'lw':1}
    default['d1.plot']  = 'plot'
    default['d1.y0']    = 'min'

    param=default
    del(default)


class visCfg:

  def __init__(self):
    '''
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
    default['d1_fill.options']  = {'lw':1,'facecolor':'#cccccc','edgecolor':'none'}
    default['d1_line.options']  = {'lw':1}
    default['d1.plot']  = 'plot'
    default['d1.y0']    = 'min'
    '''

    #self.default=default
    #self.config=default.copy()
    self.config=param.copy()


  def get_param(self,name0,name1=False):
    keys=self.config.keys()
    res={}
    for k in keys:
      k0=k.split('.')[0]
      k1='.'.join(k.split('.')[1:])
      if k0==name0:
        res[k1]=self.config[k]

    if name1:
      return res.get(name1,None)
    else:
      return res


  def set_param(self,**kargs):
    quiet=kargs.pop('quiet',0)

    def find_key(key):
      keys=self.config.keys()
      found=[]
      for k in keys:
        K=k.split('.')[-1]
        if K==key or k==key: found+=[k]
      return found

    for K in kargs.keys():
      key=find_key(K.replace('__','.'))
      if len(key)>1 and not quiet: print 'found more than one matching key!',key
      else:
       key=key[0]
       self.config[key]=kargs[K]


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
        #print 'adding sec ',section
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

  def init_handles(self):
    self.handles={}
    for k in ['mappable','quiver','quiver_key','d1']: self.handles[k]=[]


  def init_figure(self,**kargs):
    fig=kargs.get('fig',None) # use existing fig; create a new one if True
    ax=kargs.get('ax',None) # use existing axes

    self.init_handles()

    prev=kargs.get('prev',None)
    if prev:
      print 'using prev !!'
      for i in ['fig','ax','cbax','cbbg']:
        setattr(self,i,getattr(prev,i))
      return

    if isinstance(fig,pl.Figure):
      self.fig=fig
    else:
      openfig=fig is True

      if not openfig:
        # check if figure exist:
        try:
          exists=self.fig.number in pl.get_fignums()
        except: exists=False
        openfig=not exists

      # figure:
      if openfig:
#        print 'NEW FIG!'
        self.fig=pl.figure(figsize=self.config['figure.figsize'],
                           facecolor=self.config['figure.facecolor'],
                           edgecolor=self.config['figure.edgecolor'])

    # main ax:
    if ax: self.ax=ax
    else:
#      print 'new ax'
      self.ax=self.fig.add_axes(self.config['axes.position'])

    # colorbar bg:
    if  self.config['colorbar.bg_position']:
      self.cbbg=self.fig.add_axes(self.config['colorbar.bg_position'])
      self.cbbg_rec=matplotlib.patches.Rectangle((0,0),1,1,**self.config['colorbar.bg'])
      self.cbbg.add_artist(self.cbbg_rec)
      self.cbbg.set(frame_on=0,xticks=[],yticks=[],visible=False)
    else: self.cbbg=False

    # colorbar ax:
    if self.config['colorbar.ax_position']:
      self.cbax=self.fig.add_axes(self.config['colorbar.ax_position'])
      self.cbax.set(visible=False)
    else: self.cbax=False


  def clear_figure(self):
    self.fig.clf()
    self.init_figure()

  def clear_axes(self):
    self.init_handles()
    try: self.ax.cla()
    except: pass
    try: self.cbax.cla()
    except: pass
    try: self.cbbg.cla()
    except: pass

  def map_info_get(self):
    res=dict(limits=self.config['axes.axis'],proj=self.config['proj.name'],res=self.config['proj.resolution'])
    res.update(self.config['proj.options'])
    return res


  def map_info_copy(self,prevObj):
    for k in ['axes.axis','proj.name','proj.resolution']:
      self.config[k]=prevObj.config[k]

    self.map_info_current=self.map_info_get()


  def init_projection(self,**kargs):
    debug_lev=kargs.get('debug_level',0)

    # proj xy limits:
    if not self.config['axes.axis']:
      xlim=-100,20
      ylim=-30,60
      self.config['axes.axis']=xlim+ylim
    else: xlim,ylim=self.config['axes.axis'][:2],self.config['axes.axis'][2:]

    # proj name:
    if self.config['proj.name']=='auto':
       cond0=ylim[1]>80
       cond1=ylim[0]<-80
       lonrange=xlim[1]-xlim[0]
       latrange=ylim[1]-ylim[0]

       name='cyl'
       #if cond0 and cond1: name='cyl'
       if lonrange<90 and (cond0 and not cond1 or (cond1 and not cond0)): name='lcc'
       elif not cond0 and not cond1: name='merc'
       elif lonrange>330 and ylim[1]<-30: name='spstere'
       elif lonrange>330 and ylim[0]>30: name='npstere'
       self.config['proj.name']=name


    # proj resolution:
    if self.config['proj.resolution']=='auto':
       lonrange=xlim[1]-xlim[0]
       latrange=ylim[1]-ylim[0]
       Amax=360*180.
       A=lonrange*latrange
       if A/Amax>0.5: resolution='c'
       elif A/Amax>0.1: resolution='l'
       elif A/Amax>0.0001: resolution='i'
       else: resolution='h'
       self.config['proj.resolution']=resolution

   
    # proj options: 
    opts=self.config['proj.options']
    if self.config['proj.name'] in ('spstere','npstere'):
      lon_0 = opts.get('lon_0',0)
      if self.config['proj.name']=='spstere': 
        blat  = opts.get('boundinglat',ylim[1])
      else:
        blat  = opts.get('boundinglat',ylim[0])

      opts['lon_0']       = lon_0
      opts['boundinglat'] = blat
    else:
      lon_0=opts.get('lon_0',0.5*(xlim[0]+xlim[1]))
      lat_0=opts.get('lat_0',0.5*(ylim[0]+ylim[1]))

      opts['lon_0'] = lon_0
      opts['lon_0'] = lat_0

    # load proj from cache or create new:
    self.map_info_current=self.map_info_get()
    cacheLab=str(self.map_info_current)
    cch=cache.Cache()

    if cch.is_stored(cacheLab,'localfile'):
      if debug_lev==2: print ' -> loading proj from cache'
      m=cch.load(cacheLab,'localfile')
    else:
      if debug_lev==2: print ' -> creating proj'

      if self.config['proj.name'] in ('spstere','npstere'):
        m = Basemap(projection=self.config['proj.name'],
                  resolution=self.config['proj.resolution'],
                  lon_0=lon_0,boundinglat=blat)
      else:
        m = Basemap(projection=self.config['proj.name'], lat_ts=0.0,
                  resolution=self.config['proj.resolution'],
                  urcrnrlon=xlim[1], urcrnrlat=ylim[1],
                  llcrnrlon=xlim[0], llcrnrlat=ylim[0],
                  lon_0=lon_0,lat_0=lat_0)

      cch.store(cacheLab,m,'localfile')

    self.map=m


  def draw_projection(self):
    #if not self.config['proj.name']:
    #  self.map=False
    #  return
    if not self.map: return

    if self.map_info_get()!=self.map_info_current:
      self.init_projection()

    # there is a bug in set_axes_limits, called at the end of drawcoastlines,
    # fillcontinents, etc, here:
    #     first draw boundary, no fill
    #     limb1 = self.drawmapboundary(fill_color='none')
    # argument ax is missing! So, boundary will go the the latest axes, which can be self.cbax!!
    # my current version is 1.0.7
    # To avoid it:
    pl.axes(self.ax)

    # better do this before draw anything else, otherwise the drawmapboundary will be
    # called many times (by set_axes_limits, depending on the projection used)
    if self.config['proj.oceans_add']:
      self.map.drawmapboundary(ax=self.ax,**self.config['proj.oceans'])

    if self.config['proj.coast_add']:
      self.map.drawcoastlines(ax=self.ax,**self.config['proj.coast'])


    if self.config['proj.continents_add']:
      self.map.fillcontinents(ax=self.ax,**self.config['proj.continents'])


    if self.config['proj.countries_add']:
      self.map.drawcountries(ax=self.ax,**self.config['proj.countries'])

    if self.config['proj.states_add']:
      self.map.drawstates(ax=self.ax,**self.config['proj.states'])

    if self.config['proj.rivers_add']:
      self.map.drawrivers(ax=self.ax,**self.config['proj.rivers'])

    if self.config['proj.meridians_add']:
      if self.config['proj.meridians_vals']=='auto':
        if self.map_info_current['proj'].endswith('stere'):
          meridians=range(0,360,45)
        else:
          xlim=self.map.llcrnrlon,self.map.urcrnrlon
          xlim=self.map.lonmin,self.map.lonmax
          meridians=ticks.loose_label(xlim[0],xlim[1])
      else: meridians=self.config['proj.meridians_vals']
      self.map.drawmeridians(meridians,ax=self.ax,**self.config['proj.meridians'])

    if self.config['proj.parallels_add']:
      if self.config['proj.parallels_vals']=='auto':
        ylim=self.map.llcrnrlat,self.map.urcrnrlat
        ylim=self.map.latmin,self.map.latmax
        parallels=ticks.loose_label(ylim[0],ylim[1])
      else: parllels=self.config['proj.parallels_vals']
      self.map.drawparallels(parallels,ax=self.ax,**self.config['proj.parallels'])


  def _convCoord(self,x,y):
    if x.size!=y.size: x,y=np.meshgrid(x,y)
    if hasattr(self,'map') and self.map: return self.map(x,y)
    else: return x,y


  def draw_scalar_field(self,x,y,v):

    if x is None or y is None:
      ny,nx=v.shape
      x,y=np.meshgrid(np.arange(nx,dtype=v.dtype),np.arange(ny,dtype=v.dtype))

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

    if x is None or y is None:
      ny,nx=v.shape
      x,y=np.meshgrid(np.arange(nx,dtype=v.dtype),np.arange(ny,dtype=v.dtype))

####    di,dj=self.config['vfield.dij']

    x,y=self._convCoord(x,y)

    m=np.zeros(x.shape,'bool')
    dij=self.config['vfield.dij']
    dij0=self.config['vfield.dij0']
    m[dij0[0]::dij[0],dij0[1]::dij[1]]=True

    #self.handles['quiver']+=[self.ax.quiver(x[m],y[m],u[m],v[m],**self.config['vfield.options'])]

    C=self.config['vfield.C']
    if C:
      if C is 'speed':
        C=np.sqrt(u[m]**2+v[m]**2)
      else: C=C[m]
      args=x[m],y[m],u[m],v[m],C
    else:
      args=x[m],y[m],u[m],v[m]

    self.handles['quiver']+=[self.ax.quiver(*args,**self.config['vfield.options'])]

    X,Y,U=self.config['vfield.key_XYU']
    if U:
      if U=='auto': U=2*ticks.nicenum(np.sqrt(u**2+v**2).mean(),1)
      if self.config['vfield.key_label']=='auto': Label=str(U)
      else: Label=self.config['vfield.key_label']

      self.handles['quiver_key']+=[self.ax.quiverkey(self.handles['quiver'][-1],X,Y,U,Label,
                                   **self.config['vfield.key_options'])]

  def draw_colorbar(self):
    if self.cbax and len(self.handles['mappable']):
      opts=self.config['colorbar.options']
      w,h=self.config['colorbar.ax_position'][2:]
      if w>h: opts['orientation']='horizontal'
      self.cb=pl.colorbar(self.handles['mappable'][-1],cax=self.cbax,**opts)

    if  self.cbbg: self.cbbg.set(visible=True)
    if  self.cbax: self.cbax.set(visible=True)


  def draw_1d(self,x=None,y=None):
    if x is None:
      x=np.arange(y.size,dtype=y.dtype)
#      args=self.config['d1_line.options']
#      self.handles['d1']+=self.ax.plot(y,**args)
#    else:
    if 1:
      plt=self.config['d1.plot']

      if  plt=='fill':
        args=self.config['d1_fill.options']
        x,y=self._convCoord(x,y)
        self.handles['d1']+=[self.ax.fill(x,y,**args)]

      elif  plt=='fill_between':
        args=self.config['d1_fill.options']

        y0=self.config['d1.y0']
        if y0=='min': Y0=self.ax.get_ylim()[0]
        elif y0=='max': Y0=self.ax.get_ylim()[1]
        else: Y0=y0

        x,y=self._convCoord(x,y) # should also convert y0 !! how!? no need for now...
        self.handles['d1']+=[self.ax.fill_between(x,y,Y0,**args)]

      elif plt=='plot':
        args=self.config['d1_line.options']
        x,y=self._convCoord(x,y)
        self.handles['d1']+=[self.ax.plot(x,y,**args)]


   # TODO: deal with x as datetimes... choose best xticks

  def draw_label(self,labType,labStr):
    if not labStr: return
    lab=self.ax.set(**{labType:labStr})[0]
    args=self.config['axes.labels']
    lab.set(**args)

  def draw_fill(self,x,y): pass

  def draw_between(self,x,y,y0): pass


class Data(Vis):
  def __init__(self,*args,**kargs):
    self.x=None
    self.y=None
    self.d=None
    self.z=None
    self.t=None
    self.tnum=None # used with time_series (2d times), no longer needed with recent versions of mpl
    self.v=None
    self.msg=''
    self.extra=[]
    self.alias=''

    if len(args)==1: self.v=args[0]
    if len(args)==2: self.x,self.v=args
    elif len(args)==3: self.x,self.y,self.v=args
    elif len(args)>3:
      print 'bad len(args) !!'
   

    for k in kargs.keys():
      setattr(self,k,kargs[k]) 

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

  def plot(self,coords='auto',proj='auto',labels=True,extras=True,isExtra=0,**kargs):
    debug_lev=kargs.pop('debug_level',0)
    self.set_param(**kargs)

    try:
      if isinstance(self.v,tuple):
        ndim=self.v[0].ndim
        vshape=self.v[0].shape
      else:
        ndim=self.v.ndim
        vshape=self.v.shape
    except:
      ndim=0
      vshape=()

    x,y,xname,yname=None,None,None,None
    if coords=='auto':
      if ndim==2:
        if   not self.x    is None and (self.x.ndim==2    or self.x.size==vshape[1]):    x,xname=self.x,'x'
        elif not self.d    is None and (self.d.ndim==2    or self.d.size==vshape[1]):    x,xname=self.d,'d'
        elif not self.t    is None and (self.t.ndim==2    or self.t.size==vshape[1]):    x,xname=self.t,'t'
        elif not self.tnum is None and (self.tnum.ndim==2 or self.tnum.size==vshape[1]): x,xname=self.tnum,'tnum'
        elif not self.y    is None and (self.y.ndim==2    or self.y.size==vshape[1]):    x,xname=self.y,'y'

        if   not self.y is None and (self.y.ndim==2 or self.y.size==vshape[0]):  y,yname=self.y,'y'
        elif not self.z is None and (self.z.ndim==2 or self.z.size==vshape[0]):  y,yname=self.z,'z'

      elif ndim==1:
        y,yname=self.v,'v'
        if   not self.d is None and self.d.ndim==1: x,xname=self.d,'d'
        elif not self.t is None and (calc.isarray(self.t) and self.t.ndim==1): x,xname=self.t,'t'
        elif not self.x is None and self.x.ndim==1: x,xname=self.x,'x'
        elif not self.y is None and self.y.ndim==1: x,xname=self.y,'y'
        elif not self.z is None and self.z.ndim==1:
          y,yname=self.z,'z'
          x,xname=self.v,'v'


    if ndim==2 and proj=='auto' and xname=='x' and yname=='y' and\
       np.all(x>=-360) and np.all(x<=360) and np.all(y>=-90) and np.all(y<=90): proj=True
    elif proj=='auto': proj=False

    if hasattr(self,'fig') and pl.fignum_exists(self.fig.number):
      if isExtra: pass
      else:
        if debug_lev==2: print ' -> will clear axes'
###        self.init_clear()
        self.clear_axes()
    else:
      if debug_lev==2: print ' -> will create new fig'
      self.init_figure()

    if proj and not isExtra:
      if not hasattr(self,'map') or self.map_info_get()!=self.map_info_current:
        # make a new projection:
        if not self.config['axes.axis'] and not x is None: self.config['axes.axis']=x.min(),x.max(),y.min(),y.max()
        if debug_lev==2: print ' -> new projection needed'
        self.init_projection(debug_level=debug_lev)
###        useProj=1

    # set labels:
    xlab=ylab=vlab=''
    if not xname is None:
      xunits=self.info[xname]['units']
      if xunits in (None,'Unk'): xlab = '%s'%self.info[xname]['name']
      else: xlab = '%s (%s)'%(self.info[xname]['name'],xunits)

    if not yname is None:
      yunits=self.info[yname]['units']
      if yunits in (None,'Unk'): ylab = '%s'%self.info[yname]['name']
      else: ylab = '%s (%s)'%(self.info[yname]['name'],yunits)

    vunits=self.info['v']['units']
    if vunits in (None,'Unk'): vlab = '%s'%self.info['v']['name']
    else: vlab = '%s (%s)'%(self.info['v']['name'],vunits)


    if ndim==2:
      if debug_lev==2: print ' -> will draw 2d'
      if isinstance(self.v,tuple):
        self.draw_vector_field(x,y,self.v[0],self.v[1])
      else:
        self.draw_scalar_field(x,y,self.v)
        if not isExtra: self.draw_colorbar()

      if labels:
        tit=vlab
        if not self.t is None:
          if isinstance(self.t,datetime.datetime):
            tit+=' $\\bullet$ %s'%self.t.strftime('%Y-%m-%d %H:%M')
          else:
            tit+=' $\\bullet$ %s to %s'%(self.t[0,0].strftime('%Y-%m-%d %H:%M'),self.t[0,-1].strftime('%Y-%m-%d %H:%M'))

        if not self.info['v']['slice'] in (None,'Unk'):
          tit+=' $\\bullet$ %s'%self.info['v']['slice']

        if not tit=='Unk': self.draw_label('title',tit)
        if not proj:###not useProj: 
          if not xlab=='Unk': self.draw_label('xlabel',xlab)
          if not ylab=='Unk': self.draw_label('ylabel',ylab)


    elif ndim==1:
      if debug_lev==2: print ' -> will draw 1d'
      self.draw_1d(x,y)
      if labels:
        if not xlab=='Unk': self.draw_label('xlabel',xlab)
        if yname=='v':
          if not vlab=='Unk': self.draw_label('title',vlab)
        else:
          if not ylab=='Unk': self.draw_label('ylabel',ylab)

    if proj and not isExtra:
      if debug_lev==2: print ' -> will draw projection', self.info['v']
      self.draw_projection()


    if self.extra and extras:
      for e in self.extra:
        e.fig=self.fig
        e.ax=self.ax
        if 0:
          e.handles=self.handles # do I need this line!??
        else: e.init_handles()
        if hasattr(self,'map'):
           e.map=self.map
           e.map_info_copy(self)

        if debug_lev==2: print ' -> will plot extras'
        e.plot(labels=False,isExtra=1)
