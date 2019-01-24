import os
import datetime
import pylab as pl
import numpy as np
from configparser import ConfigParser
from collections import OrderedDict
from mpl_toolkits.basemap import Basemap
from . import dateu as dts, ticks, cache, calc, cookbook as cb


param={}
param['plot.label']=''
param['plot.zorder']=1

# figure and axes:
param['figure.figsize']=None#8,6
param['figure.edgecolor']=None
param['figure.facecolor']=None
param['axes.position']=.1,.1,.75,.8
param['axes.axis']=False # xlim ylim
param['axes.labels']={'fontsize':10}
param['colorbar.ax_position']=.875,.1,.03,.8
param['colorbar.options']={}
param['colorbar.bg_position']=False #.86,.08,.08,.84
param['colorbar.bg']={'edgecolor':'none','facecolor':'white','alpha':0.5,'linewidth':0}

# map projection:
param['proj.name']='auto'
param['proj.resolution']='auto'
param['proj.options']={}

param['proj.coast_add']=True
param['proj.coast']={'color': '#999999', 'linewidth': 0.5, 'zorder':999}

param['proj.continents_add'] = True
param['proj.continents']     = {'color': '#f2f2f2', 'zorder':999}

param['proj.oceans_add'] = False
param['proj.oceans'] = {'color': 'aqua', 'zorder':999}

param['proj.countries_add'] = True
param['proj.countries']     = {'color': '#d4d4d4', 'zorder':999}

param['proj.meridians_add'] = True
param['proj.meridians_vals'] = 'auto'
param['proj.meridians']     = {'labels':[0,0,0,1],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5, 'zorder':999}

param['proj.parallels_add'] = True
param['proj.parallels_vals'] = 'auto'
param['proj.parallels']     = {'labels':[1,0,0,0],'color': '#4d4d4d', 'dashes': (1, 4), 'linewidth': 0.5, 'zorder':999}

param['proj.rivers_add'] = False
param['proj.rivers']     = {'color': 'blue', 'linewidth': 0.5, 'zorder':999}

param['proj.states_add'] = False
param['proj.states']     = {'color': 'green', 'linewidth': 0.5, 'zorder':999}

# about scalar field:
param['field.plot']  = 'pcolormesh'
param['field.clim']  = False
param['field.cvals'] = False
param['field.cmap']  = False
param['field.extend']  = None
param['field.linewidths']  = None

# about vector field:
param['vfield.options']={'units':'x','scale':None}
param['vfield.key_options']={'coordinates':'figure'}
param['vfield.key_XYU']=0.05,0.05,'auto'
param['vfield.key_label']='auto'
param['vfield.dij0']  = 0,0
param['vfield.dij']  = 5,5
param['vfield.C']  = None

# about 1d plots:
param['d1_fill.options']  = {'lw':1,'facecolor':'#cccccc','edgecolor':'none'}
param['d1_line.options']  = {'lw':1}
param['d1.plot']  = 'plot'
param['d1.y0']    = 'min'


class visCfg:

  def __init__(self):
    #self.config=param.copy()
    # this is a shallow copy! It is fine for param with integers and
    # string, but not ok for a param with a dict! We need a deepcopy:
    import copy
    self.config=copy.deepcopy(param)

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
      if len(key)>1 and not quiet: print('found more than one matching key!',key)
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
         name=''.join(tmp[1:])
      else:
         section='unknown'
         name=k

      if not config.has_section(section):
        config.add_section(section)
        #print('adding sec ',section)
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
    cbax=kargs.get('cbax',None) # use existing axes for colorbar
    if ax: fig=ax.figure

    self.init_handles()

    prev=kargs.get('prev',None)
    if prev:
      print('using prev !!')
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
        self.fig=pl.figure(figsize=self.config['figure.figsize'],
                           facecolor=self.config['figure.facecolor'],
                           edgecolor=self.config['figure.edgecolor'])

    # main ax:
    if ax: self.ax=ax
    else:
      self.ax=self.fig.add_axes(self.config['axes.position'])

    # colorbar bg:
    if  self.config['colorbar.bg_position']:
      self.cbbg=self.fig.add_axes(self.config['colorbar.bg_position'])
      self.cbbg_rec=pl.matplotlib.patches.Rectangle((0,0),1,1,**self.config['colorbar.bg'])
      self.cbbg.add_artist(self.cbbg_rec)
      self.cbbg.set(frame_on=0,xticks=[],yticks=[],visible=False)
    else: self.cbbg=False

    # colorbar ax:
    if cbax: self.cbax=cbax
    else:
      if self.config['colorbar.ax_position']:
        self.cbax=self.fig.add_axes(self.config['colorbar.ax_position'])
        self.cbax.set(visible=False)
      else: self.cbax=False


  def inherit(self,parent):
    self.fig=parent.fig
    self.ax=parent.ax
    self.init_handles()
    if hasattr(parent,'map'):
#       self.map=parent.map
#       self.map_info_copy(parent)
       self.copy_projection(parent)

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


#  def map_info_copy(self,prevObj):
#    for k in ['axes.axis','proj.name','proj.resolution','proj.options']:
#      self.config[k]=prevObj.config[k]
#
#    #self.map_info_current=self.map_info_get()
#    self.map_info_current=prevObj.map_info_get()


  def copy_projection(self,prevObj):
    self.map=prevObj.map
    for k in ['axes.axis','proj.name','proj.resolution','proj.options']:
      self.config[k]=prevObj.config[k]

    self.map_info_current=prevObj.map_info_get()


  def set_projection(self,opts):
    '''initiates the projection using basemap options provided as dict'''
    if not opts: return

    self.config['proj.options']={}
    for name in opts:
      val=opts[name]
      if name=='projection': self.config['proj.name']=val
      elif name=='resolution': self.config['proj.resolution']=val
      else: self.config['proj.options'][name]=val

    self.init_projection()


  def rm_projection(self):
    self.config['proj.name']='auto'
    self.config['proj.resolution']='auto'
    self.config['proj.options']={}
    self.config['axes.axis']=False
    self.map=False


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
      opts['lon_0'] = opts.get('lon_0',0)
      if self.config['proj.name']=='spstere':
        opts['boundinglat'] = opts.get('boundinglat',ylim[1])
      else:
        opts['boundinglat'] = opts.get('boundinglat',ylim[0])

    elif self.config['proj.name'] == 'stere':
      opts['lon_0']=opts.get('lon_0',0.5*(xlim[0]+xlim[1]))
      opts['lat_0']=opts.get('lat_0',0.5*(ylim[0]+ylim[1]))
      opts['lat_ts']=opts.get('lat_ts',0.5*(ylim[0]+ylim[1]))

    elif self.config['proj.name'] in ('merc','tmerc','cyl'):
      opts['urcrnrlon']=opts.get('urcrnrlon',xlim[1])
      opts['urcrnrlat']=opts.get('urcrnrlat',ylim[1])
      opts['llcrnrlon']=opts.get('llcrnrlon',xlim[0])
      opts['llcrnrlat']=opts.get('llcrnrlat',ylim[0])

    else:
      opts['lon_0']=opts.get('lon_0',0.5*(xlim[0]+xlim[1]))
      opts['lat_0']=opts.get('lat_0',0.5*(ylim[0]+ylim[1]))


    # load proj from cache or create new:
    self.map_info_current=self.map_info_get()
    cacheLab=str(self.map_info_current)
    cch=cache.Cache()

    if cch.is_stored(cacheLab,'localfile'):
      if debug_lev==2: print(' -> loading proj from cache')
      m=cch.load(cacheLab,'localfile')
    else:
      if debug_lev==2: print(' -> creating proj')

##      print(opts)
##      print(self.config['proj.name'])
      m=Basemap(projection=self.config['proj.name'],
                resolution=self.config['proj.resolution'],
                **opts)

##      try:
##        m=Basemap(projection=self.config['proj.name'],
##                  resolution=self.config['proj.resolution'],
##                  **self.config['proj.options'])
##      except:
##
##        if self.config['proj.name'] in ('spstere','npstere'):
##          m = Basemap(projection=self.config['proj.name'],
##                    resolution=self.config['proj.resolution'],
##                    #lon_0=lon_0,boundinglat=blat)
##                    **opts)
##        elif self.config['proj.name'] == 'stere':
##          m = Basemap(projection=self.config['proj.name'],
##                    resolution=self.config['proj.resolution'],
##                    urcrnrlon=xlim[1], urcrnrlat=ylim[1],
##                    llcrnrlon=xlim[0], llcrnrlat=ylim[0],
##                    **opts)
##        else:
##          m = Basemap(projection=self.config['proj.name'], lat_ts=0.0,
##                    resolution=self.config['proj.resolution'],
##                    urcrnrlon=xlim[1], urcrnrlat=ylim[1],
##                    llcrnrlon=xlim[0], llcrnrlat=ylim[0],
##                    **opts)

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
    try:
      pl.sca(self.ax) # does not work inside tk gui with Figure
    except: pass

    # better do this before draw anything else, otherwise the drawmapboundary will be
    # called many times (by set_axes_limits, depending on the projection used)
    if self.config['proj.oceans_add']:
      self.map.drawmapboundary(ax=self.ax,**self.config['proj.oceans'])

    if self.config['proj.coast_add']:
      self.map.drawcoastlines(ax=self.ax,**self.config['proj.coast'])


    if self.config['proj.continents_add']:
      try: # may fail for high zooms
        self.map.fillcontinents(ax=self.ax,**self.config['proj.continents'])
      except: pass


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
      else: parallels=self.config['proj.parallels_vals']
      self.map.drawparallels(parallels,ax=self.ax,**self.config['proj.parallels'])

    for k in self.ax.spines:
      self.ax.spines[k].set_zorder(999.5)

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
      # accepted cmap: cmapname; list of colors, or one color for contours

      if cb.isstr(self.config['field.cmap']):
        try:
          # wont work for one color, like field.cmap='k'
          cmap=eval('pl.cm.'+self.config['field.cmap'])
        except:
          cmap=[self.config['field.cmap']]

      else: cmap=self.config['field.cmap']

      try: # list of colors for contour of contourf
        iter(cmap)
        args['colors']=cmap
      except:
        args['cmap']=cmap

    x,y=self._convCoord(x,y)
    meth=eval('self.ax.'+self.config['field.plot'])

#    if self.config['field.clim']:
#      vmin,vmax=self.config['field.clim']
#      args['vmin']=vmin
#      args['vmax']=vmax

    if self.config['field.clim'] is False:
       if self.config['field.cvals'] is False:
          self.config['field.clim']=v.min(),v.max()
       else:
         if isinstance(self.config['field.cvals'],int):
           self.config['field.clim']=v.min(),v.max()
         else:
           self.config['field.clim']=self.config['field.cvals'][0],self.config['field.cvals'][-1]

    args['vmin'],args['vmax']=self.config['field.clim']

#    if self.config['field.linecolors']:
#      args['colors']=self.config['field.linecolors']
    if self.config['field.linewidths']:
      args['linewidths']=self.config['field.linewidths']

    if self.config['field.extend']:
      args['extend']=self.config['field.extend']

    args['zorder']=self.config['plot.zorder']

    # about label: not supported by contour and contourf
    if self.config['field.plot'] in ('pcolor','pcolormesh'):
      plabel=self.config['plot.label'] if self.config['plot.label'] else self.label
      args['label']=plabel

    if self.config['field.plot'] in ('pcolor','pcolormesh'):
      self.handles['mappable']+=[meth(x,y,v,**args)]
    elif self.config['field.plot'] in ('contour','contourf'):
       # contour(f) does not like bad values in coordinates, even when masked
       # so:
       if np.ma.isMA(x): x[x.mask]=x.mean()
       if np.ma.isMA(y): y[y.mask]=y.mean()

#       if self.config['field.cvals'] is False: vals=20
#       else: vals=self.config['field.cvals']
       if self.config['field.cvals'] is False:
         self.config['field.cvals']=np.linspace(v.min(),v.max(),20)

       vals=self.config['field.cvals']

       self.handles['mappable']+=[meth(x,y,v,vals,**args)]


  def draw_vector_field(self,x,y,u,v):

    if x is None or y is None:
      ny,nx=v.shape
      x,y=np.meshgrid(np.arange(nx,dtype=v.dtype),np.arange(ny,dtype=v.dtype))


    # rotate u,v to projection angle:
    if hasattr(self,'map'):# and self.map:
###      print u.shape,v.shape,x.shape,y.shape
      u,v=self.map.rotate_vector(u,v,x,y)

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

    plabel=self.config['plot.label'] if self.config['plot.label'] else self.label
    kargs={}
    kargs['label']=plabel
    kargs['zorder']=self.config['plot.zorder']
    kargs.update(self.config['vfield.options'])
    self.handles['quiver']+=[self.ax.quiver(*args,**kargs)]
    ########self.handles['quiver']+=[self.ax.quiver(*args,**self.config['vfield.options'])]

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
      else: opts['orientation']='vertical'

      try:
        self.cb=pl.colorbar(self.handles['mappable'][-1],cax=self.cbax,**opts)
        if  self.cbbg: self.cbbg.set(visible=True)
        if  self.cbax: self.cbax.set(visible=True)
      except: pass # wont work if one contour



  def draw_1d(self,x=None,y=None):
    if x is None:
      x=np.arange(y.size,dtype=y.dtype)
#      args=self.config['d1_line.options']
#      self.handles['d1']+=self.ax.plot(y,**args)
#    else:

    plabel=self.config['plot.label'] if self.config['plot.label'] else self.label

#    if 1:
    plt=self.config['d1.plot']

    if  plt=='fill':
      args=self.config['d1_fill.options']
      args['label']=plabel
      args['zorder']=self.config['plot.zorder']
      x,y=self._convCoord(x,y)
      self.handles['d1']+=[self.ax.fill(x,y,**args)]

    elif  plt=='fill_between':
      args=self.config['d1_fill.options']
      args['label']=plabel
      args['zorder']=self.config['plot.zorder']

      y0=self.config['d1.y0']
      if y0=='min': Y0=self.ax.get_ylim()[0]
      elif y0=='max': Y0=self.ax.get_ylim()[1]
      else: Y0=y0

      x,y=self._convCoord(x,y) # should also convert y0 !! how!? no need for now...
      self.handles['d1']+=[self.ax.fill_between(x,y,Y0,**args)]

    elif plt=='plot':
      args=self.config['d1_line.options']
      args['label']=plabel
      args['zorder']=self.config['plot.zorder']
      x,y=self._convCoord(x,y)
      self.handles['d1']+=[self.ax.plot(x,y,**args)]

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
    self.label=''

    if len(args)==1: self.v=args[0]
    if len(args)==2: self.x,self.v=args
    elif len(args)==3: self.x,self.y,self.v=args
    elif len(args)>3:
      print('bad len(args) !!')


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
      try: print('%5s shape=%12s size=%d'%(i,str(a.shape),a.size))
      except: print('%5s %s'%(i,str(a)))

    if not self.msg: msg='EMPTY'
    else: msg=self.msg
    print(' == msg == %s'%msg)

  def plot(self,coords='auto',proj='auto',labels=True,extras=True,isExtra=0,**kargs):
    ax=kargs.pop('ax',None)
    cbax=kargs.pop('cbax',None)

    if ax:
      self.ax=ax
      axp=ax.get_position()
      self.config['axes.position']=axp.x0,axp.y0,axp.width,axp.height

    if cbax:
      self.cbax=cbax
      axp=cbax.get_position()
      self.config['colorbar.ax_position']=axp.x0,axp.y0,axp.width,axp.height

    parent=kargs.pop('inherit',0)
    if parent:
      self.inherit(parent)
      isExtra=1

    if isExtra: labels=False

    debug_lev=kargs.pop('debug_level',0)
    self.set_param(**kargs)

    # make vfield vars 2d to simplify the process:
    if isinstance(self.v,tuple) and self.v[0].ndim==1:
      self.x=self.x[:,np.newaxis]
      self.y=self.y[:,np.newaxis]
      self.v=self.v[0][:,np.newaxis],self.v[1][:,np.newaxis]

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

        if   xname!='y' and not self.y is None and (self.y.ndim==2 or self.y.size==vshape[0]):  y,yname=self.y,'y'
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


    try:
      isFig=pl.fignum_exists(self.fig.number)
      isAxFig=pl.fignum_exists(self.ax.figure.number)
    except:
      isFig=True
      isAxFig=True

    if hasattr(self,'fig') and isFig:
      if isExtra: pass
      else:
        if debug_lev==2: print(' -> will clear axes')
        self.clear_axes()
    else:
      if hasattr(self,'ax') and isAxFig:
        if debug_lev==2: print(' -> will use previous ax')
        if hasattr(self,'cbax'): self.init_figure(ax=self.ax,cbax=self.cbax)
        else: self.init_figure(ax=self.ax)
      else:
        if debug_lev==2: print(' -> will create new fig')
        self.init_figure()

    if proj and not isExtra:
      if not hasattr(self,'map') or not self.map or self.map_info_get()!=self.map_info_current:
        # make a new projection:
        if not self.config['axes.axis'] and not x is None:
          self.config['axes.axis']=x.min(),x.max(),y.min(),y.max()

        if debug_lev==2: print(' -> new projection needed')
        self.init_projection(debug_level=debug_lev)

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


    # vfield, 1d or 2d
    if np.iscomplexobj(self.v):
      if debug_lev==2: print(' -> will draw vfield')
      self.draw_vector_field(x,y,np.real(self.v),np.imag(self.v))

    # scalar:
    else:
      if ndim==2:
        if debug_lev==2: print(' -> will draw 2d')
        self.draw_scalar_field(x,y,self.v)
        if not isExtra:
          if debug_lev==2: print(' -> will draw colorbar')
          self.draw_colorbar()

      elif ndim==1:
        if debug_lev==2: print(' -> will draw 1d')
        self.draw_1d(x,y)

    # add labels:
    if labels:
      if ndim==2:
        tit=vlab
        if not self.t is None:
          #if isinstance(self.t,datetime.datetime):
          # now netcdftime.num2date returns array or netcdftime.datetime
          # obkjects instead of datetime.datetime
          import netcdftime
          if isinstance(self.t,datetime.datetime) or isinstance(self.t,netcdftime.datetime):
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
        if not xlab=='Unk': self.draw_label('xlabel',xlab)
        if yname=='v':
          if not vlab=='Unk': self.draw_label('title',vlab)
        else:
          if not ylab=='Unk': self.draw_label('ylabel',ylab)


    if proj and not isExtra:
      if debug_lev==2: print(' -> will draw projection', self.info['v'])
      self.draw_projection()


    if self.extra and extras:
      for e in self.extra:
        e.inherit(self)
        if debug_lev==2: print(' -> will plot extras')
        e.plot(isExtra=1,debug_level=debug_lev)

#      for e in self.extra:
#        e.fig=self.fig
#        e.ax=self.ax
#        if 0:
#          e.handles=self.handles # do I need this line!??
#        else: e.init_handles()
#        if hasattr(self,'map'):
#           e.map=self.map
#           e.map_info_copy(self)
#
#        if debug_lev==2: print ' -> will plot extras'
#        e.plot(labels=False,isExtra=1)

  def extra_find(self,lab=None):
    '''
    find extra obj by label in self.label
    '''

    if lab:
      labs=[i.label for i in self.extra]
      import fnmatch
      m=fnmatch.filter(labs,lab)
      return [self.extra[labs.index(i)] for i in m] # warning: labels must be different
    else:
      return {i.label:i for i in self.extra}


class MData:##################(Data):
  def __init__(self,obs,warnings=[]):
    self._data=obs

    self.errors=cb.unique([i.msg for i in obs])
    self.warnings=cb.unique(warnings)

    # remove empty:
    self.errors=[i for i in self.errors if len(i)]
    self.warnings=[i for i in self.warnings if len(i)]

#    self.msg=[]
#    for i in obs:
#      self.msg+=i.msg.split('\n')
#
#    # remoe repeated messages:
#    self.msg=cb.unique(self.msg)

    '''
    self.__dict__.update(obs[0].__dict__)

    # same config, also for extras, which should have the same len (or less!)
    for i in obs[1:]:
      i.config=obs[0].config
# not for extras ... test !
#      for j,e in enumerate(i.extra):
#        e.config=obs[0].extra[j].config

    #self.extra+=obs[1:] # shoud not be used cos change is done in plce
    # and self[0].extra will also be affected! thus:
    self.extra=self.extra+obs[1:]
    '''

  def plot(self,**kargs):
###    for i in self:
###      print i.config['plot.zorder']

    self[0].plot(**kargs)

    for i in self[1:]:
      i.plot(inherit=self[0],**kargs)

    return

    kargs['isExtra']=1
#    kargs['debug_level']=2
    for i in self[1:]:
      i.inherit(self[0])
      i.plot(**kargs)

  def find(self,lab=None):
    O=[]
    for i in self:
      O+=[i]+i.extra

    if not lab: return O

    import fnmatch
    labs=[i.label for i in O]
    if len(cb.unique(labs))==len(labs):
      m=fnmatch.filter(labs,lab)
      return [O[labs.index(i)] for i in m]
    else: # allowing repeated labels:
      out=[]
      for i in range(len(O)):
        if fnmatch.filter([labs[i]],lab): out+=[O[i]]

      return out


  def __getitem__(self,i):
    return self._data[i]

  def __len__(self): return len(self._data)

