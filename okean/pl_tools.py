import pylab as plt
import numpy as np


class InteractiveLine:
  '''Create and edit intercative broken line, polygon or spline

     - auxiliary points are created with left mouse button and stored in self.x and self.y
     - stop creation with other mouse button
     - spline or broken line is stored in self.xx and self.yy
     - press key "e" with figure selected to start edition:
       . click (left mouse) on aux points to move them
       . use right mouse to add a new point - click over a segment between 2 aux points
       . use 3rd mouse button to remove one point

    To use a predefined line, use the arg xy=(x,y)
    A command can be called after the line creation (arg cmd)
    To keep xlim, ylim, use arg axis = xlim+ylim

    mma, Texas A&M
  '''

  def __init__(self,ax=False,xy=False,cmd=False,cmd_disable_edit=True,axis=False,type='spline',nmax=-1):
    '''Line types:
      - broken line (type=blin)
      - closed broken line, ie, polygon (type=cblin)
      - spline (type=spline, default)
    '''

    if not ax: ax=plt.gca()

    self.x=np.array(())
    self.y=np.array(())
    self.line=False
    self.lineAux=False


    self.type=type
    self.nmax=nmax
    self.cmd=cmd
    self.axis=axis
    self.cmd_disable_edit=cmd_disable_edit

    self.ax=ax
    self.figure=ax.figure
    self.events={}

    if xy:
      self.x,self.y=xy
      self.update_line()

    self.start()

  def start(self):
    # start line creation
    self.events['mousedown'] = self.figure.canvas.mpl_connect('button_press_event', self.mouseclick)
    self.events['mousemove'] = self.figure.canvas.mpl_connect('motion_notify_event', self.mousemove)

  def edit(self):

    def onpick(event):
      thisline = event.artist

      # picking line aux:
      if thisline is self.lineAux:
        if  event.mouseevent.button==1: # move
          self.pind=event.ind
          self.events['mousemove'] = self.figure.canvas.mpl_connect('motion_notify_event', self.mousemove_edit)
          self.events['mouseup'] = self.figure.canvas.mpl_connect('button_release_event',self.stop_edit)
        elif  event.mouseevent.button==2: # remove
          if self.type=='cblin':
            rmMax=3
            event.ind=event.ind[0]
            if event.ind in (0,len(self.x)-1):
              self.x=np.delete(self.x,0)
              self.y=np.delete(self.y,0)
              self.x=np.delete(self.x,-1)
              self.y=np.delete(self.y,-1)
              self.x=np.append(self.x,self.x[0])
              self.y=np.append(self.y,self.y[0])
            else:
              self.x=np.delete(self.x,event.ind)
              self.y=np.delete(self.y,event.ind)

          else:
            rmMax=2
            self.x=np.delete(self.x,event.ind)
            self.y=np.delete(self.y,event.ind)

          if self.x.size<rmMax: # remove all line if only one point
            self.stop()
            self.remove_line()
          else:
            self.update_line()

      # picking line:
      elif thisline is self.line:
        if  event.mouseevent.button==3: # add new point
          x=self.line.get_xdata()
          y=self.line.get_ydata()
          d=(x-event.mouseevent.xdata)**2+(y-event.mouseevent.ydata)**2
          i=np.where(d==d.min())[0][0]
          if self.Iaux[i-2]==self.Iaux[i+2]==self.Iaux[i]:
            i=self.Iaux[i]
          else: return # too close from spline aux point

          self.x=np.concatenate((self.x[:i+1],[event.mouseevent.xdata],self.x[i+1:]))
          self.y=np.concatenate((self.y[:i+1],[event.mouseevent.ydata],self.y[i+1:]))
          self.update_line()



    if self.line:
      self.lineAux.set(mfc='#cccccc',mec='b',ms=5)
      self.events['pick']=self.figure.canvas.mpl_connect('pick_event', onpick)
      self.figure.canvas.draw()


  def edit_stop(self):
    if self.line:
      self.lineAux.set(mfc='w',mec='r',ms=4)
      self.figure.canvas.mpl_disconnect(self.events['pick'])
      self.figure.canvas.draw()

  def onkey(self,event):
    if event.key=='e':
      if self.lineAux.get_markersize()==4: self.edit()
      else: self.edit_stop()


  def update_line(self,calc_inds=True):
    # update line while creating of editing...
    if not self.line:
      self.lineAux,=self.ax.plot(self.x,self.y,'s',mfc='w',mec='r',picker=5,zorder=1001,ms=4)
      self.line,=self.ax.plot(self.x,self.y,'r',marker='None',lw=0.5,picker=2,zorder=1000,alpha=.7)
    else:

      n=np.unique((self.x+self.y*1j)).size
      if self.type=='spline' and n>3:
        from scipy import interpolate
        if n<self.x.size:
          tck,u = interpolate.splprep([self.x[:-1],self.y[:-1]],s=0)
        else:
          tck,u = interpolate.splprep([self.x,self.y],s=0)

        xnew= np.arange(0,1.01,0.01)
        xnew=np.linspace(0,1,10*(n-1))
        out = interpolate.splev(xnew,tck,der=0)
        xnew,ynew=out[0],out[1]

      else: # broken line
        xnew=np.array((),'f')
        ynew=np.array((),'f')
        for i in range(self.x.size-1):
          if self.x[i+1]!=self.x[i]:
            m=(self.y[i+1]-self.y[i])/(self.x[i+1]-self.x[i])
            xx=np.linspace(self.x[i],self.x[i+1],10)
            yy=m*(xx-self.x[i])+self.y[i]
          else:
            yy=np.linspace(self.y[i],self.y[i+1],10)
            xx=yy*0+self.x[i]

          xnew=np.append(xnew,xx[:-1])
          ynew=np.append(ynew,yy[:-1])

        # add last:
        xnew=np.append(xnew,xx[-1])
        ynew=np.append(ynew,yy[-1])


      if calc_inds:
        # find segment number for each point:
        iaux=[]
        for i  in range(self.x.size):
              d=(xnew-self.x[i])**2+(ynew-self.y[i])**2
              j=np.where(d==d.min())[0][0]
              iaux+=[j]

        # set segment number:
        self.Iaux=np.zeros(xnew.shape,'i')
        i=-1
        for i in range(len(iaux)-2):
              self.Iaux[iaux[i]:iaux[i+1]]=i

        self.Iaux[iaux[i+1]:]=i+1

      # update lines:
      self.line.set_data(xnew,ynew)
      self.lineAux.set_data(self.x,self.y)

    # aux points stored in self.x, self.y
    # also store in xx and yy the line points if spline, and a copy of x,y if broken line
    if self.type=='spline':
      self.xx=self.line.get_xdata()
      self.yy=self.line.get_ydata()
    else: self.xx,self.yy=self.x,self.y

    if self.axis: self.ax.axis(self.axis)
    self.figure.canvas.draw()


  def remove_line(self):
    # delete all
    self.ax.lines.remove(self.line)
    self.ax.lines.remove(self.lineAux)
    self.figure.canvas.draw()
    self.line=False
    self.lineAux=False
    self.x=np.array(())
    self.y=np.array(())

  def mouseclick(self,event):
    # line creation...
    if event.button==1 and event.inaxes: # build line:
      if len(self.x):
        self.x[-1]=event.xdata
        self.y[-1]=event.ydata
        self.x=np.append(self.x,event.xdata)
        self.y=np.append(self.y,event.ydata)
      else:
        self.x=np.append(self.x,[event.xdata]*2)
        self.y=np.append(self.y,[event.ydata]*2)

      if self.nmax>0 and self.x.size>self.nmax:
        self.x=self.x[:-1]
        self.y=self.y[:-1]
        self.update_line()
        self.stop()


    else: # stop without using last point:
      if self.x.size:
        if self.x.size==2: # just one point--> remove line
          self.remove_line()
        else:
          self.x=self.x[:-1]
          self.y=self.y[:-1]

          if self.type=='cblin':
            self.x=np.append(self.x,self.x[0])
            self.y=np.append(self.y,self.y[0])

          self.update_line()

      self.stop()


  def mousemove(self,event):
    # show line change while creating ...
    if self.x.size and event.inaxes:
      self.x[-1]=event.xdata
      self.y[-1]=event.ydata
      self.update_line(calc_inds=False)


  def mousemove_edit(self,event):
      # edit line
      if self.x.size and event.inaxes:
        self.x[self.pind]=event.xdata
        self.y[self.pind]=event.ydata
        self.update_line(calc_inds=False)


  def stop(self,run=True):
    # stop line creation
    self.figure.canvas.mpl_disconnect(self.events['mousedown'])
    self.figure.canvas.mpl_disconnect(self.events['mousemove'])

    exec_cmd=run and callable(self.cmd)
    enable_edit=False
    if not exec_cmd or (exec_cmd and not self.cmd_disable_edit): enable_edit=True

    if enable_edit:
        self.events['key']=self.figure.canvas.mpl_connect('key_press_event', self.onkey)

    # do something after line is created:
    if exec_cmd: return self.cmd()


  def stop_edit(self,event):
    # stop line editing
    self.figure.canvas.mpl_disconnect(self.events['mousemove'])


class InteractiveRect:
  def __init__(self,ax=False,cmd=False,axis=False):
    if not ax: ax=plt.gca()
    self.x=[]
    self.y=[]
    self.line=False

    self.cmd=cmd
    self.axis=axis

    self.ax=ax
    self.figure=ax.figure
    self.events={}

    self.pind=2 # bott-right corner
    self.start()

  def start(self):
    # start line creation
    self.events['mousedown'] = self.figure.canvas.mpl_connect('button_press_event', self.mouseclick)
    self.events['mouseup']   = self.figure.canvas.mpl_connect('button_release_event', self.stop)

    # enable/disable edition witj key "e"
    self.events['key']=self.figure.canvas.mpl_connect('key_press_event', self.onkey)


  def edit(self):

    def onpick(event):
      thisline = event.artist
      if thisline is self.line:
        if  event.mouseevent.button==1: # move
          self.pind=event.ind
          self.events['mousemove'] = self.figure.canvas.mpl_connect('motion_notify_event', self.mousemove)
        else: # remove
          self.remove_line()

    if self.line:
      self.line.set(mfc='#cccccc',mec='b',ms=5)
      self.events['pick']=self.figure.canvas.mpl_connect('pick_event', onpick)
      self.figure.canvas.draw()


  def edit_stop(self):
    if self.line:
      self.line.set(mfc='w',mec='r',ms=0)
      self.figure.canvas.mpl_disconnect(self.events['pick'])
      self.figure.canvas.draw()

  def onkey(self,event):
    print 'key !!'
    if event.key=='e':
      if self.line.get_markersize()==0: self.edit()
      else: self.edit_stop()

  def update_line(self):
    if not self.line:
      self.line,=self.ax.plot(self.x,self.y,'r-',picker=5,marker='s',ms=0,mfc='w',mec='r',alpha=.7)
    else: self.line.set_data(self.x,self.y)
    if self.axis: self.ax.axis(self.axis)
    self.figure.canvas.draw()


  def remove_line(self):
    self.ax.lines.remove(self.line)
    self.figure.canvas.draw()
    self.line=False
    self.x=[]
    self.y=[]


  def mouseclick(self,event):
    self.button=event.button
    self.x=np.zeros(9,'f')+event.xdata
    self.y=np.zeros(9,'f')+event.ydata

    self.events['mousemove'] = self.figure.canvas.mpl_connect('motion_notify_event', self.mousemove)

  def mousemove(self,event):
    if event.inaxes:## self.x and not event.xdata is None and not event.ydata is None:
      # pind 2
      ix  = [2,3,4]
      ix2 = [1,5]
      iy  = [0,1,2,8]
      iy2 = [3,7]
      if np.any(self.pind==3):
        iy=[]
        iy2=[]
      elif np.any(self.pind==4):
        iy=[4,5,6]
      elif np.any(self.pind==5):
        iy=[4,5,6]
        ix=[]
        ix2=[]
      elif np.any(self.pind==6):
        iy=[4,5,6]
        ix=[0,6,7,8]
      elif np.any(self.pind==7):
        ix=[0,6,7,8]
        iy=[]
        iy2=[]
      elif np.any(self.pind==0):
        ix=[0,6,7,8]
      elif np.any(self.pind==1):
        ix=[]
        ix2=[]

      self.x[ix]=event.xdata
      self.y[iy]=event.ydata
      self.x[ix2]=(self.x.max()+self.x.min())/2.
      self.y[iy2]=(self.y.max()+self.y.min())/2.


      self.update_line()

  def stop(self,e):
    events='mousedown','mousemove'
    for n in events:
      if self.events.has_key(n): self.figure.canvas.mpl_disconnect(self.events[n])

    self.x=np.array(self.x)
    self.y=np.array(self.y)

#    try:
#     if not e is False: return self.cmd()
#    except: pass

    # do something after line is created:
    if callable(self.cmd): return self.cmd()


def _trend_cmap(cmap,x):
  def original_colors(ncolors):
    mm=plt.cm.ScalarMappable(cmap=cmap)
    mm.set_clim((0,1))
    var=np.linspace(0,1,ncolors,1)
    return mm.to_rgba(var)

  Ncolors=len(x)

  Ncolors0=Ncolors*10
  x0=np.arange(0,1,1./Ncolors0)
  x0=x0/x0[-1]

  x=x/x[-1]
  colors0=original_colors(Ncolors0)
  colors=np.zeros((Ncolors,3))
  for i in range(3): colors[:,i]=np.interp(x,x0,colors0[:,i])
  return plt.cm.colors.ListedColormap(colors)


def trend_cmap(cmap,N=2):
  '''
  Ex:
  import numpy as np
  import pylab as pl

  x=np.linspace(-3,3,100)
  x,y=np.meshgrid(x,x)
  z = 3*(1-x)**2*np.exp(-(x**2)-(y+1)**2)\
      -10*(x/5-x**3-y**5)*np.exp(-x**2-y**2)\
      - 1/3*np.exp(-(x+1)**2-y**2)

  cmap0=pl.cm.jet
  cmap1=trend_cmap(cmap0,2)
  cmap2=trend_cmap(cmap0,4)

  pl.figure()
  ax0=pl.subplot(131,xticks=[],yticks=[])
  k=ax0.pcolormesh(z,cmap=cmap0)
  pl.colorbar(k)
  ax1=pl.subplot(132,xticks=[],yticks=[])
  k=ax1.pcolormesh(z,cmap=cmap1)
  pl.colorbar(k)

  ax2=pl.subplot(133,xticks=[],yticks=[])
  k=ax2.pcolormesh(z,cmap=cmap2)
  pl.colorbar(k)


  '''

  if isinstance(N,int):
    Ncolors=256
    x=np.arange(0,1,1./Ncolors)**N
  else:
    #[(a,A),(b,B),C]
    if len(N)==3: #[(a,A),(b,B),C]
      a,A=N[0]
      b,B=N[1]
      C=N[2]
      a,A,b,B=[i*255./C for i in (a,A,b,B)]
      C=256
      x=np.zeros(C,'f')

      x1=np.arange(0,a,a/A)
      x2=np.arange(a,b,(b-a)/(B-A))
      x3=np.arange(b,C,(C-b)/(C-B))
      x=np.concatenate((x1,x2,x3))
    elif len(N)==2: # [(a,B),C]
      a,A=N[0]
      C=N[1]
      a,A=[i*255./C for i in (a,A)]
      C=256
      x=np.zeros(C,'f')

      x1=np.arange(0,a,a/A)
      x2=np.arange(a,C,(C-a)/(C-A))
      x=np.concatenate((x1,x2))


  return _trend_cmap(cmap,x)


class cmap_aux:
  def __init__(self): pass

  def _store(self,names,cmaps):
    self.cmap_d={}

    for i in range(len(names)):
       self.cmap_d[names[i]]=cmaps[i]
       setattr(self,names[i],cmaps[i])
       setattr(self,names[i]+'_r',self._invert(cmaps[i]))

    for k in names: self.cmap_d[k]=getattr(self,k)

  def _invert(self,c):
    if isinstance(c,plt.matplotlib.colors.LinearSegmentedColormap):
      return invert_lsc(c) 
    elif isinstance(c,plt.matplotlib.colors.ListedColormap):
      return invert_lc(c)


  def show(self):
    '''Display colormaps'''

    a=np.outer(np.arange(0,1,0.01),np.ones(10))
    plt.figure(figsize=(8,5))
    plt.subplots_adjust(top=0.8,bottom=0.05,left=0.05,right=0.95)

    l=len(self.cmap_d)+1
    names=sorted(self.cmap_d.keys())
    for i,name in enumerate(names):
      m=self.cmap_d[name]
      plt.subplot(1,l,i+1)
      plt.axis("off")
      plt.imshow(a,aspect='auto',cmap=m,origin="lower")
      plt.text(0.5,1.01,name,rotation=50,fontsize=10,
              transform=plt.gca().transAxes,va='bottom')

    plt.show()



class ucmaps(cmap_aux):
  def __init__(self):

    # names and cmaps:
    names='mod_jet','mod_jet2','mod_jet3','freshwater','freshwater2',\
          'oxygen'

    cmaps=[getattr(self,'gen_'+n)() for n in names]

    for i in range(16):
      names+='oceano_%02d'%i,
      cmaps+=[self.gen_oceano(i)]

    # store cmap and inverted cmap:
    self._store(names,cmaps)


  def gen_mod_jet(self):
    cdict =   {'red':   ((0.00, 0.4,  0.4),
                         (0.35, 0.3,  0.3),
                         (0.66, 1.0,  1.0),
                         (0.85, 0.9,  0.9),
                         (0.93, 0.75,  0.75),
                         (1.00, 0.83, 0.83)),
               'green': ((0.00,  0.4, 0.4),
                         (0.125, 0.3, 0.3),
                         (0.375, 1.0, 1.0),
                         (0.64,  1.0, 1.0),
                         (0.75,  0.5, 0.5),
                         (0.93,  0.5, 0.5),
                         (1.00,  0.8, 0.8)),
               'blue':  ((0.00, 0.7, 0.7),
                         (0.11, 1.0, 1.0),
                         (0.34, 1.0, 1.0),
                         (0.65, 0.0, 0.0),
                         (0.85,  0.6, 0.6),
                         (1.00, 0.8, 0.8))}

    return plt.matplotlib.colors.LinearSegmentedColormap('mod_jet',cdict,256)

  def gen_mod_jet2(self):
    cdict =   {'red':   ((0.00, 0.4,  0.4),
                         (0.35, 0.3,  0.3),
                         (0.66, 1.0,  1.0),
                         (0.85, 0.9,  0.9),
                         (0.93, 0.75,  0.75),
                         (0.98, 0.83,  0.83),
                         (0.99, 0.7,  0.7),
                         (1.00, 0.83, 0.83)),
               'green': ((0.00,  0.4, 0.4),
                         (0.125, 0.3, 0.3),
                         (0.375, 1.0, 1.0),
                         (0.64,  1.0, 1.0),
                         (0.75,  0.5, 0.5),
                         (0.93,  0.5, 0.5),
                         (1.00,  0.8, 0.8)),
               'blue':  ((0.00, 0.7, 0.7),
                         (0.11, 1.0, 1.0),
                         (0.34, 1.0, 1.0),
                         (0.65, 0.0, 0.0),
                         (0.85,  0.6, 0.6),
                         (1.00, 0.8, 0.8))}

    return plt.matplotlib.colors.LinearSegmentedColormap('mod_jet2',cdict,256)

  def gen_mod_jet3(self):
    cdict =   {'red':   ((0.00, 0.4,  0.4),
                         (0.35, 0.3,  0.3),
                         (0.66, 1.0,  1.0),
                         (0.85, 0.9,  0.9),
                         (0.93, 0.75,  0.75),
                         (1.00, 0.92, 0.83)),
               'green': ((0.00,  0.4, 0.4),
                         (0.125, 0.3, 0.3),
                         (0.375, 1.0, 1.0),
                         (0.64,  1.0, 1.0),
                         (0.75,  0.5, 0.5),
                         (0.93,  0.5, 0.5),
                         (1.00,  0.88, 0.8)),
               'blue':  ((0.00, 0.7, 0.7),
                         (0.11, 1.0, 1.0),
                         (0.34, 1.0, 1.0),
                         (0.65, 0.0, 0.0),
                         (0.85,  0.6, 0.6),
                         (1.00, 0.88, 0.8))}

    return plt.matplotlib.colors.LinearSegmentedColormap('mod_jet',cdict,256)

  def gen_freshwater(self):
    cdict = {'red': ((0.0,  0.0, 0.95),
                    (1/3., 0.65, 1.0),
                    (2/3., 0.9, 0.75),
                    (1.0,  0.3, 1.0)),

           'green': ((0.0,  0.0, 0.95),
                     (1/3., 0.65, 0.85),
                     (2/3., 0.4, 0.1),
                     (1.0,  0.0, 1.0)),

           'blue': ((0.0,  0.0, 0.95),
                    (1/3., 0.65, 1.0),
                    (2/3., 0.9, 0.75),
                    (1.0,  0.3, 1.0))}

    Ccdict = {'red': ((0.0,  0.0, 0.95),
                    (1/3., 0.65, 0.725),
                    (2/3., 0.502, 0.263),
                    (1.0,  0.078, 1.0)),

           'green': ((0.0,  0.0, 0.95),
                     (1/3., 0.65, 0.788),
                     (2/3., 0.545, 0.396),
                     (1.0,  0.114, 1.0)),

           'blue': ((0.0,  0.0, 0.95),
                    (1/3., 0.65, 1.000),
                    (2/3., 0.690, 0.859),
                    (1.0,  0.251, 1.0))}

    return plt.matplotlib.colors.LinearSegmentedColormap('freshwater',cdict,256)


  def gen_freshwater2(self):
    cdict = {'red': ((0.0,  0.0, 0.835),
                    (1/3., 0.610, 1.0),
                    (2/3., 0.9, 0.75),
                    (1.0,  0.3, 1.0)),

           'green': ((0.0,  0.0, 0.784),
                     (1/3., 0.404, 0.85),
                     (2/3., 0.4, 0.1),
                     (1.0,  0.0, 1.0)),

           'blue': ((0.0,  0.0, 0.902),
                    (1/3., 0.902, 1.0),
                    (2/3., 0.9, 0.75),
                    (1.0,  0.3, 1.0))}

    return plt.matplotlib.colors.LinearSegmentedColormap('freshwater',cdict,256)


  def gen_oxygen(self,v=(0,135.,300.)):
    a=v[1]/v[2]
    cdict = {'red': ((0.0,  0.0,  1.0),
                     (a,    1.0,  0.93),
                     (1.0,  0.53, 1.0)),

           'green': ((0.0,  0.0,  0.0),
                     (a,    1.0,  0.93),
                     (1.0,  0.53, 1.0)),

           'blue': ((0.0,  0.0, 0.0),
                    (a,    0.0, 0.93),
                    (1.0,  0.53, 1.0))}

    return plt.matplotlib.colors.LinearSegmentedColormap('oxygen',cdict,256)


  def gen_oceano(self,ind):
    '''ind=0..16'''
    if ind==0:
      cdict = {'red':   ((0.00,  0.0,  0.95),
                         (0.4,   0.6,  0.6),
                         (0.7,   0.0,  0.0),
                         (1.0,   0.65, 0.0)),

               'green': ((0.00, 0.0,   0.95),
                         (0.4,  0.6,   0.6),
                         (0.7,  0.0,   0.0),
                         (1.0,  0.0,   0.0)),

               'blue':  ((0.00, 0.0,   1.0),
                         (0.3,  0.61,  0.61),
                         (0.9,  0.0,   0.0),
                         (1.0,  0.85,  0.0))}
      cdict=plt.cm.revcmap(cdict)

    else:
      i=((28,23,3),(28,28,3),(28,3,23),(28,3,28),
         (28,3,3),(3,23,23),(3,23,28),(3,23,3),
         (3,28,23),(3,28,28),(3,28,3),(3,3,23),
         (3,3,28),(3,3,3),(23,3,28))

      r,g,b=i[ind-1]
      _cm=plt.matplotlib._cm
      cdict = {
        'red': _cm.gfunc[r],
        'green': _cm.gfunc[g],
        'blue': _cm.gfunc[b]}

    cdict=plt.cm.revcmap(cdict)
    return plt.matplotlib.colors.LinearSegmentedColormap('oceano_%02d'%ind,cdict,256)


def invert_lsc(C):
  '''Invert linearSegmentedColormap'''
  r=plt.cm.revcmap(C._segmentdata)
  return plt.matplotlib.colors.LinearSegmentedColormap(C.name+'_r',r,C.N)


def invert_lc(C):
  '''Invert ListedColormap'''
  r=np.flipud(C.colors)
  return plt.matplotlib.colors.ListedColormap(r, name=C.name+'_r', N=C.N)


class ncview_cmaps(cmap_aux):
  def __init__(self,filespath='auto'):
    self.filespath=filespath

    # find ncview colormap files:
    if self.filespath=='auto':
      import os
      here=os.path.dirname(os.path.abspath(__file__))
      self.filespath=os.path.join(here,'data','ncview_cmaps')

    import glob
    files=glob.glob(os.path.join(self.filespath,'colormaps_*.h'))

    # names and cmaps:
    names=[]
    cmaps=[]
    for f in files:
      cmaps+=[self.gen_cmap(f)]
      names+=[cmaps[-1].name]

    # store cmap and inverted cmap:
    self._store(names,cmaps)


  def read_file(self,f):
    '''Read ncview colormaps_<name>.h file'''

    l=open(f).readlines()

    i=-1
    for k in l:
      i+=1
      if k.startswith('static'): break

    l=l[i:]

    i0=l[0].find('{')
    i1=l[-1].find('}')

    l[0]=l[0][i0+1:]
    l[-1]=l[-1][:i1]

    r=[]
    for i in range(len(l)):
      line=l[i].replace(',',' ').strip()
      vals=[int(j) for j in line.split()]
      r=np.hstack((r,vals))

    r.shape=r.size/3,3
    return r/255.


  def gen_cmap(self,file,name='auto',N=None):
    '''Read ncview colormaps_<name>.h file'''

    if name=='auto': name=file.split('colormaps_')[-1][:-2]
    r=self.read_file(file)
    return plt.matplotlib.colors.ListedColormap(r, name=name, N=N)


# my cmaps:
cm=ucmaps()

# other colormaps,
# basemap:
from mpl_toolkits. basemap import cm as cm_basemap

# ncview:
cm_ncview=ncview_cmaps()

