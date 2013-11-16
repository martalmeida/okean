import pylab as pl
import numpy as np
import ticks


def fillout(x,y,lims,color=False,**kargs):
  # check polygon direction... ie, check need for flip x,y !!
  xi,xe,yi,ye=lims

  x=np.flipud(x)
  y=np.flipud(y)

  i=np.where(x==x.min())[0][0]

  x=np.concatenate([x[i:], x[:i+1]])
  y=np.concatenate([y[i:], y[:i+1]])

  tmpx=np.array([xi,   xi, xe, xe, xi, xi,   x[0]])
  tmpy=np.array([y[0], ye, ye, yi, yi, y[0], y[0]])

  x=np.concatenate([tmpx, x])
  y=np.concatenate([tmpy, y])

  if color is False: return x,y
  else:
    o=pl.fill(x,y,color,**kargs)
    return x,y,o


def wind_rose(D,F,**kargs):
  '''
  Example:
  import numpy as np

  d=np.arange(0,360,10)
  D=np.array(())
  V=np.array(())
  for i in range(len(d)):
    n=d[i]/10.
    D=np.append(D,np.ones(n)*d[i])
    V=np.append(V,np.arange(n))

  wind_rose(D,V)

  '''
  if isinstance(D,list) or isinstance(D,tuple): D=np.array(D)
  if isinstance(F,list) or isinstance(F,tuple): F=np.array(F)

  if D is 0:
    d=np.arange(0,360,10)
    D=np.array([])
    F=np.array([])
    for i in range(len(d)):
      n=d[i]/10.
      D=np.append(D,d[i]+np.zeros((1,n)))
      F=np.append(F,range(1,n+1))


  dtype='standard'
  nAngles=36
  ri=1/30.
  quad=0
  legType=0
  percBg='w'
  titStr=''
  legStr=''
  cmap='jet'
  colors=[]
  Ag=False # intensity subdivs.
  ci=[] # percentage circles
  lineColors='k'
  borderColor=False
  onAxes=False
  iflip=0
  inorm=0
  parent=0
  IncHiLow=1 # include values higher and lower that the limits of Ag.
  figPos=False
  axPos=False
  FontSize=8
  LineWidth=.5
  legTickLab=True
  NSEWlab=True
  percLabels=True
  labels=True

  for k in kargs.keys():
    if   k=='dtype':     dtype       = kargs[k]
    elif k=='n':         nAngles     = kargs[k]
    elif k=='ri':        ri          = kargs[k]
    elif k=='quad':      quad        = kargs[k]
    elif k=='legtype':   legType     = kargs[k]
    elif k=='percbg':    percBg      = kargs[k]
    elif k=='labtitle':  titStr      = kargs[k]
    elif k=='lablegend': legStr      = kargs[k]
    elif k=='cmap':      cmap        = kargs[k]
    elif k=='colors':    colors      = kargs[k]
    elif k=='di':        Ag          = kargs[k]
    elif k=='ci':        ci          = kargs[k]
    elif k=='lcolor':    lineColors  = kargs[k]
    elif k=='bcolor':    borderColor = kargs[k]
    elif k=='iflip':     iflip       = kargs[k]
    elif k=='inorm':     inorm       = kargs[k]
    elif k=='parent':    parent      = kargs[k]
    elif k=='incout':    IncHiLow    = kargs[k]
    elif k=='figpos':    figPos      = kargs[k]
    elif k=='axpos':     axPos       = kargs[k]
    elif k=='fontsize':  FontSize    = kargs[k]
    elif k=='linewidth': LineWidth   = kargs[k]
    elif k=='lablegticks':legTickLab  = kargs[k]
    elif k=='NSEWlab':   NSEWlab  = kargs[k]
    elif k=='labels':    labels  = kargs[k]
    elif k=='percLabels': percLabels  = kargs[k]
    elif k=='ax':
      try:
        onAxes,onAxesX,onAxesY,onAxesR= kargs[k]
      except:
        print ':: cannot place wind rose on axes, bad argument for ax'
        return

  if not labels:
    legTickLab=NSEWlab=percLabels=False

  #  open fig:
  if parent:
    fig=parent
  else:
    if figPos: fig=pl.figure(figsize=figPos)
    else:      fig=pl.figure()

  # other options:
  # size of the full rectangle:
  rs=1.2
  rl=1.7
  # axes W/H = (rs+rl)/(2*rs)

  # directions conversion:
  if dtype=='meteo': D=np.mod(-90-D,360)


  # angles subdivisons:
  D=np.mod(D,360.)
  Ay=np.linspace(0,360,nAngles+1)-0.5*360/nAngles

  # calc instensity subdivisions:
  if Ag is False:
    Ag=ticks.loose_label_n(F.min(),F.max(),ntick=7)

  E=np.zeros([len(Ay)-1,len(Ag)-1])
  for i in range(len(Ay)-1):
    if i==0:
      I=( (D>=Ay[i]) & (D<Ay[i+1]) ) | (D>=Ay[-1])
    else:
      I=(D>=Ay[i]) & (D<Ay[i+1])

    b=F[I]

    for j in range(len(Ag)-1):
      if j==len(Ag)-2:
        J=((b>=Ag[j]) & (b<=Ag[j+1])) # include data with last Ag
      else:
        J=((b>=Ag[j]) & (b<Ag[j+1]))
      E[i,j]=len(np.where(J)[0])

    if IncHiLow:
      E[i,0]=len(np.where(b<Ag[1])[0])
      E[i,-1]=len(np.where(b>=Ag[-2])[0])

  b=np.sum(E,1)/len(D)*100

  # normalize data: TODO

  # check if has values higher or lower than the Ag limits
  hasH=len(np.where(F>Ag[-1])[0])
  hasL=len(np.where(F<Ag[0])[0])

  # calc number of percentage circles to draw:
  if not ci:
    if inorm:
      ci=[25, 50, 75]
      g=120
      ncircles=3
    else:
      dcircles=np.array([1, 2, 5, 10, 15, 20, 25, 30, 50])
      ncircles=3
      d=abs(1./(dcircles/max(b))-ncircles)

      i=np.where(d==np.min(d))[0]
      d=dcircles[i[0]]
      if d*ncircles<max(b): ncircles=ncircles+1
      ci=np.arange(1,ncircles+1)*d
      g=ncircles*d
  else:
    ncircles=len(ci)
    g=max(max(ci),max(b))

  # plot axes, percentage circles and percent. data:
  if onAxes:
    wrAx=onAxes
  else:
    if axPos: wrAx=fig.add_axes(axPos)
    else:     wrAx=fig.add_axes()

  try: pl.axes(wrAx)
  except: pass

  ri=g*ri
  pl.fill([-rs*g, rl*g, rl*g, -rs*g],[-rs*g, -rs*g, rs*g, rs*g],'w',edgecolor='w',linewidth=LineWidth,label='wr_bg')
  pl.plot([-g-ri, -ri, np.nan, ri, g+ri, np.nan, 0, 0, np.nan, 0, 0],
          [0, 0, np.nan, 0, 0, np.nan, -g-ri, -ri, np.nan, ri, g+ri],':',color=lineColors,linewidth=LineWidth)

  t0=np.arange(0,361)*np.pi/180
  labs=[]
  Ang=np.array([1, 3, 5, 7])/4.*np.pi

  Valign='top', 'top', 'bottom', 'bottom'
  Halign='center', 'left', 'left', 'right'

  for i in range(ncircles):
    x=(ci[i]+ri)*np.cos(t0)
    y=(ci[i]+ri)*np.sin(t0)

    pl.plot(x,y,':',color=lineColors,linewidth=LineWidth);

    if percLabels: pl.text((ci[i]+ri)*np.cos(Ang[quad]),(ci[i]+ri)*np.sin(Ang[quad]),str(ci[i])+'%',
      verticalalignment=Valign[quad],horizontalalignment=Halign[quad],
      backgroundcolor=percBg,fontsize=FontSize)


  # calc colors:
  if not colors:
    cor=range(len(Ag)-1)
    q=pl.cm.ScalarMappable()
    q.set_clim([Ag[0], Ag[-2]])
    for j in range(len(Ag)-1):
      cor[j]=q.to_rgba(Ag[j])
      #cor[j]=caxcolor(Ag(j),[Ag[0], Ag[-2]],cmap);
  else:
    cor=colors

  # fill data:
  n=np.sum(E,1)
  #if  iflip...
  for i in range(len(Ay)-1):
    if n[i]:
      t=np.linspace(Ay[i],Ay[i+1],20)*np.pi/180
      r1=ri
      for j in range(len(Ag)-1):
        r2=E[i,j]/n[i] *b[i] +r1

        x=np.array(r1*np.cos(t[0]))
        x=np.append(x,r2*np.cos(t))
        x=np.append(x,r1*np.cos(np.flipud(t)))

        y=np.array(r1*np.sin(t[0]))
        y=np.append(y,r2*np.sin(t))
        y=np.append(y,r1*np.sin(np.flipud(t)))

      #if iflip, jcor=length(Ag)-1-j+1;
      #else, jcor=j;
      #end

        if E[i,j]>0: pl.fill(x,y,facecolor=cor[j],linewidth=LineWidth)
        r1=r2;



  # N S E W labels:
  if NSEWlab:
    bg='none';
    args={'backgroundcolor':bg,'fontsize':FontSize}
    pl.text(-g-ri, 0,'WEST', verticalalignment='top',   horizontalalignment='left', **args);
    pl.text( g+ri, 0,'EAST', verticalalignment='top',   horizontalalignment='right',**args);
    pl.text( 0,-g-ri,'SOUTH',verticalalignment='bottom',horizontalalignment='left', **args);
    pl.text( 0, g+ri,'NORTH',verticalalignment='top',   horizontalalignment='left', **args);


  # scale legend:
  L=(g*rl-g-ri)/7
  h=(g+ri)/10
  dy=h/3

  x0=g+ri+(g*rl-g-ri)/7
  x1=x0+L
  y0=-g-ri

  if legType==1: # contimuous
    for j in range(len(Ag)-1):
      lab=str(Ag[j])
      if j==0 and hasL and IncHiLow:
        lab=''

      y1=y0+h
      pl.fill([x0, x1, x1, x0],[y0, y0, y1, y1],facecolor=cor[j],linewidth=LineWidth)
      if legTickLab: pl.text(x1+L/4,y0,lab,verticalalignment='center',fontsize=FontSize)
      y0=y1

    if legTickLab and not (hasH and IncHiLow):
      pl.text(x1+L/4,y0,str(Ag[-1]),verticalalignment='center',fontsize=FontSize)

  else:
     for j in range(len(Ag)-1):
      lab=str(Ag[j])+ ' - '+ str(Ag[j+1])
      if j==0 and hasL and  IncHiLow:
        lab='<'+str(Ag[1])

      if j==len(Ag)-2 and hasH and IncHiLow:
        lab='>='+str(Ag[j])

      y1=y0+h
      pl.fill([x0, x1, x1, x0],[y0+dy, y0+dy, y1, y1],facecolor=cor[j],linewidth=LineWidth)
      if legTickLab: pl.text(x1+L/4,(y0+dy+y1)/2,lab,verticalalignment='center',fontsize=FontSize)
      y0=y1


  # title and legend label:
  x=np.mean([-g*rs,g*rl])
  y=np.mean([g+ri,g*rs])
  pl.text(x,y,titStr,horizontalalignment='center',fontsize=FontSize)

  x=x0
  y=y1+dy
  pl.text(x,y,legStr,horizontalalignment='left',verticalalignment='bottom',fontsize=FontSize)


  #pl.show()
  pl.axis('image')
  pl.axis('off')
  return fig


def plot_ellipse(x,y,major,minor,inc,phase,scale='auto',**kargs):
  '''
  x = np.linspace(1,100,8)
  y = np.linspace(1,100,10)
  x,y = np.meshgrid(x,y)
  major = (x+y)/30.
  minor = np.sqrt(x+y)/10.
  inc   = x/2.
  phase = y

  plot_ellipse(x,y,major,minor,inc,phase)

  options:
    ax, gca()
    type: 0,1,2, different ways to show phase and orientation
    color

  mma 2013
  '''


  if 'ax' in kargs.keys():
    ax=kargs['ax']
    newax=False
  else:
    if len(pl.gcf().axes): newax=False
    else: newax=True
    ax=pl.gca()

  if not (isinstance(x,np.ndarray) or np.ma.isMA(x)):
    x=np.asarray([x])
    y=np.asarray([y])
    major=np.asarray([major])
    minor=np.asarray([minor])
    inc=np.asarray([inc])
    phase=np.asarray([phase])

  def ellipse(major,minor,inc,phase,pos=(0,0),**kargs):
    def cart2pol(x,y):
      return np.angle(x+1j*y),np.abs(x+1j*y)
    def pol2cart(a,r):
      return r*np.cos(a),r*np.sin(a)

    type=kargs.get('type',0)
    color=kargs.get('color','k')

   # ellipse:
    t=np.linspace(0,2*np.pi,100)
    x=major*np.cos(t-phase*np.pi/180)
    y=minor*np.sin(t-phase*np.pi/180)
    th,r=cart2pol(x,y)
    x,y=pol2cart(th+inc*np.pi/180,r)

    if type in (0,1): a=major/10.
    elif type==2: a=major/20.

    Fi=30 #angle
    teta=np.arctan2(y[-1]-y[-2],x[-1]-x[-2])*180/np.pi
    fi=(180-Fi+teta)*np.pi/180
    fii=(180+Fi+teta)*np.pi/180

    if type==2:
      aux=2 # extra arrow point
      teta_aux=np.arctan2(y[aux]-y[aux-1],x[aux]-x[aux-1])*180/np.pi
      fi_aux=(180-Fi+teta_aux)*np.pi/180
      fii_aux=(180+Fi+teta_aux)*np.pi/180

      aux2=5 # extra arrow point
      teta_aux2=np.arctan2(y[aux2]-y[aux2-1],x[aux2]-x[aux2-1])*180/np.pi
      fi_aux2=(180-Fi+teta_aux2)*np.pi/180
      fii_aux2=(180+Fi+teta_aux2)*np.pi/180


    Fi=Fi*np.pi/180

    P1=[x[0], y[0]]
    P2=[P1[0]+a*np.cos(fi)/np.cos(Fi),  P1[1]+a*np.sin(fi)/np.cos(Fi),  P1[-1]]
    P3=[P1[0]+a*np.cos(fii)/np.cos(Fi), P1[1]+a*np.sin(fii)/np.cos(Fi), P1[-1]]

    if type==2:
      P1_aux=[x[aux], y[aux]]
      P2_aux=[P1_aux[0]+a*np.cos(fi_aux)/np.cos(Fi),  P1_aux[1]+a*np.sin(fi_aux)/np.cos(Fi),  P1[-1]]
      P3_aux=[P1_aux[0]+a*np.cos(fii_aux)/np.cos(Fi), P1_aux[1]+a*np.sin(fii_aux)/np.cos(Fi), P1[-1]]

      P1_aux2=[x[aux2], y[aux2]]
      P2_aux2=[P1_aux2[0]+a*np.cos(fi_aux2)/np.cos(Fi),  P1_aux2[1]+a*np.sin(fi_aux2)/np.cos(Fi),  P1[-1]]
      P3_aux2=[P1_aux2[0]+a*np.cos(fii_aux2)/np.cos(Fi), P1_aux2[1]+a*np.sin(fii_aux2)/np.cos(Fi), P1[-1]]

    if type==0:
      x=pos[0]+np.concatenate(([0],x,[np.nan, P2[0], P1[0], P3[0]]))
      y=pos[1]+np.concatenate(([0],y,[np.nan, P2[1], P1[1], P3[1]]))
    elif type==1:
      x=pos[0]+np.concatenate((x,[np.nan, P2[0], P1[0], P3[0]]))
      y=pos[1]+np.concatenate((y,[np.nan, P2[1], P1[1], P3[1]]))
    elif type==2:
      x=pos[0]+x
      y=pos[1]+y

      x1=[x[0]]
      y1=[y[0]]

      x_aux=pos[0]+np.asarray([P2_aux[0], P1_aux[0], P3_aux[0]])
      y_aux=pos[1]+np.asarray([P2_aux[1], P1_aux[1], P3_aux[1]])

      x_aux2=pos[0]+np.asarray([P2_aux2[0], P1_aux2[0], P3_aux2[0]])
      y_aux2=pos[1]+np.asarray([P2_aux2[1], P1_aux2[1], P3_aux2[1]])

      x=np.concatenate((x,[np.nan],x_aux, [np.nan], x_aux2))
      y=np.concatenate((y,[np.nan],y_aux, [np.nan], y_aux2))

    res=[]
    res+=[pl.matplotlib.lines.Line2D(x,y,color=color)]
    if type == 2:
      res+=[pl.matplotlib.lines.Line2D(x1,y1,color=color,marker='.')]

    return res

  # scale:
  if x.size==1: scale=1
  else:
    if scale=='auto':
      L=np.sqrt((x.max()-x.min())**2+(y.max()-y.min())**2)
      if x.ndim==1: ny,nx=1,x.size
      else: ny,nx=x.shape
      n=np.sqrt(ny**2+nx**2)
      L1=L/n # ~width of 1 ellipse
      scale=float(L1)/major.mean()
      # remove 1/3:
      scale=scale*2/3.
      print 'scale=',scale

  major=major*scale
  minor=minor*scale

  res=[]
  major=major.ravel()
  minor=minor.ravel()
  inc=inc.ravel()
  phase=phase.ravel()
  x=x.ravel()
  y=y.ravel()
  for i in range(x.size):
    res+=ellipse(major[i],minor[i],inc[i],phase[i],pos=(x[i],y[i]),**kargs)

  [ax.add_artist(i) for i in res]

  if newax:
    # not working !!
    #ax.relim()
    #ax.autoscale_view(True,True,True)
    # then:
    xmin=np.inf
    xmax=-np.inf
    ymin=np.inf
    ymax=-np.inf
    for i in res:
      xd=i.get_xdata()
      yd=i.get_ydata()
      if np.size(xd)<=1: continue
      xmin=np.min([xd[~np.isnan(xd)].min(),xmin])
      xmax=np.max([xd[~np.isnan(xd)].max(),xmax])
      ymin=np.min([yd[~np.isnan(yd)].min(),ymin])
      ymax=np.max([yd[~np.isnan(yd)].max(),ymax])

    ax.axis([xmin,xmax,ymin,ymax])

