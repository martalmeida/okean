import pylab as pl
import numpy as np


def fillout(x,y,lims,color=False,**kargs):
  '''
  Fill outside polygon.
  Fill the region between the polygon x,y and the rectangle lims
  (x0,x1,y0,y1). If no color is provided, the final polygon (ready to
  be used with pl.fill) is returned. Otherwise, pl.fill is done
  using alguments color and kargs. In this case, will be used the
  current axes or the one provided by key ax.

  Example:
    x=np.array([0,1,1,0.5])
    y=np.array([0,-.1,1.1,1])
    fillout(x,y,[-2,2,-2,2],'g',edgecolor='none')
  '''

  from okean.calc import poly_area
  if poly_area(x,y)<0:
    x,y=x[::-1],y[::-1]

  xi,xe,yi,ye=lims
  i=np.where(x==x.min())[0][0]

  x=np.concatenate([x[i:], x[:i+1]])
  y=np.concatenate([y[i:], y[:i+1]])

  tmpx=np.array([xi,   xi, xe, xe, xi, xi,   x[0]])
  tmpy=np.array([y[0], ye, ye, yi, yi, y[0], y[0]])

  x=np.concatenate([tmpx, x])
  y=np.concatenate([tmpy, y])

  if color is False: return x,y
  else:
    ax=kargs.pop('ax',pl.gca())
    o=ax.fill(x,y,color,**kargs)
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

  D=np.asarray(D)
  F=np.asarray(F)

  ax          = kargs.get('ax',          pl.gca())
  bg          = kargs.get('bg',          'w')
  dtype       = kargs.get('dtype',       'standard')
  nAngles     = kargs.get('n',           36)
  ri          = kargs.get('ri',          1/30.)
  quad        = kargs.get('quad',        0)
  legType     = kargs.get('legtype',     2)
  percBg      = kargs.get('percbg',      'w')
  titStr      = kargs.get('title',       '')
  legStr      = kargs.get('legend',      '')
  cmap        = kargs.get('cmap',        None)
  colors      = kargs.get('colors',      [])
  Ag          = kargs.get('di',          False) # intensity subdivs
  ci          = kargs.get('ci',          []) # percentage circles
  lineColors  = kargs.get('lcolor',      'k')
  borderColor = kargs.get('bcolor',      False)
  iflip       = kargs.get('iflip',       0)
  inorm       = kargs.get('inorm',       0)
  IncHiLow    = kargs.get('incout',      1) # include values higher and lower that the limits of Ag
  FontSize    = kargs.get('fontsize',    8)
  LineWidth   = kargs.get('linewidth',   0.5)
  legTickLab  = kargs.get('lablegticks', True)
  NSEWlab     = kargs.get('NSEWlab',     True)
  labels      = kargs.get('labels',      True)
  percLabels  = kargs.get('percLabels',  True)

  if not labels:
    legTickLab=NSEWlab=percLabels=False

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
    from okean import ticks
    Ag,sAg=ticks.loose_label_n(F.min(),F.max(),ntick=7,labels=True)
  else:
    sAg=['%g'%i for i in Ag]

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
  ri=g*ri
  if bg: ax.fill([-rs*g, rl*g, rl*g, -rs*g],[-rs*g, -rs*g, rs*g, rs*g],
                 bg,edgecolor='w',linewidth=LineWidth,label='wr_bg')

  ax.plot([-g-ri, -ri, np.nan, ri, g+ri, np.nan, 0, 0, np.nan, 0, 0],
          [0, 0, np.nan, 0, 0, np.nan, -g-ri, -ri, np.nan, ri, g+ri],
          ':',color=lineColors,linewidth=LineWidth)

  t0=np.arange(0,361)*np.pi/180
  labs=[]
  Ang=np.array([1, 3, 5, 7])/4.*np.pi

  Valign='top', 'top', 'bottom', 'bottom'
  Halign='center', 'left', 'left', 'right'

  for i in range(ncircles):
    x=(ci[i]+ri)*np.cos(t0)
    y=(ci[i]+ri)*np.sin(t0)

    ax.plot(x,y,':',color=lineColors,linewidth=LineWidth);

    if percLabels: ax.text((ci[i]+ri)*np.cos(Ang[quad]),(ci[i]+ri)*np.sin(Ang[quad]),str(ci[i])+'%',
      verticalalignment=Valign[quad],horizontalalignment=Halign[quad],
      backgroundcolor=percBg,fontsize=FontSize)


  # calc colors:
  if not colors:
    cor=list(range(len(Ag)-1))
    q=pl.cm.ScalarMappable(cmap=cmap)
    q.set_clim([Ag[0], Ag[-2]])
    for j in range(len(Ag)-1):
      cor[j]=q.to_rgba(Ag[j])
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

        if E[i,j]>0: ax.fill(x,y,facecolor=cor[j],linewidth=LineWidth)
        r1=r2;


  # N S E W labels:
  if NSEWlab:
    bg='none';
    args={'backgroundcolor':bg,'fontsize':FontSize}
    ax.text(-g-ri, 0,'WEST', verticalalignment='top',   horizontalalignment='left', **args);
    ax.text( g+ri, 0,'EAST', verticalalignment='top',   horizontalalignment='right',**args);
    ax.text( 0,-g-ri,'SOUTH',verticalalignment='bottom',horizontalalignment='left', **args);
    ax.text( 0, g+ri,'NORTH',verticalalignment='top',   horizontalalignment='left', **args);


  # scale legend:
  L=(g*rl-g-ri)/7.
  h=(g+ri)/10.
  dy=h/3

  x0=g+ri+(g*rl-g-ri)/7.
  x1=x0+L
  y0=-g-ri

  if legType==1: # contimuous
    for j in range(len(Ag)-1):
      #lab=str(Ag[j])
      lab=sAg[j]
      if j==0 and hasL and IncHiLow:
        lab=''

      y1=y0+h
      ax.fill([x0, x1, x1, x0],[y0, y0, y1, y1],facecolor=cor[j],linewidth=LineWidth)
      if legTickLab: ax.text(x1+L/4,y0,lab,verticalalignment='center',fontsize=FontSize)
      y0=y1

    if legTickLab and not (hasH and IncHiLow):
      ax.text(x1+L/4,y0,str(Ag[-1]),verticalalignment='center',fontsize=FontSize)

  elif legType==2:
     for j in range(len(Ag)-1):
      #lab=str(Ag[j])+ ' - '+ str(Ag[j+1])
      lab=sAg[j]+ ' - '+ sAg[j+1]
      if j==0 and hasL and  IncHiLow:
        #lab='<'+str(Ag[1])
        lab='<'+sAg[1]

      if j==len(Ag)-2 and hasH and IncHiLow:
        #lab='>='+str(Ag[j])
        lab='>='+sAg[j]

      y1=y0+h
      ax.fill([x0, x1, x1, x0],[y0+dy, y0+dy, y1, y1],facecolor=cor[j],linewidth=LineWidth)
      if legTickLab: ax.text(x1+L/4,(y0+dy+y1)/2,lab,verticalalignment='center',fontsize=FontSize)
      y0=y1


  # title and legend label:
  x=np.mean([-g*rs,g*rl])
  y=np.mean([g+ri,g*rs])
  ax.text(x,y,titStr,horizontalalignment='center',fontsize=FontSize)

  if legType in [1,2]:
    x=x0
    y=y1+dy
    ax.text(x,y,legStr,horizontalalignment='left',verticalalignment='bottom',fontsize=FontSize)


  ax.axis('equal') #'image')
  ax.axis('off')



def plot_ellipse(x,y,major,minor,inc,phase,scale='auto',fill={},**kargs):
  '''
  Draw ellipses with semi-major axis, semi-minor axis, inclination
  and phase. Deppending on type, ellipses orientations and phase  can be
  illustrated (in different ways).

  Options:
    scale: 'auto' by default or some number
    ax, gca()
    type: 0,1,2,3,4 different ways to show phase and orientation
    line options (color, lw, etc)
    marker options for type 2 (dict inside kargs with key marker). Ex:
      plot_ellipse(...,**{'marker':{'ms':20,'marker':'*'}})
    default color is 'k'
    default marker is '.'
    arrow_scale, arrow sclae for type 0, 1 and 2. Default is 1/10 (of major)
    for type 0,1 and 1/20 for type 2.
    arrow_angle, angle of arrow used in types 0, 1 and 2 (default is 30).
    relim: True, recalc axis limits (True if new axes is created)
    fill: fill ellipse, use this argument to provide fill color, alpha, etc
          Ex: ...fill={'color','b',alpha=.5}. Default: facecolor: 'b',
          lw:0, alpha:.5

  Examples:
    x = np.linspace(1,100,8)
    y = np.linspace(1,100,10)
    x,y = np.meshgrid(x,y)
    major = (x+y)/30.
    minor = np.sqrt(x+y)/10.
    inc   = x/2.
    phase = y

    pl.figure()
    plot_ellipse(x,y,major,minor,inc,phase)
    pl.figure()
    opts=dict(color='g',marker=dict(color='r',marker='s',ms=5))
    plot_ellipse(x,y,major,minor,inc,phase,type=2,**opts)


  mma 2013
  '''

  type=kargs.get('type',0)
  relim=kargs.get('relim',True)
  try: kargs.pop('type')
  except: pass
  try: kargs.pop('relim')
  except: pass

  if not 'color' in kargs.keys(): kargs['color']='k'

  # marker options for type 2:
  if 'marker' in kargs.keys():
    mkargs=kargs.pop('marker')
    if not 'color'  in mkargs.keys(): mkargs['color']=kargs['color']
    if not 'marker' in mkargs.keys(): mkargs['marker']='.'
  else:
    mkargs={'color': kargs['color'],'marker':'.'}

  if 'ax' in kargs.keys():
    ax=kargs.pop('ax')
    newax=False
  else:
    if len(pl.gcf().axes): newax=False
    else: newax=True
    ax=pl.gca()

  if 'arrow_scale' in kargs.keys():
    sarrow_01=sarrow_2=kargs['arrow_scale']
    kargs.pop('arrow_scale')
  else:
    sarrow_01=1/10.
    sarrow_2=1/20.

  aarrow=kargs.pop('arrow_angle',30)

  if not (isinstance(x,np.ndarray) or np.ma.isMA(x)):
    x=np.asarray([x])
    y=np.asarray([y])
    major=np.asarray([major])
    minor=np.asarray([minor])
    inc=np.asarray([inc])
    phase=np.asarray([phase])

  def ellipse(major,minor,inc,phase,pos=(0,0)):
    def cart2pol(x,y):
      return np.angle(x+1j*y),np.abs(x+1j*y)
    def pol2cart(a,r):
      return r*np.cos(a),r*np.sin(a)

    # ellipse:
    t=np.linspace(0,2*np.pi,100)
    x=major*np.cos(t-phase*np.pi/180)
    y=minor*np.sin(t-phase*np.pi/180)
    th,r=cart2pol(x,y)
    x,y=pol2cart(th+inc*np.pi/180,r)

    if type in (0,1): a=major*sarrow_01
    elif type==2: a=major*sarrow_2
    else: a=0

    Fi=aarrow #angle
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

    x_fill=pos[0]+x
    y_fill=pos[0]+y
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
    elif type==3:
      x=pos[0]+x
      y=pos[1]+y
    elif type==4:
      x=pos[0]+np.concatenate(([0],x))
      y=pos[1]+np.concatenate(([0],y))

    res={}
    res['line']=pl.matplotlib.lines.Line2D(x,y,**kargs)
    if type == 2:
      res['marker']=pl.matplotlib.lines.Line2D(x1,y1,**mkargs)

    if fill:
      if not 'lw'        in fill.keys(): fill['lw']        = 0
      if not 'facecolor' in fill.keys(): fill['facecolor'] = 'b'
      if not 'alpha'     in fill.keys(): fill['alpha']     = .5
      res['fill']=pl.matplotlib.patches.Ellipse(pos,2*major,2*minor,inc,**fill)

    return res

  # scale:
  if x.size==1 and scale=='auto': scale=1
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
      #print 'scale=',scale

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
    res+=[ellipse(major[i],minor[i],inc[i],phase[i],pos=(x[i],y[i]))]

  for k in res[0].keys():
    [ax.add_artist(i[k]) for i in res]

  if newax: relim=True
  if relim:
    # not working !!
    #ax.relim()
    #ax.autoscale_view(True,True,True)
    # then:
    xmin=np.inf
    xmax=-np.inf
    ymin=np.inf
    ymax=-np.inf
    for i in ax.artists:#res:
      try:
        xd=i.get_xdata()
      except: continue
      yd=i.get_ydata()
      if np.size(xd)<=1: continue
      xmin=np.min([xd[~np.isnan(xd)].min(),xmin])
      xmax=np.max([xd[~np.isnan(xd)].max(),xmax])
      ymin=np.min([yd[~np.isnan(yd)].min(),ymin])
      ymax=np.max([yd[~np.isnan(yd)].max(),ymax])

    ax.axis([xmin,xmax,ymin,ymax])

  return res


def parrow(x,y,n=4,d=0,d0=0,**kargs):
  '''
  Along path arrows.
  Draws arrows connecting x,y points (1d). Each arrow will have n
  points, d points between. Use d=-1 to have arrows connected. The
  input d0 controls the beginning of the x,y path. It should be in
  the range  d+1 >= d0 > 2-n. If negative, d0 points will be subtracted
  in the first arrow.

  Options:
    head_scale, 1/5., arrow head length, relative to speed, as
      defined by head_speed
    head_speed, mean, sum or last. The speed to scale arrow head can
      be the mean of the n-1 segments, the sum, or the length of only
      the last segment
    head_angle, 30 (deg), arrow head angle
    ax, pl.gca()


  Example:
    r = np.arange(0, .075, 0.0005)
    theta = 90*np.pi*r
    y = r*np.sin(theta)
    x = r*np.cos(theta)
    parrow(x,y)
    pl.axis('equal')

  mma 2014
  '''

  head_scale = kargs.pop('head_scale',1./5)
  head_speed = kargs.pop('head_speed','mean') # sum, last
  head_angle = kargs.pop('head_angle',30) # deg
  ax=kargs.pop('ax',pl.gca())
  c=kargs.pop('C',None)

  if c is None:
    c=np.ma.zeros(x.shape,x.dtype)
    noC=True
  else: noC=False
  try:
    len(c.mask)
  except:
    c=np.ma.array(c)
    c.mask=False*c.data

  if d0:
    if d0<0:
      if n+d0<2: print('Warning: d0 is too low!') # arrow is just one point...
      x=np.hstack((x[0]+np.zeros(-d0),x))
      y=np.hstack((y[0]+np.zeros(-d0),y))
      c=np.ma.hstack((c[0]+np.zeros(-d0),c))
      c.mask[:-d0]=True
    else:
      if d0>d+1: print('Warning: d0 too high!') # a new arrow should be on the way...
      x=x[d0:]
      y=y[d0:]
      c=c[d0:]

  if d>=0:
    N=n+d
    X=np.reshape(x[:x.size/N*N],(x.size/N,N))
    Y=np.reshape(y[:y.size/N*N],(y.size/N,N))
    C=np.reshape(c[:c.size/N*N],(c.size/N,N))
    X=X[:,:n]
    Y=Y[:,:n]
    C=C[:,:n]

    nend=x.size%N
    if nend>1:
      xadd=np.zeros(X.shape[1],X.dtype)
      xadd[:nend]=x[-nend:]
      xadd[nend:]=x[-1]

      yadd=np.zeros(Y.shape[1],Y.dtype)
      yadd[:nend]=y[-nend:]
      yadd[nend:]=y[-1]

      cadd=np.zeros(C.shape[1],C.dtype)
      cadd[:nend]=c[-nend:]
      cadd[nend:]=c[-1]

      X=np.vstack((X,xadd))
      Y=np.vstack((Y,yadd))
      C=np.ma.vstack((C,cadd))
      C.mask[-1,nend:]=True

    X,Y,C=X.T,Y.T,C.T



  else: # d=-1, no space between arrows
    N=n-1
    X=np.zeros((x.size/N,n),x.dtype)
    Y=np.zeros((y.size/N,n),y.dtype)
    C=np.ma.zeros((c.size/N,n),c.dtype)

    X[:,:-1]=np.reshape(x[:x.size/N*N],(x.size/N,N))
    Y[:,:-1]=np.reshape(y[:y.size/N*N],(y.size/N,N))
    C[:,:-1]=np.reshape(c[:c.size/N*N],(c.size/N,N))

    X[:-1,-1]=X[1:,0]
    Y[:-1,-1]=Y[1:,0]
    C[:-1,-1]=C[1:,0]
    C.mask[:-1,-1]=True

    # add last val if available or remove last arrow!
    try:
      X[-1,-1]=x[x.size/N*N]
      Y[-1,-1]=y[y.size/N*N]
      C[-1,-1]=c[c.size/N*N]
    except:
      X=X[:-1,:]
      Y=Y[:-1,:]
      C=C[:-1,:]

    nend=X.size-x.size-X.shape[0]+1 # points not used
    if nend:
      xadd=np.zeros(X.shape[1],X.dtype)
      xadd[:-nend+1]=x[nend-1:]
      xadd[-nend+1:]=x[-1]

      yadd=np.zeros(Y.shape[1],Y.dtype)
      yadd[:-nend+1]=y[nend-1:]
      yadd[-nend+1:]=y[-1]

      cadd=np.zeros(C.shape[1],C.dtype)
      cadd[:-nend+1]=c[nend-1:]
      cadd[-nend+1:]=c[-1]

      X=np.vstack((X,xadd))
      Y=np.vstack((Y,yadd))
      C=np.ma.vstack((C,cadd))
      C.mask[-1,-nend+1:]=True

    X,Y,C=X.T,Y.T,C.T

  if head_speed=='sum':
    # use the whole n points speed to define arrow head:
    s=np.sqrt((X[1:]-X[:-1])**2+(Y[1:]-Y[:-1])**2).sum(0)
    a=s*head_scale
  elif head_speed=='mean':
    # or use avg instead:
    s=np.sqrt((X[1:]-X[:-1])**2+(Y[1:]-Y[:-1])**2).mean(0)
    a=s*head_scale*(n-1)
  elif head_speed=='last':
    s=np.sqrt((X[-1]-X[-2])**2+(Y[-1]-Y[-2])**2)
    a=s*head_scale*(n-1)

  if head_angle>0 and head_scale>0:
    teta=np.arctan2(Y[-1]-Y[-2],X[-1]-X[-2])
    Fi=head_angle*np.pi/180
    fi=(np.pi-Fi+teta)
    fii=(np.pi+Fi+teta)

    Xa=np.vstack((X[-1]+a*np.cos(fi)/np.cos(Fi),X[-1],X[-1]+a*np.cos(fii)/np.cos(Fi)))
    Ya=np.vstack((Y[-1]+a*np.sin(fi)/np.cos(Fi),Y[-1],Y[-1]+a*np.sin(fii)/np.cos(Fi)))
    X=np.vstack((X,Xa))
    Y=np.vstack((Y,Ya))

    if (d<0 and nend) or (d>=0 and nend>1): # d<0 means d==-1 !
      X[-3:,-1]=X[-4,-1]
      Y[-3:,-1]=Y[-4,-1]


  cm=None
  if 'cmap' in kargs:
    cm=kargs.pop('cmap')
    norm=kargs.pop('norm',pl.matplotlib.colors.Normalize)
    if noC: C=s # by default use speed to set colors
    else: C=C.mean(0)
  elif not 'color' in kargs:
    kargs['color']='k'

  p=ax.plot(X,Y,**kargs)
  if cm:
    cNorm  = norm(vmin=C.min(), vmax=C.max())
    scalarMap = pl.cm.ScalarMappable(norm=cNorm, cmap=cm)
    [p[i].set_color(scalarMap.to_rgba(C[i])) for i in range(len(p))]
    scalarMap._A=[]
    return p,scalarMap
  else: return p
