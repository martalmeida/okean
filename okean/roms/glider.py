from okean import roms, dateu, netcdf, calc
import numpy as np


class RomsGlider:
  def __init__(self,f,x,y,t,quiet=1):
    self.fname=f
    self.x=x
    self.y=y
    self.t=t
    if not quiet: print 'loading ',f
    self.roms=roms.His(f)
    self.nc=self.roms.nc

    self.inds0={}
    self.uinds={}

  def plot(self,**kargs):
    o=self.roms.grid.plot(**kargs)
    x,y=o.map(self.x,self.y)
    o.ax.plot(x,y,'r')
    o.ax.plot(x[0],y[0],'bo')
    o.ax.plot(x[-1],y[-1],'ro')

  def model_lon_lat(self,varname):
    return self.roms.grid.vars(ruvp=varname)[:2]

  def horiz_var_type(self,varname):
    tag=self.roms.var_at(varname)
    #if self.roms.hast(varname): tag=tag+'_t'
    return tag

  def locate(self,varname,quiet=1):

    # find points inside grid:
    inds=self.roms.grid.ingrid(self.x,self.y)
    inds=np.where(inds)[0]

    # lon, lat:
    lon,lat=self.model_lon_lat(varname)

    ntimes=len(self.t)
    V=np.ma.array(())

    Iind,Jind,Tind=[],[],[]
    INDS=[]
    for c in range(ntimes):
      if not quiet: print '  %s   %d of %d' % (self.t[c].isoformat(' '),c,ntimes)
      if not c in inds:
        if not quiet: print ' - outside domain in time %s'%t[c].isoformat()
        INDS+=[(-1,-1,-1)]
        continue

      # find time indices:
      if self.t[c]>=self.roms.time[0] and self.t[c]<=self.roms.time[-1]:
        i0=np.where(self.roms.time<=self.t[c])[0][-1]
      else: i0=-1

      # find i,j bbox:
      d=(lon-self.x[c])**2+(lat-self.y[c])**2
      i,j=np.where(d==d.min())
      i=i[0]
      j=j[0]
      isIn=False
      for I in [i-1,i]:
        if I+1>=lon.shape[0] or I<0: continue
        for J in [j-1,j]:
          if J+1>=lon.shape[1] or J<0: continue

          xp=[lon[I,J],lon[I,J+1],lon[I+1,J+1],lon[I+1,J]]
          yp=[lat[I,J],lat[I,J+1],lat[I+1,J+1],lat[I+1,J]]
          if calc.inpolygon(self.x[c],self.y[c],xp,yp):
            i,j=I,J
            isIn=True
            break
        else:
          # Continue if the inner loop wasn't broken.
          continue
        # Inner loop was broken, break the outer.
        break

      INDS+=[(i,j,i0)]

      Iind+=[ i,   i, i+1, i+1 ]
      Jind+=[ j, j+1, j+1,   j ]
      Tind+=[i0,  i0,  i0,  i0 ]

      if i0<0: ii0=i0
      else: ii0=i0+1

      Iind+=[  i,    i,  i+1, i+1 ]
      Jind+=[  j,  j+1,  j+1,   j ]
      Tind+=[ii0,  ii0,  ii0, ii0 ]


    Iind=np.asarray(Iind)
    Jind=np.asarray(Jind)
    Tind=np.asarray(Tind)
    INDS=np.asarray(INDS)

    # deal with points at boundary ! just in case...
    imax,jmax=[k-1 for k in lon.shape]
    tmax=self.roms.time.size-1
    Iind=np.where(Iind>imax,imax,Iind)
    Jind=np.where(Jind>jmax,jmax,Jind)
    Tind=np.where(Tind>tmax,tmax,Tind)

    inds=[(Iind[i],Jind[i],Tind[i]) for i in range(len(Iind))]
    uinds=np.asarray(list(frozenset(inds)))

    tag=self.horiz_var_type(varname)
    self.inds0[tag]=INDS
    self.uinds[tag]=uinds

    return INDS, uinds


  def extract(self,varname,method='fast',nfast=1,quiet=1):
    if not self.roms.hast(varname):
      method,nfast='fast',1

    tag=self.horiz_var_type(varname)
    if tag in self.uinds.keys():
      inds0=self.inds0[tag]
      uinds=self.uinds[tag]
    else:
      inds0,uinds=self.locate(varname,quiet)

    J,I,T = uinds.T
    J0,I0,T0 = inds0.T

    if not self.roms.hast(varname): no_time=True
    else: no_time=False


    if method=='fast': # faster but will download more data than needed!
      II=range(I.min(),I.max()+1)
      JJ=range(J.min(),J.max()+1)
      TT=range(T[T>=0].min(),T.max()+1)


      if not calc.isiterable(nfast) and nfast==1:
        if not quiet: print 'loading %s : ijt= %d %d %d'%(varname,len(II),len(JJ),len(TT))
        v=netcdf.use(self.nc,varname,xiSEARCH=II,etaSEARCH=JJ,ocean_time=TT)
      else:
        v=False
        if not calc.isiterable(nfast):
          nfast=min(nfast,len(TT)-1)
          tmp=range(TT[0],TT[-1]+2,(TT[-1]-TT[0])/nfast)
          if tmp[-1]<TT[-1]+1: tmp+=[TT[-1]+1]
        else: tmp=nfast

        for k in range(len(tmp)-1):
          tt=range(tmp[k],tmp[k+1])
          if not quiet: print 'loading %s : ijt= %d %d %d (t %d to %d)'%(varname,len(II),len(JJ),len(tt), tt[0],tt[-1])
          vtmp=netcdf.use(self.nc,varname,xiSEARCH=II,etaSEARCH=JJ,ocean_time=tt)
          if not v is False:
            if vtmp.ndim<v.ndim: vtmp=vtmp[np.newaxis,:] # if len of last tt is 1 !
            v=np.vstack((v,vtmp))
          else: v=vtmp

      
      if v.ndim>3:
        V=np.ma.zeros((I.size,v.shape[1]),'f')
      else:
        V=np.ma.zeros(I.size,'f')

      for i in range(I.size):
        xi   = I[i]-I.min()
        eta  = J[i]-J.min()
        if T[i]>=0: tind = T[i]-T[T>=0].min()
        else: tind=T[i] # negative means data ouside model time

        if v.ndim==4:
          V[i]=v[tind,:,eta,xi]
        elif v.ndim==3:
          V[i]=v[tind,eta,xi]
        elif v.ndim==2:
          V[i]=v[eta,xi]

    else:
      V=False
      for i in range(len(I)):
        if T[i]<0: continue
        if not quiet: print 'loading %s  (%d of %d): ijt= %d %d %d'%(varname,i,len(I),I[i],J[i],T[i])
        v=netcdf.use(self.nc,varname,xiSEARCH=I[i],etaSEARCH=J[i],ocean_time=T[i])
        if V is False:
          if v.ndim>1: shape=(len(I),)+v.shape
          else: shape=len(I)
          V=np.ma.zeros(shape,'f')

        V[i]=v

    lon,lat=self.model_lon_lat(varname)
    U=np.array(())
    for i in range(len(self.t)):
      xi   = I0[i]
      eta  = J0[i]
      tind = T0[i]
      if tind<0: continue # data outside model time

      # rotate cell before interp:
      xp=np.asarray([lon[eta,xi],lon[eta,xi+1],lon[eta+1,xi+1],lon[eta+1,xi],self.x[i]])
      yp=np.asarray([lat[eta,xi],lat[eta,xi+1],lat[eta+1,xi+1],lat[eta+1,xi],self.y[i]])
      xp,yp=calc.rot2d(xp,yp,self.roms.grid.angle[eta,xi])
      x=xp[:-1].min(),xp[-1],xp[:-1].max()
      y=yp[:-1].min(),yp[-1],yp[:-1].max()

      A = x[1]-x[0]
      a = x[2]-x[1]
      B = y[1]-y[0]
      b = y[2]-y[1]

      tcond=(T==tind)#|(tind<0)
      tcond1=(T==tind+1)#|(tind+1<0)
      j0=np.where((I==xi)&(J==eta)&tcond)[0][0]
      j1=np.where((I==xi+1)&(J==eta)&tcond)[0][0]
      j2=np.where((I==xi+1)&(J==eta+1)&tcond)[0][0]
      j3=np.where((I==xi)&(J==eta+1)&tcond)[0][0]
      u0=(V[j0]*a*b+V[j1]*A*b+V[j2]*A*B+V[j3]*a*B)/(a*b+A*b+A*B+a*B)
      
      if not no_time:
        dt0=self.t[i]-self.roms.time[tind]
        dt1=self.roms.time[tind+1]-self.t[i]
        dt0=dt0.days*86400+dt0.seconds
        dt1=dt1.days*86400+dt1.seconds

        k0=np.where((I==xi)&(J==eta)&tcond1)[0][0]
        k1=np.where((I==xi+1)&(J==eta)&tcond1)[0][0]
        k2=np.where((I==xi+1)&(J==eta+1)&tcond1)[0][0]
        k3=np.where((I==xi)&(J==eta+1)&tcond1)[0][0]

        u1=(V[k0]*a*b+V[k1]*A*b+V[k2]*A*B+V[k3]*a*B)/(a*b+A*b+A*B+a*B)
        u=(u0*dt1+u1*dt0)/(dt1+dt0)
      else:
        u=u0

      if not U.size:
        U=np.ma.zeros((len(self.t),u.size),'f')
        U=np.ma.masked_where(U==0,U)

      U[i]=u

    U=np.squeeze(U)
    return U


  def depth(self,varname='r',quiet=1):
    h=self.extract('h',quiet=quiet)
    zeta=self.extract('zeta',quiet=quiet)
    z,zw=roms.roms_tools.s_levels(h,zeta,self.roms.s_params)
    isW=varname[0]=='w'
    if isW:z=zw
    return z.T
