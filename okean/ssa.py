from numpy import zeros, array, outer, rot90, mean, diag, unique, any, sum
import numpy as np
from numpy.linalg import svd

import pylab as pl


class ssa:
  '''
  Singular Spectrum Analysis

  Methods/Usage:
    1. s=ssa(<filename>)
    2. s.truncate(iStart,iEnd)
    3. s.calc(winWidth)
      3.1 s.embed(winWidth)
      3.2 s.decomp()
      3.3 s.group()
      3.4 s.diagavg()
    4. s.reconst(eigentriples)
    5. s.forecast(dataType,startPoint,nPoints)
    6. plots:
      6.1 data: s.plot('data')
      6.2 reconstructed: s.plot('rec')
      6.3 residual: s.plot('resid')
      6.4 eigenvector left:      s.plot('evl',inds)
      6.5 eigenvector right:     s.plot('evr',inds)
      6.6 dispersion plot left:  s.plot('dpl',inds)
      6.7 dispersion plot right: s.plot('dpr',inds)
      6.8 forecast: s.plot('forec')
  '''


  def __init__(self,fname):
    '''
    Input:
      fname, series text file
    '''

    if isinstance(fname,np.ndarray):
      self.fname=''
      self.__data=fname
      self.range=0,len(self.__data)
      self.dtype=self.__data.dtype
    else:
      self.fname=fname
      self.load()


  def load(self):
    '''
    Load text data file to be used by SSA
    '''
    s=open(self.fname).readlines()
    self.__data=array(map(float,s))
    self.range=0,len(self.__data)
    self.dtype=self.__data.dtype


  def data(self):
    '''
    Returns original SSA data series (may be truncated)
    '''
    
    return self.__data[self.range[0]:self.range[1]]


  def truncate(self,i1=0,i2=-1):
    '''
    Truncate the data series

    Inputs:
      i1, i2, first and last data index
    '''
    
    self.range=i1,i2


  def embed(self,win):
    '''
    SSA embedding step: create the trajectory matrix

    Input:
      win, window width (<= data length /2)
    '''

    data=self.data()
    N=len(data)
    if win>N/2:
      print('error: win width max = N/2 = ',N/2, ' N=',N)
      return

    if win<2:
      print('error: win width min = 2')
      return

    K=N-win+1
    X=zeros((win,K),self.dtype)
    for j in range(K):
      for i in range(win):
         X[i,j]=data[i+j]
     
    self.win=win
    self.traject=X


  def decomp(self):
    '''
    SSA singular value decomposition step:
    create the left and right singular vector matrice and the
    diagonal singular values matrix    
    '''
    
    U,Lambd,V=svd(self.traject)
    self.svd_U=U
    self.svd_lambda=Lambd
    self.svd_V=V.T # here differs from matlab !!!


  def group(self):
    '''
    SSA grouping step: create the matrix A [ L K L]
    '''
    
    L=self.win
    K=self.svd_V.shape[0]

    Last  = min(L,K) # L, always
    Kast  = max(L,K) # K, always

    A=zeros((Last,Kast,Last),self.dtype)
    for i in range(Last):
      A[:,:,i] = self.svd_lambda[i]*outer(self.svd_U[:,i],self.svd_V[:,i])

    self.group_A=A


  def diagavg(self):
    '''
    SSA diagonal averaging reconstruction step:
    create the matrix G [ Nt L]
    '''

    L,K=self.group_A.shape[:-1]
    Last  = min(L,K)
    Kast  = max(L,K)
    F=self.data()
    Nt=len(F)

    G=zeros((Nt,Last),self.dtype)
    for i in range(Last):
      aux=self.group_A[:,:,i]
      if L<K: # always !
        aux=rot90(aux)
      else:
        aux=rot90(aux.T)

      j=-1
      for n in range(-Kast+1,Last):
        j+=1
        G[j,i] = mean(diag(aux,n))

    self.diagAvg=G


  def calc(self,win):
    '''
    Perform required SSA steps:
    embed, decomp, group and diagavg

    Input:
      win, window width, input of the embed method

    Example:
       obj=ssa('file.dat')
       win=20

       # The four steps:
       ob.embed(win)
       ob.decomp()
       ob.group()
       ob.diagavg()

       # Can be done as:
       ob=calc(win)
    '''
    self.embed(win)
    self.decomp()
    self.group()
    self.diagavg()


  def reconst(self,eigentriples):
    '''
    SSA reconstruction step

    Input:
      eigentriples, reconstruction eigentriples indice

    Example:
      obj.reconst([1,2,3])
    '''
	  
    L,K=self.group_A.shape[:-1]
    Last  = min(L,K)
    Kast  = max(L,K)
    F=self.data()
    Nt=len(F)

    eigentriples=unique(eigentriples)
    if any(map(lambda i: i>Last or i<0, eigentriples)):
       print('error: eigentriples indice must be < %d' % Last)
       return

    auxi=self.diagAvg[:,eigentriples]
    FR    = sum(auxi,1)
    Resid = F-FR

    self.eigentriples = eigentriples
    self.reconstr     = FR
    self.resid        = Resid


  def forecast(self,dataType,startPoint, nPoints):
    import numpy as np
    F=self.data()
    U=self.svd_U
    V=self.svd_U

    if dataType==0: F=self.data()
    else: F=self.reconstr

    L,K=self.group_A.shape[:-1]
    Last  = min(L,K)
    Kast  = max(L,K)

    PI=self.svd_U[Last-1,self.eigentriples]
    vert_coef = (PI**2).sum()
    Raux      = np.tile(PI,(Last-1,1) )* self.svd_U[:Last-1,self.eigentriples]
    R         = 1/(1-vert_coef) * np.sum(Raux,1)
    Ff        = np.tile(np.nan,(nPoints))

    for i in range(nPoints):
      moveWin = np.array([F[startPoint-Last+i+1:startPoint-1+1].tolist()+ Ff[max(0,i-Last+1):i-1+1].tolist()])
####      #print i,F[startPoint-Last+i:startPoint-1].shape, Ff[max(0,i-Last+1):i-1].shape
#####     print startPoint,Last, startPoint-Last+i+1,startPoint-1+1,max(0,i-Last+1),i-1
      R=np.squeeze(R)
      moveWin=np.squeeze(moveWin)
####      print i, moveWin.shape, R.shape
      Ff[i]=np.dot(R,moveWin)

    self.fff=Ff


  def plot(self,what='data',inds=False):
    '''
    Plot data, reconstructed and residual series
    Plot eigenvector or dispersion plots from the left and right
    singular vector matrice

    Inputs:
      what,  data, for original (or truncated) series
             rec, for reconstructed series
	     resid, for residual series

	     evl, for eigenvector left plot
	     evr, for eigenvector right plot
	     dpl, for dispersion plot left
	     dpr, for dispersion plot right

	     forc, for forecast data 
	     
      inds, indice for eigenvector and dispersion plots

    Example:
      obj.plot(evr,[1,2,3])       
    '''
    
    pl.figure()

    if inds:
      try: iter(inds)
      except: inds=[inds]

    if what=='data':
      pl.plot(self.data())
    elif what=='rec':
      pl.plot(self.data())
      pl.plot(self.reconstr)
    elif what=='resid':
      pl.plot(self.data())
      pl.plot(self.resid)
    elif what=='evl': # eigenvectors left
      #pl.plot(self.svd_U[:,inds])
      for i in inds:
        pl.plot(self.svd_U[:,i],'.-')
	
      pl.title('Eigenvector (left) %s' % str(inds))
      
    elif what=='evr': # eigenvectors right
      #pl.plot(self.svd_V[:,inds])
      for i in inds:
        pl.plot(self.svd_V[:,i],'.-')
	
      pl.title('Eigenvector (right) %s' % str(inds))
    
    elif what=='dpl': # dispersion plot left
      for i in inds:
        pl.plot(self.svd_U[:,i],self.svd_U[:,i+1],'.-')
	
      pl.title('Disperion plot (left) %s' % str(inds))
      
    elif what=='dpr': # dispersion plot right
      for i in inds:
        pl.plot(self.svd_V[:,i],self.svd_V[:,i+1],'.-')

      pl.title('Disperion plot (right) %s' % str(inds))

    elif what=='forec':
	    pass # TODO !!!!!!
