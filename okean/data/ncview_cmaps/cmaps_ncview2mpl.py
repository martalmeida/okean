'''Reads ncview colormap files (ncview/src/colormaps_*.h')
and converts to matplotlib

mma mar 2014
'''

import numpy as np
import pylab as pl

def read_file(f):
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


def gen_cmap(file,name='auto',N=None):
  '''Read ncview colormaps_<name>.h file'''

  if name=='auto': name=file.split('colormaps_')[-1][:-2]
  r=read_file(file)
  return pl.matplotlib.colors.ListedColormap(r, name=name, N=N)


def show():
  '''Display ncview colormaps'''

  a=np.outer(np.arange(0,1,0.01),np.ones(10))
  pl.figure(figsize=(8,5))
  pl.subplots_adjust(top=0.8,bottom=0.05,left=0.05,right=0.95)

  import glob
  files=glob.glob('./*.h')
  files.sort()
  maps=[gen_cmap(f) for f in files]

  l=len(maps)+1
  for i, m in enumerate(maps):
    pl.subplot(1,l,i+1)
    pl.axis("off")
    pl.imshow(a,aspect='auto',cmap=m,origin="lower")
    pl.text(0.5,1.01,m.name,rotation=50,fontsize=10,
            transform=pl.gca().transAxes,va='bottom')

  pl.savefig("colormaps.png",dpi=100,facecolor='gray')
  #pl.show()

if __name__=='__main__': show()
