'''
Ticks and ticklabels for graphic axes
from:

GRAPHIC GEMS I, Andrew S. Glassner, 1990
Nice numbers for graph labels, p 61

Online book:
http://www.scribd.com/doc/7075218/Graphics-Gems-I

Repository:
http://tog.acm.org/resources/GraphicsGems/

Author homepage:
http://www.glassner.com


Martinho MA, oct2009
mma@ua.pt
'''

import numpy as np

def nicenum(x,rnd):
  '''
  Find a "nice" number approximately equal to x.
  Round the number if rnf is True, take ceiling if rnd is False.
  '''

  exp=np.floor(np.log10(x))
  f=x/10**exp
  if rnd:
    if f<1.5: nf=1.
    elif f<3.: nf=2.
    elif f<7.: nf=5.
    else:     nf=10.

  else:
    if f<=1.:  nf=1.
    elif f<=2.: nf=2.
    elif f<=5.: nf=5.
    else:       nf=10.

  return nf*10**exp


def loose_label(min,max,ntick=5,labels=False):
  '''
  Label the data range from min to max loosely
  (tight method is similar)

  Inputs:

  min, max: data range

  ntick: desired number of tick marks (may not be real number of
         ticks returned

  labels: tick labels with the "best" number of fractional digits

  Example:

  >>loose_label(105,543)
  array([ 100.,  200.,  300.,  400.,  500.,  600.])
  >>
  >>loose_label(2.03,2.17,labels=True)
  (array([ 2.  ,  2.05,  2.1 ,  2.15,  2.2 ]),
  ['2.00', '2.05', '2.10', '2.15', '2.20'])
  '''

  range=nicenum(max-min,False)
  d=nicenum(range/(ntick-1),True)
  graphmin=np.floor(min/d)*d
  graphmax=np.ceil(max/d)*d
  nfrac=int(np.max((-np.floor(np.log10(d)),0))) # number of fractional digits to show

  if not labels:
    return np.arange(graphmin,graphmax+d/2.,d)
  else:
    format='%.'+str(nfrac)+'f'
    ret=np.arange(graphmin,graphmax+d/2.,d)
    return ret,[format % i for i in ret]


def tight(min,max,ntick=5,labels=False):
  '''
   Returns output of loose_label inside [min,max]
   See loose_label
  '''

  if labels:
    tk,tks=loose_label(min,max,ntick,labels)
  else:
    tk=loose_label(min,max,ntick,labels)

  i=np.where(tk>min)[0]
  j=np.where(tk<max)[0]

  if len(i): tk=tk[i[0]:]
  if len(j): tk=tk[:j[-1]]

  if labels:
    if len(i): tks=tks[i[0]:]
    if len(j): tks=tks[:j[-1]]


  if labels: return tk,tks
  else: return tk


def loose_label_n(min,max,ntick=5,labels=False):
  '''
  loose_label with at least n ticks
  '''

  tk=[]
  n0=ntick
  while len(tk)< ntick:
    tk=loose_label(min,max,ntick=n0)
    n0+=1

  return tk

def nice_bounds(arr):
  '''nice bounds for n-d array values'''
  x=loose_label(arr.min(),arr.max())
  return x[0],x[-1]
