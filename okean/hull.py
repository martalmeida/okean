'''
Concave and convex hull
Several methods are available for convex hull. Some may be faster but may fail
(convex_hull_mpl and convex_hull_test). The one to try first is convex_hull, which
is based on Clarkson's hull code.
Concave hull/alpha shape is also based on Clarkson's hull code

mma jan 2012
m.martalmeida@gmail.com
'''

import pylab as pl
import numpy as np
import sys, os
import subprocess
import tempfile

hull_exec = "hull.exe"

# concave and convex hull based on Clarkson's hull code:

def convex_hull(x,y):
  '''Convex hull
     Returns indices i of the closed path x[i],y[i]
     Based on Clarkson's hull code
     (http://www.netlib.org/voronoi/hull.html)

     Ex:
     np.random.seed(0)
     x=np.random.rand(100)
     y=np.random.rand(100)

     pl.plot(x,y,'bo')
     i=convex_hull(x,y)
     pl.plot(x[i],y[i],'r-x')
  ''' 
  return _connect(_hull(x,y,'convex'))

def concave_hull(x,y):
  '''Concave hull (alpha shape)
     Returns indices i of the closed path x[i],y[i]
     Based on Clarkson's hull code
     (http://www.netlib.org/voronoi/hull.html)

     See convex_hull for additional help
  '''
  return _connect(_hull(x,y,'concave'))

def _connect(ind):
  '''Connects the indices returned bu _hull'''
  ind=ind[:]
  res=()
  k=False
  found=0
  while ind:
    if not found: # in case of more than one shape!!
      if k: res=res+(k,)
      k=ind[0]
      ind.pop(0)

    for i in range(len(ind)):
      found=0
      if ind[i][0]==k[-1]:
        k+=[ind[i][1]]
        ind.pop(i)
        found=1
        break
      elif ind[i][1]==k[-1]:
        k+=[ind[i][0]]
        ind.pop(i)
        found=1
        break

  if res: res=res+(k,) # last shape
  if not res: res=k # only one shape
  return res

def _hull(x,y,what,alpha='auto'):
    '''Runs Clarkson's hull
    Based on:
    http://stackoverflow.com/questions/6833243/how-can-i-find-the-alpha-shape-concave-hull-of-a-2d-point-cloud
    '''
    # Write points to tempfile
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    for point in np.vstack((x,y)).T.tolist():
        tmpfile.write("%0.7f %0.7f\n" % tuple(point))
    tmpfile.close()

    # Run hull
    if what=='concave':
      if alpha=='auto':
        command = "%s -A -m1000000 -oN < %s 2>/dev/null" % (hull_exec, tmpfile.name)
      else:
        command = "%s -aa%s -m1000000 -oN < %s 2>/dev/null" % (hull_exec, str(alpha),tmpfile.name)

    elif what=='convex':
      command = "%s -m1000000  < %s > hout-alf 2>/dev/null" % (hull_exec, tmpfile.name)
    #print >> sys.stderr, "Running command: %s" % command
    retcode = subprocess.call(command, shell=True)
    if retcode != 0:
        print >> sys.stderr, "Warning: bad retcode returned by hull.  Retcode value:" % retcode
    os.remove(tmpfile.name)

    # Parse results
    results_file = open("hout-alf")
    results_file.next() # skip header
    results_indices = [[int(i) for i in line.rstrip().split()] for line in results_file]
    #print "results length = %d" % len(results_indices)
    results_file.close()
    os.remove(results_file.name)

    return results_indices


# several convex hull methods: (try them if convex_hull fails!!)
# - using delaunay from mayplotlib
# - using another method, the quickhull
# - using the method from scipy

def convex_hull_mpl(x,y):
  '''Convex hull based on delaunay triangulation
     Based on http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg20889.html
     Returns indices i of the closed convex hull x[i],y[i]

     Based on matplotlib delaunay triangulations, which may fail!! If so, use convex_hull or
     convex_hull_scipy

     Ex:
     np.random.seed(0)
     x = np.random.rand(10)
     y = np.random.rand(10)
     pl.figure()
     pl.plot(x, y, 'bo')
     b=convex_hull_i(x,y)
     pl.plot(x[b],y[b])
  '''

  if 0:
    triang = pl.matplotlib.tri.Triangulation(x, y)

    # Determine which triangulation edges are boundary edges; these
    # correspond to edges that have no neighboring triangle.  Store them as
    # pairs of (start,end) point indices.
    boundaries = []
    for i in range(len(triang.triangles)):
      for j in range(3):
        if triang.neighbors[i,j] < 0:
          # Triangle edge (i,j) has no neighbor so is a boundary edge.
          boundaries.append((triang.triangles[i,j],
                           triang.triangles[i,(j+1)%3]))


    # make it continuous:
    b=list(boundaries[0])
    boundaries.pop(0)
    while boundaries:
      for i in range(len(boundaries)):
        k=boundaries[i]
        if k[0]==b[-1]:
          b+=[k[1]]
          boundaries.pop(i)
          break

  
  else: # or just:
      triang = pl.matplotlib.tri.delaunay.Triangulation(x, y)
      b=triang.hull+[triang.hull[0]]
 
  return np.asarray(b)


def convex_hull_test(x,y):
  '''Convex hull... another method
     Based on: http://en.literateprograms.org/Quickhull_%28Python,_arrays%29
     Failed in some tests!!
     Returns the convex hull (x,y), not the indices

     Ex:
     np.random.seed(0)
     x = np.random.rand(10)
     y = np.random.rand(10)
     pl.figure()
     pl.plot(x, y, 'bo')
     xc,yc=convex_hull_test(x,y)
     pl.plot(xc,yc)
  '''
  sample=np.vstack((x,y)).T
  link = lambda a,b: np.concatenate((a,b[1:]))
  edge = lambda a,b: np.concatenate(([a],[b]))

  def dome(sample,base):
    h, t = base
    dists = np.dot(sample-h, np.dot(((0,-1),(1,0)),(t-h)))
    outer = np.repeat(sample, dists>0, 0)

    if len(outer):
      pivot = sample[np.argmax(dists)]
      return link(dome(outer, edge(h, pivot)),dome(outer, edge(pivot, t)))
    else:
      return base

  if len(sample) > 2:
    axis = sample[:,0]
    base = np.take(sample, [np.argmin(axis), np.argmax(axis)], 0)
    res=link(dome(sample, base),dome(sample, base[::-1]))
    return res[:,0],res[:,1]
  else:
    return x,y


def convex_hull_scipy(x,y):
  '''
     Calculate subset of points that make a convex hull around points
     Recursively eliminates points that lie inside two neighbouring points until only convex hull is remaining.

     based on: http://www.scipy.org/Cookbook/Finding_Convex_Hull

     Returns the convex hull (x,y), not the indices
     For the indices, convex_hull_mpl may be used (may be much faster, but may fail)

     Ex:
     np.random.seed(0)
     x = np.random.rand(10)
     y = np.random.rand(10)
     pl.figure()
     pl.plot(x, y, 'bo')
     xc,yc=convex_hull_scipy(x,y)
     pl.plot(xc,yc)
  '''
  points=np.vstack((x,y))

  graphic=False #  use pylab to show progress?
  smidgen=0.0075 # offset for graphic number labels - useful values depend on your data range

  def _angle_to_point(point, centre):
    '''calculate angle in 2-D between points and x axis'''
    delta = point - centre
    res = np.arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += np.pi
    return res


  def _draw_triangle(p1, p2, p3, **kwargs):
    tmp = np.vstack((p1,p2,p3))
    x,y = [x[0] for x in zip(tmp.transpose())]
    p.fill(x,y, **kwargs)
    #time.sleep(0.2)


  def area_of_triangle(p1, p2, p3):
    '''calculate area of any triangle given coordinates of the corners'''
    return np.linalg.norm(np.cross((p2 - p1), (p3 - p1)))/2.

  if graphic:
    p.clf()
    p.plot(points[0], points[1], 'ro')

  n_pts = points.shape[1]
  assert(n_pts > 5)
  centre = points.mean(1)
  if graphic: p.plot((centre[0],),(centre[1],),'bo')
  angles = np.apply_along_axis(_angle_to_point, 0, points, centre)
  pts_ord = points[:,angles.argsort()]
  if graphic:
    for i in xrange(n_pts): 
      p.text(pts_ord[0,i] + smidgen, pts_ord[1,i] + smidgen, '%d' % i)

  pts = [x[0] for x in zip(pts_ord.transpose())]
  prev_pts = len(pts) + 1
  k = 0
  while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        if graphic: p.gca().patches = []
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                   pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if graphic:
                _draw_triangle(centre, pts[i], pts[(i + 1) % n_pts], \
                               facecolor='blue', alpha = 0.2)
                _draw_triangle(centre, pts[(i + 1) % n_pts], \
                               pts[(i + 2) % n_pts], \
                               facecolor='green', alpha = 0.2)
                _draw_triangle(centre, pts[i], pts[(i + 2) % n_pts], \
                               facecolor='red', alpha = 0.2)
            if Aij + Ajk < Aik:
                if graphic: p.plot((pts[i + 1][0],),(pts[i + 1][1],),'go')
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1

  pts+=[pts[0]]
  pts=np.asarray(pts)
  return pts[:,0],pts[:,1]

def simplify(x,y,lev=0):
  '''ugly way to reduce the number of x,y points'''
  from okean import ticks
  nnx=ticks.nicenum(x.mean(),True)
  nny=ticks.nicenum(y.mean(),True)

  expx=np.floor(np.log10(nnx))
  expy=np.floor(np.log10(nny))

  x=np.round(x,lev+int(expx))
  y=np.round(y,lev+int(expy))

  xy=x+1j*y
  np.random.shuffle(xy) # assure i is always different!
  xy,ind=np.unique(xy,True)

  return xy.real,xy.imag,ind
  



