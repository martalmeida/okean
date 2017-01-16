import os
try: import memcache
except: memcache=False
import tempfile
import cPickle
import hashlib

class Cache:
  def __init__(self,dir='auto'):
    if dir is 'auto':
      self.dir=os.path.join(tempfile.gettempdir(),'pycache')
    else:
      self.dir=dir

    self.data={}

  def local_fname(self,label):
    label=hashlib.md5(label).hexdigest()
    f=os.path.join(self.dir,label)
    return f

  def is_stored(self,label,type):
    if type=='remote': 
      self.r=memcache.Client([kargs['client']])
      return self.r.has_key(label)
    elif type=='localmem':
      return self.data.has_key(label)
    elif type=='localfile':
      f=self.local_fname(label)
      return os.path.isfile(f)


  def store(self,label,data,type,*kargs):
    # remote, localmem, locafile
    if type=='remote':

      if not memcache:
        print 'error: cannot use memcache'
        return 1

      if kargs.has_key('time'): secs=kargs['time']
      else: secs=0

      self.r=memcache.Client([kargs['client']])
      s=self.r.set(label,data,time=secs)
      self.data[label]=None
      return s

    elif type=='localmem':
      self.data[label]=data
      return 0

    elif type=='localfile':
      f=self.local_fname(label)
      if not os.path.isdir(self.dir): os.makedirs(self.dir)
      s=cPickle.dump(data,open(f,'w+'))
      self.data[label]=None
      return s    
 
  def load(self,label,type):
    if type=='remote': 
      self.r=memcache.Client([kargs['client']])
      return self.r.get(label)
    elif type=='localmem':
      return self.data[label]
    elif type=='localfile':
      return cPickle.load(open(self.local_fname(label)))


