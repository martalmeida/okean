"""
Usefull stuff:
From http://aspn.activestate.com/ASPN/Cookbook/Python and other places

python 2x, 3x
"""

import os

##class _GetMaxIndex:
##  # sys.maxint may not be INT_MAX but LONG_MAX
##  # For 2.4.x, it is *just* that the documentation should be fixed; the
##  # code needs to stay as it is even though there is no convenient way
##  # to get at the maximum index. Just trying one time at startup would
##  # be the recommended way:
##  def __getslice__(self, a, mi):
##    #import sys
##    #sys.maxindex = mi
##    return mi
##
##def maxint():
##  return _GetMaxIndex()[:]

def maxint(warning=1):
  import sys
  try: return sys.maxint
  except:
    # note that the is no maxint in python 3. By maxsize can be used
    # as the max practical index to use
    if warning: print('warning: maxsize is not maxint!')
    return sys.maxsize


def relativepath(s1,s2):
  '''
  ex:
    >>relativepath('/home/user/','/home/')
    >>'../'
    >>relativepath('/home/user/','/home/user/path1/')
    >>'path1/'
  '''
  s1=s1.split(os.path.sep)
  s2=s2.split(os.path.sep)
  n=-1
  for i in range(min(len(s1),len(s2))):
    if s1[i]==s2[i]: n+=1

  s=''
  for i in range(len(s1)-n-2): s=os.path.join('..',s)
  s=[s[:-1]]
  s.extend(s2[n+1:])
  return os.path.join(*s)


def run0(s):
  '''
  exec full cmd as str, ex: run('ls -l')
  See also run
  '''
  out=os.popen(s)
  out=out.read()
  return out.split('\n')[:-1]


def run(*args):
  '''
  ex:
    >>run('ls','-l')
  '''
  import subprocess
  return subprocess.Popen(args,stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,close_fds=False).communicate()[0].split(b'\n')[:-1]
  # why close_fds=False --> http://bramp.net/blog/python-close_fds-issue
  # split(b is needed for python 3x


def parameters(only=None, exclude=None, ignore='self'):
    """Returns a dictionary of the calling functions 
       parameter names and values.

       The optional arguments can be used to filter the result:

           only           use this to only return parameters 
                          from this list of names.

           exclude        use this to return every parameter 
                          *except* those included in this list
                          of names.

           ignore         use this inside methods to ignore 
                          the calling object's name. For 
                          convenience, it ignores 'self' 
                          by default.

    """
    import inspect
    args, varargs, varkw, defaults = \
        inspect.getargvalues(inspect.stack()[1][0])
    if only is None:
        only = args[:]
        if varkw:
            only.extend(defaults[varkw].keys())
            defaults.update(defaults[varkw])
    if exclude is None:
        exclude = []
    exclude.append(ignore)
    return dict([(attrname, defaults[attrname])
        for attrname in only if attrname not in exclude])


def unique(seq, keepstr=True):
  '''
  ex:
    >>unique((1,2,3,4,4))
    >>(1, 2, 3, 4)
    >>unique('12344')
    >>'1234'
    >>unique('12344',0)
    >>['1', '2', '3', '4']
  '''
  t = type(seq)
  if t==str:
    t = (list, ''.join)[bool(keepstr)]
  seen = []
  return t(c for c in seq if not (c in seen or seen.append(c)))


def search(file,paths=False):
  """Searches a path for a specified file.
  Searches the paths specified in Environment for all files matching
  Filespec. If Environment is not specified, the system PATH is used.
  Original from Bill McNeill <billmcn@speakeasy.net>
  """

  import glob

  if not paths:
    paths = os.environ["PATH"]
  else:
    paths = os.environ[paths]

  for path in paths.split(os.path.pathsep):
    for match in glob.glob(os.path.join(path, file)):
      return match


def tree(dir='.', padding='', print_files=False,last=0,noext=()):
    '''Print tree structure of path specified.
    '''

    if last:
      print(padding[:-1] + ' -- ' + os.path.basename(os.path.abspath(dir)) + '/')
    else:
      print(padding[:-1] + '|-- ' + os.path.basename(os.path.abspath(dir)) + '/')
    padding = padding + '  '
    files = []
    if print_files:
        files = os.listdir(dir)
    else:
        files = [x for x in os.listdir(dir) if os.path.isdir(dir + os.sep + x)]

    files.sort()
    count = 0
    for file in files:
        count += 1
        path = dir + os.sep + file
        if os.path.isdir(path):
            if count == len(files):
                tree(path, padding + ' ', print_files,last=1)
            else:
                tree(path, padding + '|', print_files)
        else:
          ext=os.path.splitext(file)[1]
          if ext not in noext:
            if count == len(files):
              print(padding + ' -- ' + file)
            else:
              print(padding + '|-- ' + file)


def username(name=False):
  if name is False:
    name=os.environ['USER']
    if not name: name=run('whoami')[0]

  cmd="cat /etc/passwd | grep %s | awk -F : '{print $5}'"%name
  fname=run0(cmd)
  if fname: fname=fname[0]
  else: fname=name
  return name,fname


def machinename():
  if os.name is 'posix':
    res=run0('hostname --fqdn 2>/dev/null')
    if len(res): return res[0]
    else: return run0('hostname')[0]
  else:
    from socket import gethostname
    return gethostname()


def report_memory(i=''):
  pid = os.getpid()
  a2 = os.popen('ps -p %d -o rss,vsz,%%mem' % pid).readlines()
  print(i, ' ', a2[1],)
  return int(a2[1].split()[1])


def hsize(size):
  '''Bit sizes in human readable format
  '''

  if size<8:
    s='b'
    sz=size
  elif size<1024:
    s='B'
    sz=size/8.
  elif size<1024**2:
    s='Kb'
    sz=size/1024.
  elif size<1024**3:
    s='Mb'
    sz=size/(1024.**2)
  elif size<1024**4:
    s='Gb'
    sz=size/(1024.**3)
  else:
    s='Tb'
    sz=size/(1024.**4)

  return sz,s


def tar_gz_dir(source_dir, destination,infolder='auto'):
  '''
  source_dir: Source directory name.
  destination: Destination filename.
  (TAR-GZ-Archive *.tar.gz)

  default inFolder is basename(realpath(source_dir))
  for none use inFolder=''
  '''
  import tarfile

  t = tarfile.open(name = destination, mode = 'w:gz')
  if isinstance(source_dir,basestring):
    if infolder=='auto': infolder=os.path.basename(os.path.realpath(source_dir))
    t.add(source_dir, infolder)
  else:
    for sd in source_dir:
      in_folder=os.path.basename(os.path.realpath(sd))
      if infolder=='auto':  t.add(sd, in_folder)
      else: t.add(sd, infolder)

  t.close()


class odict(dict):
    '''
    Old, use instead (python >=2.7)
    from collections import OrderedDict
    '''

    def __init__(self, d={}):
        self._keys = d.keys()
        dict.__init__(self, d)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        dict.__setitem__(self, key, item)
        # a peculiar sharp edge from copy.deepcopy
        # we'll have our set item called without __init__
        if not hasattr(self, '_keys'):
            self._keys = [key,]
        if key not in self._keys:
            self._keys.append(key)

    def clear(self):
        dict.clear(self)
        self._keys = []

    def items(self):
        items = []
        for i in self._keys:
            items.append((i, self[i]))
        return items

    def keys(self):
        return self._keys

    def pop(self,key):
      self.__delitem__(key)

    def popitem(self):
        if len(self._keys) == 0:
            raise KeyError('dictionary is empty')
        else:
            key = self._keys[-1]
            val = self[key]
            del self[key]
            return key, val

    def setdefault(self, key, failobj = None):
        dict.setdefault(self, key, failobj)
        if key not in self._keys:
            self._keys.append(key)

    def update(self, d):
        for key in d.keys():
            if not self.has_key(key):
                self._keys.append(key)
        dict.update(self, d)

    def values(self):
        v = []
        for i in self._keys:
            v.append(self[i])
        return v

    def move(self, key, index):

        """ Move the specified to key to *before* the specified index. """

        try:
            cur = self._keys.index(key)
        except ValueError:
            raise KeyError(key)
        self._keys.insert(index, key)
        # this may have shifted the position of cur, if it is after index
        if cur >= index: cur = cur + 1
        del self._keys[cur]

    def index(self, key):
        if not self.has_key(key):
            raise KeyError(key)
        return self._keys.index(key)

    def __iter__(self):
        for k in self._keys:
            yield k

    def rename_key(self,key,newname):
      if not self.has_key(newname):
        self._keys[self._keys.index(key)]=newname
        val = self[key]
        dict.__delitem__(self, key)
        self[newname]=val
      else: raise KeyError('Key %s already exists' % newname)


def most_common(l):
  '''Most common element in a sequence
  ex:
    >>most_common([1,2,2,2,2,2,3,4])
    >>(5, 2)
  '''

  d = {}
  for elm in l:
    d[elm] = d.get(elm, 0) + 1

  counts = [(j,i) for i,j in d.items()]
  count, max_elm = max(counts)
  return count, max_elm

def isstr(s):
  '''
  String check compatible with python 2x and 3x
  Same as isinstance(s,basestring) in python 2x
  '''
  return isinstance(s,(''.__class__,u''.__class__))
  #or
  #try: basestring
  #except NameError: basestring=str
  #return isinstance(s,basestring)
