if 0:
  f='okean/version.py'
  from okean.cookbook import run0
  v=run0('cd okean; git rev-list HEAD | wc -l; cd ..')
  v=int(v[0])
  l=open(f).readlines()
  l[1]='rev=%d\n'%v
  open(f,'w').writelines(l)
else: # just use current date
  f='okean/__init__.py'
  import datetime
  s="__version__='%s'\n"% datetime.datetime.today().isoformat(' ')

  lines=open(f).readlines()
  for i,l in enumerate(lines):
    if l.strip().startswith('__version__='):
      lines[i]=s

  open(f,'w').writelines(lines)
