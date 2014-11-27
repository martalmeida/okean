f='okean/version.py'
if 0:
  from okean.cookbook import run0
  v=run0('cd okean; git rev-list HEAD | wc -l; cd ..')
  v=int(v[0])
  l=open(f).readlines()
  l[1]='rev=%d\n'%v
  open(f,'w').writelines(l)
else: # just use current date
  import datetime
  open(f,'w').write('%s\n'%datetime.datetime.today().isoformat(' '))
