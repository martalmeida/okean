from okean.cookbook import run0

v=run0('cd okean; git rev-list HEAD | wc -l; cd ..')
v=int(v[0])

f='okean/version.py'
l=open(f).readlines()
print l[:2]
l[1]='rev=%d\n'%v
print l[:2]
open(f,'w').writelines(l)
