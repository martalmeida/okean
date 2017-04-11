# f2py compilation, ex:
# f2py -c some.f -m some --fcompiler=gnu95
# f2py -c some.f -m some --f90flags==-32
# f2py -c some.f -m some --fcompiler=pg

FCname=pg
FCname=gnu95
#FCname=intelem

# check with: f2py -c --help-fcompiler

all: alga rtools

check:
	f2py -c --help-fcompiler

defaul:
	f2py -c alga.f -m alga
	f2py -c rtools.f -m rtools

alga:
	f2py -c alga_.f -m alga_ --fcompiler=$(FCname)

rtools:
	f2py -c rtools.f -m rtools --fcompiler=$(FCname)

rtools_test:
	f2py -c rtools_test.f -m rtools_test --fcompiler=$(FCname)
