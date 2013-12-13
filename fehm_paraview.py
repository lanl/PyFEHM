
from optparse import OptionParser

usage = 'Usage: fehm_paraview.py [options] input grid'
parser = OptionParser(usage=usage)
parser.add_option('-i','--incon', dest = 'inconname', metavar = 'FILE', help = 'FEHM restart/initial conditions file')
parser.add_option('-c','--cont', dest = 'contname', metavar = 'FILE', help = 'FEHM contour output file, supports wildcard, e.g., *.csv')
parser.add_option('-e','--exe', dest = 'exename', metavar = 'PATH', help = 'path to paraview executable, default is paraview.exe')
parser.add_option('-d',action='store_true',dest='diff', default=False, help = 'PyFEHM will calculate differences in contour output')

(options, args) = parser.parse_args()

from fdata import*

if options.inconname != None:
	dat = fdata(args[0],args[1],options.inconname)
else:
	dat = fdata(args[0],args[1])

c = None
if options.contname: c = fcontour(options.contname)
exe = 'paraview.exe'
if options.exename: exe = options.exename
dat.paraview(exe=exe,contour=c,diff=options.diff)