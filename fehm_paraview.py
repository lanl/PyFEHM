
import os
from optparse import OptionParser

usage = 'Usage: fehm_paraview.py [options] input grid'
parser = OptionParser(usage=usage)
parser.add_option('-i','--incon', dest = 'inconname', metavar = 'FILE', help = 'FEHM restart/initial conditions file')
parser.add_option('-c','--cont', dest = 'contname', metavar = 'FILE', help = 'FEHM contour output file, supports wildcard, e.g., *.csv')
parser.add_option('-e','--exe', dest = 'exename', metavar = 'PATH', help = 'path to paraview executable, default is paraview.exe')
parser.add_option('-z','--zscale', dest = 'zscale', metavar = 'FL64', help = 'factor by which to scale view of z-axis')
parser.add_option('-d',action='store_true',dest='diff', default=False, help = 'PyFEHM will calculate differences in contour output')

(options, args) = parser.parse_args()

from fdata import*

if len(args) < 2:
	# look for fehmn.files
	if len(args) == 1: flname = args[0]
	elif os.path.isfile('fehmn.files'): flname = 'fehmn.files'
	
	dat = fdata(flname)
	# assess whether contour data should be loaded as well
elif len(args) == 1:
	# assume some form of fehmn.files - open it
else:
	# assume input and grid file have been provided
	if options.inconname != None:
		dat = fdata(args[0],args[1],options.inconname)
	else:
		dat = fdata(args[0],args[1])

c = None
if options.contname: c = fcontour(options.contname)
exe = 'paraview.exe'
if options.exename: exe = options.exename
zscale = 1.
if options.zscale: zscale = float(options.zscale)
dat.paraview(exe=exe,contour=c,diff=options.diff,zscale=zscale)