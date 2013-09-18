"""For reading FEHM output files."""

import numpy as np
import os
try:
	from matplotlib import pyplot as plt
	from mpl_toolkits.mplot3d import axes3d
	from matplotlib import cm
	import matplotlib
except ImportError:
	'placeholder'
from copy import copy
from ftool import*
import platform

from fdflt import*
dflt = fdflt()

WINDOWS = platform.system()=='Windows'
if WINDOWS: copyStr = 'copy'; delStr = 'del'
else: copyStr = 'cp'; delStr = 'rm'

if True: 					# output variable dictionaries defined in here, indented for code collapse
	cont_var_names_avs=dict([
	('X coordinate (m)','x'),
	('Y coordinate (m)','y'),
	('Z coordinate (m)','z'),
	('node','n'),
	('Liquid Pressure (MPa)','P'),
	('Vapor Pressure (MPa)','P_vap'),
	('Capillary Pressure (MPa)','P_cap'),
	('Saturation','saturation'),
	('Temperature (deg C)','T'),
	('Porosity','por'),
	('X Permeability (log m**2)','perm_x'),
	('Y Permeability (log m**2)','perm_y'),
	('Z Permeability (log m**2)','perm_z'),
	('X displacement (m)','disp_x'),
	('Y displacement (m)','disp_y'),
	('Z displacement (m)','disp_z'),
	('X stress (MPa)','strs_xx'),
	('Y stress (MPa)','strs_yy'),
	('Z stress (MPa)','strs_zz'),
	('XY stress (MPa)','strs_xy'),
	('XZ stress (MPa)','strs_xz'),
	('YZ stress (MPa)','strs_yz'),
	('Youngs Mod (MPa)','E'),
	('Excess Shear (MPa)','tau_ex'),
	('Shear Angle (deg)','phi_dil'),
	('Zone','zone'),
	('Liquid Density (kg/m**3)','density'),
	('Vapor Density (kg/m**3)','density_vap'),
	('Source (kg/s)','flow'),
	('Liquid Flux (kg/s)','flux'),
	('Vapor Flux (kg/s)','flux_vap'),
	('Volume Strain','strain'),
	('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
	('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
	('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
	('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
	('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
	('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
	]) 

	cont_var_names_tec=dict([
	('X coordinate (m)','x'),
	('Y coordinate (m)','y'),
	('Z coordinate (m)','z'),
	('X Coordinate (m)','x'),
	('Y Coordinate (m)','y'),
	('Z Coordinate (m)','z'),
	('node','n'),
	('Node','n'),
	('Liquid Pressure (MPa)','P'),
	('Vapor Pressure (MPa)','P_vap'),
	('Capillary Pressure (MPa)','P_cap'),
	('Saturation','saturation'),
	('Water Saturation','saturation'),
	('Super-Critical/Liquid CO2 Saturation','co2_liquid'),
	('Gaseous CO2 Saturation','co2_gas'),
	('Dissolved CO2 Mass Fraction','co2_aq'),
	('CO2 Phase State','co2_phase'),
	('Temperature (<sup>o</sup>C)','T'),
	('Temperature (deg C)','T'),
	('Porosity','por'),
	('X Permeability (log m**2)','perm_x'),
	('Y Permeability (log m**2)','perm_y'),
	('Z Permeability (log m**2)','perm_z'),
	('X displacement (m)','disp_x'),
	('Y displacement (m)','disp_y'),
	('Z displacement (m)','disp_z'),
	('X stress (MPa)','strs_xx'),
	('Y stress (MPa)','strs_yy'),
	('Z stress (MPa)','strs_zz'),
	('XY stress (MPa)','strs_xy'),
	('XZ stress (MPa)','strs_xz'),
	('YZ stress (MPa)','strs_yz'),
	('Youngs Mod (MPa)','E'),
	('Excess Shear (MPa)','tau_ex'),
	('Shear Angle (deg)','phi_dil'),
	('Zone','zone'),
	('Liquid Density (kg/m**3)','density'),
	('Vapor Density (kg/m**3)','density_vap'),
	('Source (kg/s)','flow'),
	('Liquid Flux (kg/s)','flux'),
	('Vapor Flux (kg/s)','flux_vap'),
	('Volume Strain','strain'),
	('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
	('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
	('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
	('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
	('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
	('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
	])

	cont_var_names_surf=dict([
	('X coordinate (m)','x'), 
	('X (m)','x'), 
	('Y coordinate (m)','y'), 
	('Y (m)','y'), 
	('Z coordinate (m)','z'), 
	('Z (m)','z'),
	('node','n'),
	('Liquid Pressure (MPa)','P'),
	('Vapor Pressure (MPa)','P_vap'),
	('Capillary Pressure (MPa)','P_cap'),
	('Saturation','saturation'),
	('Temperature (deg C)','T'),
	('Porosity','por'),
	('X Permeability (log m**2)','perm_x'),
	('Y Permeability (log m**2)','perm_y'),
	('Z Permeability (log m**2)','perm_z'),
	('X displacement (m)','disp_x'),
	('Y displacement (m)','disp_y'),
	('Z displacement (m)','disp_z'),
	('X stress (MPa)','strs_xx'),
	('Y stress (MPa)','strs_yy'),
	('Z stress (MPa)','strs_zz'),
	('XY stress (MPa)','strs_xy'),
	('XZ stress (MPa)','strs_xz'),
	('YZ stress (MPa)','strs_yz'),
	('Youngs Mod (MPa)','E'),
	('Excess Shear (MPa)','tau_ex'),
	('Shear Angle (deg)','phi_dil'),
	('Zone','zone'),
	('Liquid Density (kg/m**3)','density'),
	('Vapor Density (kg/m**3)','density_vap'),
	('Source (kg/s)','flow'),
	('Liquid Flux (kg/s)','flux'),
	('Vapor Flux (kg/s)','flux_vap'),
	('Volume Strain','strain'),
	('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
	('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
	('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
	('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
	('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
	('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
	('Water Saturation','saturation'),
	('Super-Critical/Liquid CO2 Saturation','co2_liquid'),
	('Gaseous CO2 Saturation','co2_gas'),
	('Dissolved CO2 Mass Fraction','co2_aq'),
	('CO2 Phase State','co2_phase'),
	('Aqueous_Species_001','species001_aq'),
	('Aqueous_Species_002','species002_aq'),
	('Aqueous_Species_003','species003_aq'),
	('Aqueous_Species_004','species004_aq'),
	('Aqueous_Species_005','species005_aq'),
	('Aqueous_Species_006','species006_aq'),
	('Aqueous_Species_007','species007_aq'),
	('Aqueous_Species_008','species008_aq'),
	('Aqueous_Species_009','species009_aq'),
	('Aqueous_Species_010','species010_aq'),
	('Aqueous_Species_011','species011_aq'),
	('Aqueous_Species_012','species012_aq'),
	('Aqueous_Species_013','species013_aq'),
	('Aqueous_Species_014','species014_aq'),
	('Aqueous_Species_015','species015_aq'),
	('Aqueous_Species_016','species016_aq'),
	('Aqueous_Species_017','species017_aq'),
	('Aqueous_Species_018','species018_aq'),
	('Aqueous_Species_019','species019_aq'),
	('Aqueous_Species_020','species020_aq'),
	])

	hist_var_names=dict([
	('denAIR','density_air'),
	('disx','disp_x'),
	('disy','disp_y'),
	('disz','disp_z'),
	('enth','enthalpy'),
	('glob','global'),
	('humd','humidity'),
	('satr','saturation'),
	('strain','strain'),
	('strx','strs_xx'),
	('stry','strs_yy'),
	('strz','strs_zz'),
	('strxy','strs_xy'),
	('strxz','strs_xz'),
	('stryz','strs_yz'),
	('wcon','water_content'),
	('denWAT','density'),
	('flow','flow'),
	('visAIR','viscosity_air'),
	('visWAT','viscosity'),
	('wt','water_table'),
	('presCAP','P_cap'),
	('presVAP','P_vap'),
	('presWAT','P'),
	('presCO2','P_co2'),
	('temp','T'),
	('co2md','massfrac_co2_aq'),
	('co2mf','massfrac_co2_free'),
	('co2mt','mass_co2'),
	('co2sg','saturation_co2g'),
	('co2sl','saturation_co2l'),
	])
	
	flxz_water_names = [
	'water_source',
	'water_sink',
	'water_net',
	'water_boundary',]
	
	flxz_vapor_names = [
	'vapor_source',
	'vapor_sink',
	'vapor_net',
	'vapor_boundary',]
	
	flxz_co2_names = [
	'co2_source',
	'co2_sink',
	'co2_in',
	'co2_out',
	'co2_boundary',
	'co2_sourceG',
	'co2_sinkG',
	'co2_inG',
	'co2_outG']
class fcontour(object): 					# Reading and plotting methods associated with contour output data.
	'''Contour output information object.
	
	'''
	def __init__(self,filename=None,latest=False,first=False,nearest=None):
		self._filename=filename
		self._times=[]   
		self._format = ''
		self._data={}
		self._row=None
		self._variables=[]  
		self.key_name=[]
		self._keyrows={}
		self.column_name=[]
		self.num_columns=0
		self._latest = latest
		self._first = first
		self._nearest = nearest
		self._nkeys=1
		if self._filename: self.read(filename,self._latest,self._first,self._nearest)
	def __getitem__(self,key):
		if key in self.times:
			return self._data[key]
		else: return None
	def read(self,filename,latest=False,first=False,nearest=[]): 						# read contents of file
		'''Read in FEHM contour output information.
		
		:param filename: File name for output data, can include wildcards to define multiple output files.
		:type filename: str
		:param latest: 	Boolean indicating PyFEHM should read the latest entry in a wildcard search.
		:type latest: bool
		:param first: Boolean indicating PyFEHM should read the first entry in a wildcard search.
		:type first: bool
		:param nearest: Read in the file with date closest to the day supplied. List input will parse multiple output files.
		:type nearest: fl64,list
		'''
		from glob import glob
		if isinstance(filename,list):
			files = filename
		else:
			files=glob(filename)
			if len(files)==0: print 'ERROR: '+filename+' not found'; return
			if self._nearest or latest or first:
				files = filter(os.path.isfile, glob(filename))
				files.sort(key=lambda x: os.path.getmtime(x))
				files2 = []
				if first:
					files2.append(files[0])
				if self._nearest:
					ts = [fl.split('_node')[0] for fl in files]
					ts = [fl.split('_sca')[0] for fl in ts]
					ts = [fl.split('_con')[0] for fl in ts]
					ts = [fl.split('_days')[0] for fl in ts]
					ts = [fl.split('.')[-2:] for fl in ts]
					ts = [float(tsi[0]+'.'+tsi[1]) for tsi in ts]
					if isinstance(self._nearest,(float,int)):
						tsi = min(enumerate(ts), key=lambda x: abs(x[1]-self._nearest))[0]
						files2.append(files[tsi])
					elif isinstance(self._nearest,(list,tuple)):
						for near in self._nearest:
							tsi = min(enumerate(ts), key=lambda x: abs(x[1]-near))[0]
							files2.append(files[tsi])
				if latest:
					files2.append(files[-1])
				files = []
				for file in files2:
					if file not in files: files.append(file)
		configured=False
		for i,fname in enumerate(files):
			print fname
			self._file=open(fname,'rU')
			if not configured:
				header=self._file.readline()
				self._detect_format(header)
				if self._format=='tec':
					self._setup_headers_tec(header)
				elif self._format=='avs':
					self._setup_headers_avs(header)
				elif self._format=='avsx':
					self._setup_headers_avsx(header)
				elif self._format=='surf':
					self._setup_headers_surf(header)
				else:
					print 'Unrecognised format'; return
				self.num_columns = len(self.variables)+1
			if self.format == 'tec': self._read_data_tec()
			elif self.format == 'surf': self._read_data_surf(fname,configured)
			elif self.format == 'avs': self._read_data_avs(fname,configured)
			elif self.format == 'avsx': 
				if configured: self._read_data_avsx()
				else: self._read_data_avsx(header)
			self._file.close()
			configured = True
		if dflt.parental_cont:
			print ''
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			print 'WARNING:'
			print ''
			print 'Contour data is indexed using the Pythonic convention in which the first index is 0. FEHM node numbering convention begins at 1.'
			print ''
			print 'THEREFORE, to get the correct contour value for a particular node, you need to pass the node index MINUS 1. Using node index to access contour data will return incorrect values.'
			print ''
			print 'For example:'
			print '>>> node10 = dat.grid.node[10]'
			print '>>> c = fcontour(\'*.csv\')'
			print '>>> T_node10 = c[c.times[-1]][\'T\'][node10.index - 1]'
			print '  or'
			print '>>> T_node10 = c[c.times[-1]][\'T\'][9]'
			print 'will return the correct value for node 10.'
			print ''
			print 'Do not turn off this message unless you understand how to correctly access nodal values from contour data.' 
			print 'To turn off this message, open the environment file \'fdflt.py\' and set self.parental_cont = False'
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			print ''
	def _detect_format(self,header):
		if header.startswith('TITLE ='):		# check for TEC output
			self._format = 'tec'
		elif header.startswith('node, '):		# check for SURF output
			self._format = 'surf'
		elif header.startswith('nodes at '):	# check for AVSX output
			self._format = 'avsx'
		elif header.split()[0].isdigit():			# check for AVS output
			self._format = 'avs'
	def _setup_headers_avsx(self,header): 		# headers for the AVSX output format
		header = header.strip().split(' : ')
		self._variables.append('n')
		for key in header[1:]: 
			if key in cont_var_names_avs.keys():
				var = cont_var_names_avs[key]
			else: var = key
			self._variables.append(var)
	def _read_data_avsx(self,header=''):				# read data in AVSX format
		if not header: header = self._file.readline()
		header = header.split('nodes at ')[1]
		header = header.split('days')[0]
		time = float(header)*24*2600
		self._times.append(time)		
		lns = self._file.readlines()
		data = []
		for ln in lns: data.append([float(d) for d in ln.strip().split(':')])
		data = np.array(data)
		self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
	def _setup_headers_avs(self,header): 		# headers for the AVS output format
		lns_num = int(header.strip().split(' ')[0])
		self._variables.append('n')
		for i in range(lns_num): 
			ln = self._file.readline().strip().split(',')[0]
			var = cont_var_names_avs[ln]
			self._variables.append(var)
	def _read_data_avs(self,fname,configured):		# read data in AVS format
		lni = fname.split('.',1)[1]
		lni = lni.split('_',1)
		time = float(lni[0])
		if lni[1].startswith('days'):
			time = time*24*2600
		self._times.append(time)		
		if configured: 
			for var in self.variables: lni=self._file.readline()
		lns = self._file.readlines()
		data = []
		for ln in lns: data.append([float(d) for d in ln.strip().split()])
		data = np.array(data)
		self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
	def _setup_headers_surf(self,header): 		# headers for the SURF output format
		header = header.strip().split(', ')
		for key in header: 
			varname = key.split('"')[0]
			if varname in cont_var_names_surf.keys():
				var = cont_var_names_surf[varname]
			else: var = varname
			self._variables.append(var)
	def _read_data_surf(self,fname,configured):		# read data in SURF format
		lni = fname.split('.',1)[1]
		if '_sca_node' in fname: 
			lni = fname.split('_sca_node',1)[0]
		elif '_con_node' in fname: 
			lni = fname.split('_con_node',1)[0] 
		day_flag = False
		if lni.endswith('days'): day_flag = True; lni = lni[:-5]
		lni = lni.split('.')
		time = float(lni[-2]+'.'+lni[-1])
		self._times.append(time)		
		if configured: lni=self._file.readline()
		lns = self._file.readlines()
		data = []
		for ln in lns: data.append([float(d) for d in ln.strip().split(',')])
		data = np.array(data)
		self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
	def _setup_headers_tec(self,header): 		# headers for the TEC output format
		while not header[0].startswith('VARIABLE'):
			header = self._file.readline()
			header = header.split(' "')
		for key in header[1:]: 
			varname = key.split('"')[0]
			if varname in cont_var_names_surf.keys():
				var = cont_var_names_tec[varname]
			else: var = varname
			self._variables.append(var)
	def _read_data_tec(self):						# read data in TEC format
		ln = self._file.readline()
		while not ln.startswith('ZONE'):
			ln = self._file.readline()
		lni = ln.split('"')[1]
		time = lni.split('days')[0].strip()
		time = float(time.split()[-1].strip())
#		unit = lni.split()[3]
#		if unit == 'days':
#			time = time*24*2600
		self._times.append(time)
		lni = ln.split(',')[1]
		nds = int(lni.split('=')[1])
		lns = self._file.readlines()
		data = []
		for ln in lns[:nds]: data.append([float(d) for d in ln.strip().split()])
		data = np.array(data)
		if len(data[0])< len(self.variables): 		# insert xyz data from previous read
			x0 = self._data[self.times[0]]['x']
			y0 = self._data[self.times[0]]['y']
			z0 = self._data[self.times[0]]['z']
			j = 0
			data2 = []
			for var in self.variables:
				if var == 'x': data2.append(x0)
				elif var == 'y': data2.append(y0)
				elif var == 'z': data2.append(z0)
				else: data2.append(data[:,j]); j +=1
			data = np.transpose(np.array(data2))
		
		self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
	def _check_inputs(self,variable, time, slice):	# assesses whether sufficient input information for slice plot
		if not variable: 
			print 'Error: no plot variable specified.'
			print 'Options are'
			for var in self.variables: print var
			return True
		if time==None: 
			print 'Error: no plot time specified.'
			print 'Options are'
			for time in self.times: print time
			return True
		if not slice: 
			print 'Error: slice orientation undefined.'
			print 'Options are'
			print '[\'x\',float] - slice parallel to y-axis at x=float'
			print '[\'y\',float] - slice parallel to x-axis at y=float'
			print '[\'theta\',float] - angle measured anti-clockwise from +x'
			print '[[float,float],[float,float]] - point to point'
			return True
		return False
	def slice(self, variable, slice, divisions, time=None, method='nearest'):
		'''Returns mesh data for a specified slice orientation from 3-D contour output data.
		
		:param variable: Output data variable, for example 'P' = pressure. Alternatively, variable can be a five element list, first element 'cfs', remaining elements fault azimuth (relative to x), dip, friction coefficient and cohesion. Will return coulomb failure stress.
		:type variable: str
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data.
		:type time: fl64
		:param slice: List specifying orientation of output slice, e.g., ['x',200.] is a vertical slice at ``x = 200``, ['z',-500.] is a horizontal slice at ``z = -500.``, [point1, point2] is a fixed limit vertical or horizontal domain corresponding to the bounding box defined by point1 and point2.
		:type slice: lst[str,fl64]
		:param divisions: Resolution to supply mesh data.
		:type divisions: [int,int]
		:param method: Method of interpolation, options are 'nearest', 'linear'.
		:type method: str:	
		:returns: X -- x-coordinates of mesh data.
		
		'''
		if time==None: 
			if np.min(self.times)<0: time = self.times[0]
			else: time = self.times[-1]
		from scipy.interpolate import griddata		
		delta = False
		if isinstance(time,list) or isinstance(time,np.ndarray):
			if len(time)>1: 
				time0 = np.min(time)
				time = np.max(time)
				delta=True		
		dat = self[time]
		
		# check to see if cfs plot requested
		cfs = False
		if isinstance(variable,list):
			if variable[0] in ['cfs','CFS']: cfs = True
		
		if not cfs:
			if delta: dat0 = self[time0]
			if isinstance(slice[0],str):
				if slice[0].startswith('y'):
					xmin = np.min(dat['x']);xmax = np.max(dat['x'])
					ymin = np.min(dat['z']);ymax = np.max(dat['z'])		
					if slice[1] == None:
						points = np.transpose(np.array([dat['x'],dat['z'],np.ones((1,len(dat['z'])))[0]]))
						slice[1] = 1
					else:
						points = np.transpose(np.array([dat['x'],dat['z'],dat['y']]))
				elif slice[0].startswith('x'):
					xmin = np.min(dat['y']);xmax = np.max(dat['y'])
					ymin = np.min(dat['z']);ymax = np.max(dat['z'])		
					if slice[1] == None:
						points = np.transpose(np.array([dat['y'],dat['z'],np.ones((1,len(dat['z'])))[0]]))
						slice[1] = 1
					else:
						points = np.transpose(np.array([dat['y'],dat['z'],dat['x']]))
				elif slice[0].startswith('z'):
					xmin = np.min(dat['x']);xmax = np.max(dat['x'])
					ymin = np.min(dat['y']);ymax = np.max(dat['y'])		
					if slice[1] == None:
						points = np.transpose(np.array([dat['x'],dat['y'],np.ones((1,len(dat['y'])))[0]]))
						slice[1] = 1
					else:
						points = np.transpose(np.array([dat['x'],dat['y'],dat['z']]))
				elif slice[0].startswith('theta'): 
					print 'Not ready yet'; return
				xrange = np.linspace(xmin,xmax,divisions[0])
				yrange = np.linspace(ymin,ymax,divisions[1])
				X,Y = np.meshgrid(xrange,yrange)
				Z = (X+np.sqrt(1.757))/(X+np.sqrt(1.757))*slice[1]
				pointsI = np.transpose(np.reshape((X,Y,Z),(3,X.size)))
				vals = np.transpose(np.array(dat[variable]))
				valsI = griddata(points,vals,pointsI,method=method)
				valsI =  np.reshape(valsI,(X.shape[0],X.shape[1]))
				if delta:
					vals = np.transpose(np.array(dat0[variable]))
					valsI0 = griddata(points,vals,pointsI,method=method)
					valsI0 =  np.reshape(valsI0,(X.shape[0],X.shape[1]))
					valsI = valsI - valsI0
			elif isinstance(slice[0],list):
				# check if horizontal or vertical slice
				dx,dy,dz = abs(slice[0][0]-slice[1][0]),abs(slice[0][1]-slice[1][1]),abs(slice[0][2]-slice[1][2])
				if 100*dz<dx and 100*dz<dy: 	#horizontal
					xmin,xmax = np.min([slice[0][0],slice[1][0]]),np.max([slice[0][0],slice[1][0]])
					ymin,ymax = np.min([slice[0][1],slice[1][1]]),np.max([slice[0][1],slice[1][1]])
					xrange = np.linspace(xmin,xmax,divisions[0])
					yrange = np.linspace(ymin,ymax,divisions[1])
					X,Y = np.meshgrid(xrange,yrange)
					Z = (X+np.sqrt(1.757))/(X+np.sqrt(1.757))*(slice[0][2]+slice[1][2])/2
				else: 							#vertical 
					xmin,xmax = 0,np.sqrt((slice[0][0]-slice[1][0])**2+(slice[0][1]-slice[1][1])**2)
					ymin,ymax = np.min([slice[0][2],slice[1][2]]),np.max([slice[0][2],slice[1][2]])
					xrange = np.linspace(xmin,xmax,divisions[0])
					yrange = np.linspace(ymin,ymax,divisions[1])
					X,Z = np.meshgrid(xrange,yrange)
					Y = X/xmax*abs(slice[0][1]-slice[1][1]) + slice[0][1]
					X = X/xmax*abs(slice[0][0]-slice[1][0]) + slice[0][0]
				points = np.transpose(np.array([dat['x'],dat['y'],dat['z']]))
				pointsI = np.transpose(np.reshape((X,Y,Z),(3,X.size)))
				vals = np.transpose(np.array(dat[variable]))
				valsI = griddata(points,vals,pointsI,method=method)
				valsI =  np.reshape(valsI,(X.shape[0],X.shape[1]))
				if delta:
					vals = np.transpose(np.array(dat0[variable]))
					valsI0 = griddata(points,vals,pointsI,method=method)
					valsI0 =  np.reshape(valsI0,(X.shape[0],X.shape[1]))
					valsI = valsI - valsI0
			
		else:
			if delta: time0 = time[0]; time = time[-1]
			X,Y,Z,sxx = self.slice('strs_xx', slice, divisions, time, method)
			X,Y,Z,syy = self.slice('strs_yy', slice, divisions, time, method)
			X,Y,Z,szz = self.slice('strs_zz', slice, divisions, time, method)
			X,Y,Z,sxy = self.slice('strs_xy', slice, divisions, time, method)
			X,Y,Z,sxz = self.slice('strs_xz', slice, divisions, time, method)
			X,Y,Z,syz = self.slice('strs_yz', slice, divisions, time, method)
			X,Y,Z,sp  = self.slice('P',       slice, divisions, time, method)
			
			dip = variable[2]/180.*math.pi
			azi = variable[1]/180.*math.pi+3.14159/2.
			nhat = np.array([np.cos(azi)*np.sin(dip),np.sin(azi)*np.sin(dip),np.cos(dip)])
			mu = variable[3]
			cohesion = variable[4]
			
			px = sxx*nhat[0]+sxy*nhat[1]+sxz*nhat[2]
			py = sxy*nhat[0]+syy*nhat[1]+syz*nhat[2]
			pz = sxz*nhat[0]+syz*nhat[1]+szz*nhat[2]
			
			sig = px*nhat[0]+py*nhat[1]+pz*nhat[2]
			tau = np.sqrt(px**2+py**2+pz**2 - sig**2)
			valsI = tau - mu*(sig-sp) - cohesion
			if delta:
				X,Y,Z,sxx = self.slice('strs_xx', slice, divisions, time0, method)
				X,Y,Z,syy = self.slice('strs_yy', slice, divisions, time0, method)
				X,Y,Z,szz = self.slice('strs_zz', slice, divisions, time0, method)
				X,Y,Z,sxy = self.slice('strs_xy', slice, divisions, time0, method)
				X,Y,Z,sxz = self.slice('strs_xz', slice, divisions, time0, method)
				X,Y,Z,syz = self.slice('strs_yz', slice, divisions, time0, method)
				X,Y,Z,sp  = self.slice('P',       slice, divisions, time0, method)
				
				px = sxx*nhat[0]+sxy*nhat[1]+sxz*nhat[2]
				py = sxy*nhat[0]+syy*nhat[1]+syz*nhat[2]
				pz = sxz*nhat[0]+syz*nhat[1]+szz*nhat[2]
				
				sig = px*nhat[0]+py*nhat[1]+pz*nhat[2]
				tau = np.sqrt(px**2+py**2+pz**2 - sig**2)
				valsI = valsI - (tau - mu*(sig-sp) - cohesion)
				
		return X, Y, Z, valsI
	def slice_plot_line(self,variable=None,time=None,slice='',divisions=[20,20],labels=False, label_size=10.,levels=10,xlims=[],	
		ylims=[],colors='k',linestyle='-',save='',	xlabel='x / m',ylabel='y / m',title='', font_size='medium', method='nearest',
		equal_axes=True):	
		'''Returns a line plot of contour data. Invokes the ``slice_plot_data`` function to interpolate slice data for plotting.
		
		:param variable: Output data variable, for example 'P' = pressure.
		:type variable: str
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data. If a list of two times is passed, the difference between the first and last is plotted.
		:type time: fl64
		:param slice: List specifying orientation of output slice, e.g., ['x',200.] is a vertical slice at ``x = 200``, ['z',-500.] is a horizontal slice at ``z = -500.``, [point1, point2] is a fixed limit vertical or horizontal domain corresponding to the bounding box defined by point1 and point2.
		:type slice: lst[str,fl64]
		:param divisions: Resolution to supply mesh data.
		:type divisions: [int,int]
		:param method: Method of interpolation, options are 'nearest', 'linear'.
		:type method: str
		:param labels: Specify whether labels should be added to contour plot.
		:type labels: bool
		:param label_size: Specify text size of labels on contour plot, either as an integer or string, e.g., 10, 'small', 'x-large'.
		:type label_size: str, int
		:param levels: Contour levels to plot. Can specify specific levels in list form, or a single integer indicating automatic assignment of levels. 
		:type levels: lst[fl64], int
		:param xlims: Plot limits on x-axis.
		:type xlims: [fl64, fl64]
		:param ylims: Plot limits on y-axis.
		:type ylims: [fl64, fl64]
		:param linestyle: Style of contour lines, e.g., 'k-' = solid black line, 'r:' red dotted line.
		:type linestyle: str
		:param save: Name to save plot. Format specified extension (default .png if none give). Supported extensions: .png, .eps, .pdf.
		:type save: str
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param title: Plot title.
		:type title: str
		:param font_size: Specify text size, either as an integer or string, e.g., 10, 'small', 'x-large'.
		:type font_size: str, int
		:param equal_axes: Specify equal scales on axes.
		:type equal_axes: bool
		
		'''	
		# at this stage, only structured grids are supported
		if time==None: time = self.times[-1]
		delta = False
		if isinstance(time,list) or isinstance(time,np.ndarray):
			if len(time)>1: 
				time0 = np.min(time)
				time = np.max(time)
				delta=True
		return_flag = self._check_inputs(variable,time,slice)
		if return_flag: return
		# gather plot data
		X, Y, Z, valsI = self.slice(variable=variable, time=time, slice=slice, divisions=divisions, method=method)
		if delta:
			X, Y, Z, valsIi = self.slice(variable=variable, time=time0, slice=slice, divisions=divisions, method=method)
			valsI = valsI - valsIi
		plt.clf()
		plt.figure(figsize=[8,8])
		ax = plt.axes([0.15,0.15,0.75,0.75])
		if xlims: ax.set_xlim(xlims)
		if ylims: ax.set_ylim(ylims)
		if equal_axes: ax.set_aspect('equal', 'datalim')
		CS = plt.contour(X,Y,valsI,levels,colors=colors,linestyle=linestyle)
		if labels: plt.clabel(CS,incline=1,fontsize=label_size)
		if xlabel: plt.xlabel(xlabel,size=font_size)		
		if ylabel: plt.ylabel(ylabel,size=font_size)
		if title: plt.title(title,size=font_size)
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
					
		extension, save_fname, pdf = save_name(save,variable=variable,time=time)
		plt.savefig(save_fname, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
		format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: os.system('epstopdf ' + save_fname); os.system(delStr+' ' + save_fname)	
	def slice_plot(self,variable=None,time=None,slice='',divisions=[20,20],levels=10,cbar=False,xlims=[],
		ylims=[],colors='k',linestyle='-',save='',	xlabel='x / m',ylabel='y / m',title='', font_size='medium', method='nearest',
		equal_axes=True,mesh_lines = None,perm_contrasts=None, scale = 1.):		
		'''Returns a filled plot of contour data. Invokes the ``slice_plot_data`` function to interpolate slice data for plotting.
		
		:param variable: Output data variable, for example 'P' = pressure.
		:type variable: str
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data. If a list of two times is passed, the difference between the first and last is plotted.
		:type time: fl64
		:param slice: List specifying orientation of output slice, e.g., ['x',200.] is a vertical slice at ``x = 200``, ['z',-500.] is a horizontal slice at ``z = -500.``, [point1, point2] is a fixed limit vertical or horizontal domain corresponding to the bounding box defined by point1 and point2.
		:type slice: lst[str,fl64]
		:param divisions: Resolution to supply mesh data.
		:type divisions: [int,int]
		:param method: Method of interpolation, options are 'nearest', 'linear'.
		:type method: str
		:param levels: Contour levels to plot. Can specify specific levels in list form, or a single integer indicating automatic assignment of levels. 
		:type levels: lst[fl64], int
		:param cbar: Add colorbar to plot.
		:type cbar: bool
		:param xlims: Plot limits on x-axis.
		:type xlims: [fl64, fl64]
		:param ylims: Plot limits on y-axis.
		:type ylims: [fl64, fl64]
		:param colors: Specify color string for contour levels.
		:type colors: lst[str]
		:param linestyle: Style of contour lines, e.g., 'k-' = solid black line, 'r:' red dotted line.
		:type linestyle: str
		:param save: Name to save plot. Format specified extension (default .png if none give). Supported extensions: .png, .eps, .pdf.
		:type save: str
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param title: Plot title.
		:type title: str
		:param font_size: Specify text size, either as an integer or string, e.g., 10, 'small', 'x-large'.
		:type font_size: str, int
		:param equal_axes: Specify equal scales on axes.
		:type equal_axes: bool
		:param mesh_lines: Superimpose mesh on the plot (line intersections correspond to node positions) according to specified linestyle, e.g., 'k:' is a dotted black line.
		:type mesh_lines: bool
		:param perm_contrasts: Superimpose permeability contours on the plot according to specified linestyle, e.g., 'k:' is a dotted black line. A gradient method is used to pick out sharp changes in permeability.
		:type perm_contrasts: bool
		
		'''	
		# at this stage, only structured grids are supported
		if time==None: 
			if np.min(self.times)<0: time = self.times[0]
			else: time = self.times[-1]
		delta = False
		if isinstance(time,list) or isinstance(time,np.ndarray):
			if len(time)>1: 
				time0 = np.min(time)
				time = np.max(time)
				delta=True
		# if data not available for one coordinate, assume 2-D simulation, adjust slice accordingly
		if 'x' not in self.variables: slice = ['x',None]
		if 'y' not in self.variables: slice = ['y',None]
		if 'z' not in self.variables: slice = ['z',None]
		return_flag = self._check_inputs(variable=variable,time=time,slice=slice)
		if return_flag: return
		# gather plot data
		X, Y, Z, valsI = self.slice(variable=variable, time=time, slice=slice, divisions=divisions, method=method)
		if delta:
			X, Y, Z, valsIi = self.slice(variable=variable, time=time0, slice=slice, divisions=divisions, method=method)
			valsI = valsI - valsIi
		plt.clf()
		plt.figure(figsize=[8,8])
		ax = plt.axes([0.15,0.15,0.75,0.75])
		if xlims: ax.set_xlim(xlims)
		if ylims: ax.set_ylim(ylims)
		if equal_axes: ax.set_aspect('equal', 'datalim')
		if not isinstance(scale,list):
			CS = plt.contourf(X,Y,valsI*scale,levels)
		elif len(scale) == 2:
			CS = plt.contourf(X,Y,valsI*scale[0]+scale[1],levels)
		if xlabel: plt.xlabel(xlabel,size=font_size)		
		if ylabel: plt.ylabel(ylabel,size=font_size)
		if title: plt.title(title,size=font_size)
		if cbar:			
			cbar=plt.colorbar(CS)
			for t in cbar.ax.get_yticklabels():
				t.set_fontsize(font_size)
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
			
		if perm_contrasts:
			if 'perm_x' not in self.variables: 
				print 'WARNING: No permeability data to construct unit boundaries.'
			else:
				X, Y, Z, k = self.slice(variable='perm_z', time=time, slice=slice, divisions=divisions, method=method)
				# calculate derivatives in X and Y directions
				dkdX = np.diff(k,1,0)#/np.diff(Y,1,0)
				dkdX = (dkdX[1:,1:-1]+dkdX[:-1,1:-1])/2
				dkdY = np.diff(k,1,1)#/np.diff(X,1,1)
				dkdY = (dkdY[1:-1,1:]+dkdY[1:-1,:-1])/2
				dk = (abs((dkdX+dkdY)/2)>0.2)*1.
				col = 'k'; ln = '-'
				for let in perm_contrasts:
					if let in ['k','r','g','b','m','c','y','w']: col = let
					if let in ['-','--','-.',':']: ln = let
				CS = plt.contour(X[1:-1,1:-1],Y[1:-1,1:-1],dk,[0.99999999999],colors=col,linestyles=ln)
				
		xlims = ax.get_xlim()
		ylims = ax.get_ylim()
		if mesh_lines:
			# add grid lines
			ax.set_xlim(xlims[0],xlims[1])
			ax.set_ylim(ylims[0],ylims[1])
			if slice[0] == 'z':
				for t in np.unique(self[self.times[0]]['x']):
					ax.plot([t,t],[ylims[0],ylims[1]],mesh_lines,zorder=100)
				for t in np.unique(self[self.times[0]]['y']):
					ax.plot([xlims[0],xlims[1]],[t,t],mesh_lines,zorder=100)
			elif slice[0] == 'x':
				for t in np.unique(self[self.times[0]]['y']):
					ax.plot([t,t],[ylims[0],ylims[1]],mesh_lines,zorder=100)
				for t in np.unique(self[self.times[0]]['z']):
					ax.plot([xlims[0],xlims[1]],[t,t],mesh_lines,zorder=100)
			elif slice[0] == 'y':
				for t in np.unique(self[self.times[0]]['x']):
					ax.plot([t,t],[ylims[0],ylims[1]],mesh_lines,zorder=100)
				for t in np.unique(self[self.times[0]]['z']):
					ax.plot([xlims[0],xlims[1]],[t,t],mesh_lines,zorder=100)
					
		extension, save_fname, pdf = save_name(save=save,variable=variable,time=time)
		plt.savefig(save_fname, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
		format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: os.system('epstopdf ' + save_fname); os.system(delStr+' ' + save_fname)	
	def profile(self, variable, profile, time=None, divisions=30, method='nearest'):
		'''Return variable data along the specified line in 3-D space. If only two points are supplied,
		the profile is assumed to be a straight line between them.
		
		:param variable: Output data variable, for example 'P' = pressure. Can specify multiple variables with a list.
		:type variable: str, lst[str]
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data.
		:type time: fl64
		:param profile: Three column array with each row corresponding to a point in the profile.
		:type profile: ndarray
		:param divisions: Number of points in profile. Only relevant if straight line profile being constructed from two points.
		:type divisions: int
		:param method: Interpolation method, options are 'nearest' (default) and 'linear'.
		:type method: str
		
		:returns: Multi-column array. Columns are in order x, y and z coordinates of profile, followed by requested variables.
		
		'''
		if isinstance(profile,list): profile = np.array(profile)
		if divisions: divisions = int(divisions)
		if time==None: time = self.times[-1]		
		from scipy.interpolate import griddata
		if not isinstance(variable,list): variable = [variable,]
		
		dat = self[time]
		points = np.transpose(np.array([dat['x'],dat['y'],dat['z']]))
		
		if profile.shape[0]==2:
			# construct line profile		
			xrange = np.linspace(profile[0][0],profile[1][0],divisions)
			yrange = np.linspace(profile[0][1],profile[1][1],divisions)
			zrange = np.linspace(profile[0][2],profile[1][2],divisions)
			profile = np.transpose(np.array([xrange,yrange,zrange]))
		
		outpoints = [list(profile[:,0]),list(profile[:,1]),list(profile[:,2])]
		for var in variable:
			vals = np.transpose(np.array(dat[var]))
			valsI = griddata(points,vals,profile,method=method)
			outpoints.append(list(valsI))
		
		return np.array(outpoints).transpose()
	def profile_plot(self,variable=None,time=None, profile=[],divisions = 30,xlims=[],ylims=[],
		color='k',marker='x-',save='',xlabel='distance / m',ylabel='',title='',font_size='medium',method='nearest',
		verticalPlot=False,elevationPlot=False):
		'''Return a plot of the given variable along a specified profile. If the profile comprises two points, 
		these are interpreted as the start and end points of a straight line profile.
		
		:param variable: Output data variable, for example 'P' = pressure. Can specify multiple variables with a list.
		:type variable: str, lst[str]
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data. If a list of two times is passed, the difference between the first and last is plotted.
		:type time: fl64
		:param profile: Three column array with each row corresponding to a point in the profile.
		:type profile: ndarray
		:param divisions: Number of points in profile. Only relevant if straight line profile being constructed from two points.
		:type divisions: int
		:param method: Interpolation method, options are 'nearest' (default) and 'linear'.
		:type method: str
		:param xlims: Plot limits on x-axis.
		:type xlims: [fl64, fl64]
		:param ylims: Plot limits on y-axis.
		:type ylims: [fl64, fl64]
		:param color: Color of profile.
		:type color: str
		:param marker: Style of line, e.g., 'x-' = solid line with crosses, 'o:' dotted line with circles.
		:type marker: str
		:param save: Name to save plot. Format specified extension (default .png if none give). Supported extensions: .png, .eps, .pdf.
		:type save: str
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param title: Plot title.
		:type title: str
		:param font_size: Specify text size, either as an integer or string, e.g., 10, 'small', 'x-large'.
		:type font_size: str, int
		
		:param verticalPlot: Flag to plot variable against profile distance on the y-axis.
		:type verticalPlot: bool
		:param elevationPlot: Flag to plot variable against elevation on the y-axis.
		:type elevationPlot: bool
		
		'''
		if time==None: time = self.times[-1]	
		delta = False
		if isinstance(time,list) or isinstance(time,np.ndarray):
			if len(time)>1: 
				time0 = np.min(time)
				time = np.max(time)
				delta=True
		if not variable: 
			print 'Error: no plot variable specified.'
			print 'Options are'
			for var in self.variables: print var
			return
		if not ylabel: ylabel = variable
		
		plt.clf()
		plt.figure(figsize=[8,8])
		ax = plt.axes([0.15,0.15,0.75,0.75])
		outpts = self.profile(variable=variable,profile=profile,time=time,divisions=divisions,method=method)
		if delta:
			outptsI = self.profile(variable=variable,profile=profile,time=time,divisions=divisions,method=method)
			outpts[:,3] = outpts[:,3] - outptsI[:,3]
		
		x0,y0,z0 = outpts[0,:3]
		x = np.sqrt((outpts[:,0]-x0)**2+(outpts[:,1]-y0)**2+(outpts[:,2]-z0)**2)
		y = outpts[:,3]
		if verticalPlot:
			temp = x; x = y; y = temp
			temp = xlabel; xlabel = ylabel; ylabel = temp
			temp = xlims; xlims = ylims; ylims = temp
		if elevationPlot:
			x = outpts[:,3]
			y = outpts[:,2]
			temp = xlabel; xlabel = ylabel; ylabel = temp
			temp = xlims; xlims = ylims; ylims = temp
		plt.plot(x,y,marker,color=color)
		
		if xlims: ax.set_xlim(xlims)
		if ylims: ax.set_ylim(ylims)
		if xlabel: plt.xlabel(xlabel,size=font_size)		
		if ylabel: plt.ylabel(ylabel,size=font_size)
		if title: plt.title(title,size=font_size)
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
		
		extension, save_fname, pdf = save_name(save,variable=variable,time=time)
		plt.savefig(save_fname, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
		format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: os.system('epstopdf ' + save_fname); os.system(delStr+' ' + save_fname)	
	def cutaway_plot(self,variable=None,time=None,divisions=[20,20,20],levels=10,cbar=False,angle=[45,45],xlims=[],method='nearest',
		ylims=[],zlims=[],colors='k',linestyle='-',save='',	xlabel='x / m', ylabel='y / m', zlabel='z / m', title='', 
		font_size='medium',equal_axes=True,grid_lines=None):
		'''Returns a filled plot of contour data on each of 3 planes in a cutaway plot. Invokes the ``slice_plot_data`` function to interpolate slice data for plotting.
		
		:param variable: Output data variable, for example 'P' = pressure.
		:type variable: str
		:param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data. If a list of two times is passed, the difference between the first and last is plotted.
		:type time: fl64
		:param divisions: Resolution to supply mesh data in [x,y,z] coordinates.
		:type divisions: [int,int,int]
		:param levels: Contour levels to plot. Can specify specific levels in list form, or a single integer indicating automatic assignment of levels. 
		:type levels: lst[fl64], int
		:param cbar: Include colorbar.
		:type cbar: bool
		:param angle: 	View angle of zone. First number is tilt angle in degrees, second number is azimuth. Alternatively, if angle is 'x', 'y', 'z', view is aligned along the corresponding axis.
		:type angle: [fl64,fl64], str
		:param method: Method of interpolation, options are 'nearest', 'linear'.
		:type method: str
		:param xlims: Plot limits on x-axis.
		:type xlims: [fl64, fl64]
		:param ylims: Plot limits on y-axis.
		:type ylims: [fl64, fl64]
		:param zlims: Plot limits on z-axis.
		:type zlims: [fl64, fl64]
		:param colors: Specify color string for contour levels.
		:type colors: lst[str]
		:param linestyle: Style of contour lines, e.g., 'k-' = solid black line, 'r:' red dotted line.
		:type linestyle: str
		:param save: Name to save plot. Format specified extension (default .png if none give). Supported extensions: .png, .eps, .pdf.
		:type save: str
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param zlabel: Label on z-axis.
		:type zlabel: str
		:param title: Plot title.
		:type title: str
		:param font_size: Specify text size, either as an integer or string, e.g., 10, 'small', 'x-large'.
		:type font_size: str, int
		:param equal_axes: Force plotting with equal aspect ratios for all axes.
		:type equal_axes: bool
		:param grid_lines: Extend tick lines across plot according to specified linestyle, e.g., 'k:' is a dotted black line.
		:type grid_lines: bool
		 
		'''	
		# check inputs
		if time==None: time = self.times[-1]
		delta = False
		if isinstance(time,list) or isinstance(time,np.ndarray):
			if len(time)>1: 
				time0 = np.min(time)
				time = np.max(time)
				delta=True
		return_flag = self._check_inputs(variable=variable,time=time,slice=slice)
		if return_flag: return
		# set up axes
		fig = plt.figure(figsize=[11.7,8.275])
		ax = plt.axes(projection='3d')
		ax.set_aspect('equal', 'datalim')
		# make axes equal
		if 'x' not in self.variables or 'y' not in self.variables or 'z' not in self.variables: 
			print 'No xyz data, skipping'; return
		xmin,xmax = np.min(self[time]['x']),np.max(self[time]['x'])
		ymin,ymax = np.min(self[time]['y']),np.max(self[time]['y'])
		zmin,zmax = np.min(self[time]['z']),np.max(self[time]['z'])		
		if equal_axes:
			MAX = np.max([xmax-xmin,ymax-ymin,zmax-zmin])/2
			C = np.array([xmin+xmax,ymin+ymax,zmin+zmax])/2
			for direction in (-1, 1):
				for point in np.diag(direction * MAX * np.array([1,1,1])):
					ax.plot([point[0]+C[0]], [point[1]+C[1]], [point[2]+C[2]], 'w')
		if not xlims: xlims = [xmin,xmax]
		if not ylims: ylims = [ymin,ymax]
		if not zlims: zlims = [zmin,zmax]
		# set view angle
		ax.view_init(angle[0],angle[1])
		
		ax.set_xlabel(xlabel,size=font_size)
		ax.set_ylabel(ylabel,size=font_size)
		ax.set_zlabel(zlabel,size=font_size)

		plt.title(title+'\n\n\n\n',size=font_size)
		
		scale = 1e6
		levels = [l/scale for l in levels]
	
		X, Y, Z, valsI = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[1],ylims[1],zlims[0]]], [divisions[0],divisions[1]], time, method)
		if delta:
			X, Y, Z, valsIi = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[1],ylims[1],zlims[0]]], [divisions[0],divisions[1]], time0, method)
			valsI = valsI - valsIi
		cset = ax.contourf(X, Y, valsI/scale, zdir='z', offset=zlims[0], cmap=cm.coolwarm,levels=levels)
		
		X, Y, Z, valsI = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[0],ylims[1],zlims[1]]], [divisions[1],divisions[2]], time, method)
		if delta:
			X, Y, Z, valsIi = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[0],ylims[1],zlims[1]]], [divisions[1],divisions[2]], time0, method)
			valsI = valsI - valsIi
		cset = ax.contourf(valsI/scale, Y, Z,  zdir='x', offset=xlims[0], cmap=cm.coolwarm,levels=levels)
		
		X, Y, Z, valsI = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[1],ylims[0],zlims[1]]], [divisions[0],divisions[2]], time, method)
		if delta:
			X, Y, Z, valsIi = self.slice(variable, [[xlims[0],ylims[0],zlims[0]],[xlims[1],ylims[0],zlims[1]]], [divisions[0],divisions[2]], time0, method)
			valsI = valsI - valsIi
		cset = ax.contourf(X, valsI/scale, Z,  zdir='y', offset=ylims[0], cmap=cm.coolwarm,levels=levels)

		if cbar:
			cbar=plt.colorbar(cset)
			tick_labels = [str(float(t*scale)) for t in levels]
			cbar.locator = matplotlib.ticker.FixedLocator(levels)
			cbar.formatter = matplotlib.ticker.FixedFormatter(tick_labels)
			cbar.update_ticks()

		if grid_lines:
			# add grid lines
			ax.set_xlim(xlims[0],xlims[1])
			ax.set_ylim(ylims[0],ylims[1])
			ax.set_zlim(zlims[0],zlims[1])
			xticks = ax.get_xticks()
			yticks = ax.get_yticks()
			zticks = ax.get_zticks()

			off = 0.
			for t in xticks:
				ax.plot([t,t],[ylims[0],ylims[1]],[zlims[0]+off,zlims[0]+off],grid_lines,zorder=100)
				ax.plot([t,t],[ylims[0]+off,ylims[0]+off],[zlims[0],zlims[1]],grid_lines,zorder=100)
			for t in yticks:
				ax.plot([xlims[0],xlims[1]],[t,t],[zlims[0]+off,zlims[0]+off],grid_lines,zorder=100)
				ax.plot([xlims[0]+off,xlims[0]+off],[t,t],[zlims[0],zlims[1]],grid_lines,zorder=100)
			for t in zticks:
				ax.plot([xlims[0],xlims[1]],[ylims[0]+off,ylims[0]+off],[t,t],grid_lines,zorder=100)
				ax.plot([xlims[0]+off,xlims[0]+off],[ylims[0],ylims[1]],[t,t],grid_lines,zorder=100)
		
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)	
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)	
		for t in ax.get_zticklabels():
			t.set_fontsize(font_size)	
		
		ax.set_xlim(xlims)
		ax.set_ylim(ylims)
		ax.set_zlim(zlims)
		
		extension, save_fname, pdf = save_name(save=save,variable=variable,time=time)
		plt.savefig(save_fname, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
		format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
	def node(self,node,time=None,variable=None):
		'''Returns all information for a specific node.
		
		If time and variable not specified, a dictionary of time series is returned with variables as the dictionary keys.
		
		If only time is specified, a dictionary of variable values at that time is returned, with variables as dictionary keys.
		
		If only variable is specified, a time series vector is returned for that variable.
		
		If both time and variable are specified, a single value is returned, corresponding to the variable value at that time, at that node.
		
		:param node: Node index for which variable information required.
		:type node: int 
		:param time: Time at which variable information required. If not specified, all output.
		:type time: fl64
		:param variable: Variable for which information requested. If not specified, all output.
		:type variable: str
		
		'''
		if 'n' not in self.variables: print 'Node information not available'; return
		nd = np.where(self[self.times[0]]['n']==node)[0][0]
		if time == None and variable == None:
			ks = copy(self.variables); ks.remove('n')
			outdat = dict([(k,[]) for k in ks])
			for t in self.times:
				dat = self[t]
				for k in outdat.keys():
					outdat[k].append(dat[k][nd])
		elif time == None:
			if variable not in self.variables: print 'ERROR: no variable by that name'; return
			outdat = []
			for t in self.times:
				dat = self[t]
				outdat.append(dat[variable][nd])
		elif variable == None:
			ks = copy(self.variables); ks.remove('n')
			outdat = dict([(k,self[time][k][nd]) for k in ks])			
		else:
			outdat = self[time][variable][nd]
		return outdat
	def _get_variables(self): return self._variables
	variables = property(_get_variables)#: (*lst[str]*) List of variables for which output data are available.
	def _get_format(self): return self._format
	format = property(_get_format) #: (*str*) Format of output file, options are 'tec', 'surf', 'avs' and 'avsx'.
	def _get_filename(self): return self._filename
	filename = property(_get_filename)  #: (*str*) Name of FEHM contour output file. Wildcards can be used to define multiple input files.
	def _get_times(self): return np.sort(self._times)
	times = property(_get_times)	#: (*lst[fl64]*) List of times (in seconds) for which output data are available.
	def _get_information(self):
		print 'FEHM contour output - format '+self._format
		print '    call format: fcontour[time][variable][node_index-1]'
		prntStr =  '    times ('+str(len(self.times))+'): '
		for time in self.times: prntStr += str(time)+', '
		print prntStr[:-2]+' days'
		prntStr = '    variables: '
		for var in self.variables: prntStr += str(var)+', '
		print prntStr
	what = property(_get_information) #:(*str*) Print out information about the fcontour object.
class fhistory(object):						# Reading and plotting methods associated with history output data.
	'''History output information object.
	
	'''
	def __init__(self,filename=None,verbose=True):
		self._filename=None	
		self._times=[]	
		self._verbose = verbose
		self._data={}
		self._row=None
		self._nodes=[]	
		self._variables=[] 
		self._keyrows={}
		self.column_name=[]
		self.num_columns=0
		self._nkeys=1
		if filename: self._filename=filename; self.read(filename)
	def __getitem__(self,key):
		if key in self.variables:
			return self._data[key]
		else: return None
	def __repr__(self): 
		retStr =  'History output for variables '
		for var in self.variables:
			retStr += var+', '
		retStr = retStr[:-2] + ' at '
		if len(self.nodes)>10:
			retStr += str(len(self.nodes)) + ' nodes.'
		else:
			if len(self.nodes)==1:
				retStr += 'node '
			else:
				retStr += 'nodes '
			for nd in self.nodes:
				retStr += str(nd) + ', '
			retStr = retStr[:-2] + '.'
		return retStr
	def read(self,filename): 						# read contents of file
		'''Read in FEHM history output information.
		
		:param filename: File name for output data, can include wildcards to define multiple output files.
		:type filename: str
		'''
		from glob import glob
		files=glob(filename)
		configured=False
		for i,fname in enumerate(files):
			if self._verbose: print fname
			self._file=open(fname,'rU')
			header=self._file.readline()
			if header.strip()=='': continue				# empty file
			self._detect_format(header)
			if self.format=='tec': 
				header=self._file.readline()
				if header.strip()=='': continue 		# empty file
				i = 0; sum_file = False
				while not header.startswith('variables'): 
					header=self._file.readline()
					i = i+1
					if i==10: sum_file=True; break
				if sum_file: continue
				self._setup_headers_tec(header)
			elif self.format=='surf': 
				self._setup_headers_surf(header)
			elif self.format=='default': 
				header=self._file.readline()
				header=self._file.readline()
				if header.strip()=='': continue 		# empty file
				i = 0; sum_file = False
				while not header.startswith('Time '): 
					header=self._file.readline()
					i = i+1
					if i==10: sum_file=True; break
				if sum_file: continue
				self._setup_headers_default(header)
			else: print 'Unrecognised format'; return
			if not configured:
				self.num_columns = len(self.nodes)+1
			if self.num_columns>0: configured=True
			if self.format=='tec':
				self._read_data_tec(fname.split('_')[-2])
			elif self.format=='surf':
				self._read_data_surf(fname.split('_')[-2])
			elif self.format=='default':
				self._read_data_default(fname.split('_')[-1].split('.')[0])
			self._file.close()
	def _detect_format(self,header):
		if header.startswith('TITLE'):
			self._format = 'tec'
		elif header.startswith('Time '):
			self._format = 'surf'
		else:
			self._format = 'default'
	def _setup_headers_tec(self,header):
		header=header.split('" "Node')
		if self.nodes: return
		for key in header[1:-1]: self._nodes.append(int(key))
		self._nodes.append(int(header[-1].split('"')[0]))
	def _setup_headers_surf(self,header):
		header=header.split(', Node')
		if self.nodes: return
		for key in header[1:]: self._nodes.append(int(key))
	def _setup_headers_default(self,header):
		header=header.split(' Node')
		if self.nodes: return
		for key in header[1:]: self._nodes.append(int(key))
	def _read_data_tec(self,var_key):
		self._variables.append(hist_var_names[var_key])
		lns = self._file.readlines()
		i = 0
		while lns[i].startswith('text'): i+=1
		data = []
		for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
		data = np.array(data)
		if data[-1,0]<data[-2,0]: data = data[:-1,:]
		self._times = np.array(data[:,0])
		self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
	def _read_data_surf(self,var_key):
		self._variables.append(hist_var_names[var_key])
		lns = self._file.readlines()
		data = []
		for ln in lns: data.append([float(d) for d in ln.strip().split(',')])
		data = np.array(data)
		if data[-1,0]<data[-2,0]: data = data[:-1,:]
		self._times = np.array(data[:,0])
		self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
	def _read_data_default(self,var_key):
		self._variables.append(hist_var_names[var_key])
		lns = self._file.readlines()
		data = []
		for ln in lns: data.append([float(d) for d in ln.strip().split()])
		data = np.array(data)
		if data[-1,0]<data[-2,0]: data = data[:-1,:]
		self._times = np.array(data[:,0])
		self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
	def time_plot(self, variable=None, node=0, t_lim=[],var_lim=[],marker='x-',color='k',save='',xlabel='',ylabel='',
		title='',font_size='medium',scale=1.,scale_t=1.): 		# produce a time plot
		'''Generate and save a time series plot of the history data.
		
		:param variable: Variable to plot.
		:type variable: str
		:param node: Node number to plot.
		:type node: int
		:param t_lim: Time limits on x axis.
		:type t_lim: lst[fl64,fl64]
		:param var_lim: Variable limits on y axis.
		:type var_lim: lst[fl64,fl64]
		:param marker: String denoting marker and linetype, e.g., ':s', 'o--'. Default is 'x-' (solid line with crosses).
		:type marker: str
		:param color: String denoting color. Default is 'k' (black).
		:type color: str
		:param save: Name to save plot.
		:type save: str
		:param xlabel: Label on x axis.
		:type xlabel: str
		:param ylabel: Label on y axis.
		:type ylabel: str
		:param title: Title of plot.
		:type title: str
		:param font_size: Font size for axis labels.
		:type font_size: str
		:param scale: If a single number is given, then the output variable will be multiplied by this number. If a two element list is supplied then the output variable will be transformed according to y = scale[0]*x+scale[1]. Useful for transforming between coordinate systems.
		:type scale: fl64
		:param scale_t: As for scale but applied to the time axis.
		:type scale_t: fl64
		'''
		if not node: print 'Error: no plot node specified.'; return
		if not variable: 
			print 'Error: no plot variable specified.'
			print 'Options are'
			for var in self.variables: print var
			return
		if not node: 
			print 'Error: no plot time specified.'
			print 'Options are'
			for node in self.nodes: print node
			return
		
		plt.clf()
		plt.figure(figsize=[8,8])
		ax = plt.axes([0.15,0.15,0.75,0.75])
		if not isinstance(scale,list):
			if not isinstance(scale_t,list):
				plt.plot(self.times*scale_t,self[variable][node]*scale,marker)
			elif len(scale_t) == 2:
				plt.plot(self.times*scale_t[0]+scale_t[1],self[variable][node]*scale,marker)
		elif len(scale) == 2:
			if not isinstance(scale_t,list):
				plt.plot(self.times*scale_t,self[variable][node]*scale[0]+scale[1],marker)
			elif len(scale_t) == 2:
				plt.plot(self.times*scale_t[0]+scale_t[1],self[variable][node]*scale[0]+scale[1],marker)
		if t_lim: ax.set_xlim(t_lim)
		if var_lim: ax.set_ylim(var_lim)
		if xlabel: plt.xlabel(xlabel,size=font_size)		
		if ylabel: plt.ylabel(ylabel,size=font_size)
		if title: plt.title(title,size=font_size)
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
		
		extension, save_fname, pdf = save_name(save,variable=variable,node=node)
		plt.savefig(save_fname, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
		format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: os.system('epstopdf ' + save_fname); os.system(delStr+' ' + save_fname)	
	def _get_variables(self): return self._variables
	variables = property(_get_variables)#: (*lst[str]*) List of variables for which output data are available.
	def _get_format(self): return self._format
	format = property(_get_format) #: (*str*) Format of output file, options are 'tec', 'surf', 'avs' and 'avsx'.
	def _get_filename(self): return self._filename
	filename = property(_get_filename)  #: (*str*) Name of FEHM contour output file. Wildcards can be used to define multiple input files.
	def _get_times(self): return np.sort(self._times)
	times = property(_get_times)	#: (*lst[fl64]*) List of times (in seconds) for which output data are available.
	def _get_nodes(self): return self._nodes
	nodes = property(_get_nodes)	#: (*lst[fl64]*) List of node indices for which output data are available.
	def _get_information(self):
		print 'FEHM history output - format '+self._format
		print '    call format: fhistory[variable][node][time_index]'
		prntStr = '    nodes: '
		for nd in self.nodes: prntStr += str(nd)+', '
		print prntStr
		prntStr =  '    times ('+str(len(self.times))+'): '
		for time in self.times: prntStr += str(time)+', '
		print prntStr[:-2]+' days'
		prntStr = '    variables: '
		for var in self.variables: prntStr += str(var)+', '
		print prntStr
	what = property(_get_information) #:(*str*) Print out information about the fhistory object.
class fzoneflux(fhistory): 					# Derived class of fhistory, for zoneflux output
	'''Zone flux history output information object.
	'''
	def __init__(self,filename=None,verbose=True):
		super(fzoneflux,self).__init__(filename=None, verbose = True)
		self._filename=None	
		self._times=[]	
		self._verbose = verbose
		self._data={}
		self._row=None
		self._zones=[]	
		self._variables=[] 
		self._keyrows={}
		self.column_name=[]
		self.num_columns=0
		self._nkeys=1
		if filename: self._filename=filename; self.read(filename)
	def _setup_headers_tec(self,header):
		'placeholder'
	def _read_data_tec(self,var_key):
		zn = int(var_key[-5:])
		if var_key.startswith('c'):
			if zn not in self._zones: self._zones.append(zn)
			if 'co2_source' not in self._variables:
				self._variables += flxz_co2_names
				for var in flxz_co2_names: self._data[var] = {}
				
			lns = self._file.readlines()
			i = 0
			while lns[i].startswith('text'): i+=1
			data = []
			for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
			data = np.array(data)
			if data[-1,0]<data[-2,0]: data = data[:-1,:]
			self._times = np.array(data[:,0])
			for j,var_key in enumerate(flxz_co2_names):
				self._data[var_key].update(dict([(zn,data[:,j+1])]))
		elif var_key.startswith('w'):
			if zn not in self._zones: self._zones.append(zn)
			if 'water_source' not in self._variables:
				self._variables += flxz_water_names
				for var in flxz_water_names: self._data[var] = {}
				
			lns = self._file.readlines()
			i = 0
			while lns[i].startswith('text'): i+=1
			data = []
			for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
			data = np.array(data)			
			if data[-1,0]<data[-2,0]: data = data[:-1,:]
			self._times = np.array(data[:,0])
			for j,var_key in enumerate(flxz_water_names):
				self._data[var_key].update(dict([(zn,data[:,j+1])]))
		elif var_key.startswith('v'):
			if zn not in self._zones: self._zones.append(zn)
			if 'vapor_source' not in self._variables:
				self._variables += flxz_vapor_names
				for var in flxz_vapor_names: self._data[var] = {}
				
			lns = self._file.readlines()
			i = 0
			while lns[i].startswith('text'): i+=1
			data = []
			for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
			data = np.array(data)			
			if data[-1,0]<data[-2,0]: data = data[:-1,:]
			self._times = np.array(data[:,0])
			for j,var_key in enumerate(flxz_vapor_names):
				self._data[var_key].update(dict([(zn,data[:,j+1])]))
	def _read_data_surf(self,var_key):
		self._variables.append(hist_var_names[var_key])
		lns = self._file.readlines()
		data = []
		for ln in lns[i:]: data.append([float(d) for d in ln.strip().split(',')])
		data = np.array(data)			
		if data[-1,0]<data[-2,0]: data = data[:-1,:]
		self._times = np.array(data[:,0])
		self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
	def _read_data_default(self,var_key):
		self._variables.append(hist_var_names[var_key])
		lns = self._file.readlines()
		data = []
		for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
		data = np.array(data)
		if data[-1,0]<data[-2,0]: data = data[:-1,:]
		self._times = np.array(data[:,0])
		self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
	def _get_zones(self): return self._zones
	def _set_zones(self,value): self._zones = value
	zones = property(_get_zones, _set_zones) #: (*lst[int]*) List of zone indices for which output data are available.
class multi_pdf(object):
	'''Tool for making a single pdf document from multiple eps files.'''
	def __init__(self,combineString = 'gswin64',
		save='multi_plot.pdf',files = [],delete_files = True):
		self.combineString = combineString
		self._save = save
		self._delete_files = delete_files
		self._assign_files(files)
	def _assign_files(self,files):
		if files == []: self._files = ImmutableDict({})
		if isinstance(files,list):
			self._files = ImmutableDict(dict([(i+1,file) for i,file in enumerate(files)]))
		elif isinstance(files,dict):
			ks = files.keys()
			for k in ks:
				if not isinstance(k,int):print 'ERROR: Dictionary keys must be integers.'; return
			self._files = ImmutableDict(files)
		elif isinstance(files,str):
			self._files = ImmutableDict(dict(((1,files),)))
	def add(self,filename,pagenum=None):
		'''Add a new page. If a page number is specified, the page will replace the current. 
		Otherwise it will be appended to the end of the document.
		
		:param filename: Name of .eps file to be added.
		:type filename: str
		:param pagenum: Page number of file to be added.
		:type pagenum: int
		
		'''
		if len(filename.split('.'))==1: filename += '.eps'
		if not os.path.isfile(filename): print 'WARNING: '+filename+' not found.'; return		
		if not filename.endswith('.eps'): print 'WARNING: Non EPS format not supported.'
		
		if pagenum and pagenum in self.files.keys():
			print 'WARNING: Replacing '+self.files[pagenum]
			self.files[pagenum] = filename
		else: 
			if not pagenum: pagenum = self._pagemax+1
			self._files.update(dict(((pagenum,filename),)))			
	def insert(self,filename,pagenum):
		'''Insert a new page at the given page number.
		
		:param filename: Name of .eps file to be inserted.
		:type filename: str
		:param pagenum: Page number of file to be inserted.
		:type pagenum: int
		'''
		if len(filename.split('.'))==1: filename += '.eps'
		if not os.path.isfile(filename): print 'WARNING: '+filename+' found.'; return
		if not filename.endswith('.eps'): print 'WARNING: Non EPS format not supported.'
			
		if pagenum > self._pagemax: self.add(filename); return
		ks = self._files.keys()
		self._files = ImmutableDict(dict([(k,self._files[k]) for k in ks if k < pagenum]+
		[(pagenum,filename)]+[(k+1,self._files[k]) for k in ks if k >= pagenum]))			
	def make(self):
		'''Construct the pdf.'''
		cs = self.combineString + ' -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile='+self.save
		for i in np.sort(self.files.keys()):
			if not self.files[i].endswith('.eps'): print 'WARNING: Cannot combine '+self.files[i]+'. EPS format required. Skipping...'; continue
			if len(self.files[i].split()) != 1:
				cs += ' "'+self.files[i]+'"'
			else:
				cs += ' '+self.files[i]
		os.system(cs)
		for i in np.sort(self.files.keys()): 
			if len(self.files[i].split()) != 1:
				os.system(delStr+' "'+self.files[i]+'"')
			else:
				os.system(delStr+' '+self.files[i])
	def _get_combineString(self): return self._combineString
	def _set_combineString(self,value): self._combineString = value
	combineString = property(_get_combineString, _set_combineString) #: (*str*)	Command line command, with options, generate pdf from multiple eps files. See manual for further instructions.
	def _get_files(self): return self._files
	files = property(_get_files) #: (*lst[str]*) List of eps files to be assembled into pdf.
	def _get_pagemax(self): 
		ks = self._files.keys()
		for k in ks:
			if not isinstance(k,int): print 'ERROR: Non integer dictionary key'; return
		if len(ks) == 0: return 0
		return np.max(ks)
	_pagemax = property(_get_pagemax)
	def _get_save(self): return self._save
	def _set_save(self,value): self._save = value
	save = property(_get_save, _set_save) #: (*str*) Name of the final pdf to output.