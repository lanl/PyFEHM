"""For manipulating FEHM grids."""

import numpy as np
import os,math
from ftool import*
from time import time
try:
	from matplotlib import pyplot as plt
	from mpl_toolkits.mplot3d import axes3d
except ImportError:
	'placeholder'
from fpost import *
from copy import deepcopy
import platform

WINDOWS = platform.system()=='Windows'
if WINDOWS: copyStr = 'copy'; delStr = 'del'; slash = '\\'
else: copyStr = 'cp'; delStr = 'rm'; slash = '/'

node_props = ('kx','ky','kz','cond_x','cond_y','cond_z','density','specific_heat','porosity','thermal_expansion','pressure_coupling',
'youngs_modulus','poissons_ratio')

nbr_update = dict((
	(0,[[0,2,4],[1,2,4]]),
	(1,[[1,2,4],[0,3,5]]),
	(2,[[0,3,4],[3,0,6]]),
	(3,[[1,3,4],[1,2,7]]),
	(4,[[0,2,5],[5,6,0]]),
	(5,[[1,2,5],[4,7,1]]),
	(6,[[0,3,5],[7,4,2]]),
	(7,[[1,3,5],[6,5,3]])
	))

class fnode(object):				#Node object.
	""" FEHM grid node object.
	"""
	__slots__ = ['_index','_position','_connections','_elements','_variable','_material','_generator','_zone',
		'_permeability','_conductivity','_density','_specific_heat','_porosity','_youngs_modulus','_poissons_ratio','_thermal_expansion',
		'_pressure_coupling','_Pi','_Ti','_Si','_S_co2gi','_S_co2li','_co2_aqi','_strsi','_dispi','_P','_T','_S',
		'_S_co2g','_S_co2l','_co2_aq','_strs','_disp']
	def __init__(self,index=None,position=None):		
		self._index = index			
		self._position=position	
		self._connections = []	
		self._elements = []		
		self._variable = ImmutableDict({}) 		
		self._material = ImmutableDict({}) 		
		self._generator = ImmutableDict({}) 		
		self._zone = ImmutableDict({}) 			
		self._permeability = None
		self._conductivity = None
		self._density = None
		self._specific_heat = None
		self._porosity = None
		self._youngs_modulus = None
		self._poissons_ratio = None
		self._thermal_expansion = None
		self._pressure_coupling = None
		self._Pi = None
		self._Ti = None
		self._Si = None
		self._S_co2gi = None
		self._S_co2li = None
		self._co2_aqi = None
		self._strsi = None
		self._dispi = None	
		self._P = None
		self._T = None
		self._S = None
		self._S_co2g = None
		self._S_co2l = None
		self._co2_aq = None
		self._strs = None
		self._disp = None		
	def __repr__(self): return 'n'+str(self.index)
	def _get_index(self): return self._index
	index = property(_get_index) #: (*int*) Integer number denoting the node.	
	def _get_position(self): return self._position
	position = property(_get_position) #: (*lst[fl64]*) List of the node's coordinates in 2- or 3-D space.
	def _get_zone(self): return self._zone
	def _set_zone(self,value): self._zone = value
	zone = property(_get_zone, _set_zone) #: (*dict*) Dictionary of zones to which the node belongs.
	def _get_generator(self): return self._generator
	def _set_generator(self,value): self._generator = value
	generator = property(_get_generator, _set_generator) #: (*dict*) Dictionary of generator properties associated with node.
	def _get_variable(self): return self._variable
	def _set_variable(self,value): self._variable = value
	variable = property(_get_variable, _set_variable) #: (*dict*) Dictionary of thermodynamic properties associated with node. Only accessible when initial conditions files are loaded in conjunction with the input file and grid.
	def _get_material(self): return self._material
	def _set_material(self,value): self._material = value
	material = property(_get_material, _set_material) #: (*dict*) Dictionary of material properties associated with node.
	def _get_info(self):
		prntStr='\n Node number: '+str(self.index)+'\n'
		prntStr+='\nGeometric properties.......\n'
		prntStr+='          x = '+str(self.position[0])+'\n'
		prntStr+='          y = '+str(self.position[1])+'\n'
		prntStr+='          z = '+str(self.position[2])+'\n'
		prntStr+=' neighbours = ['
		for nd in self.connected_nodes:	prntStr+=str(nd.index)+','
		prntStr = prntStr[:-1]+']\n'
		if self.zone:
			prntStr+='\nZone membership............\n'
			for zn in self.zonelist:
				prntStr+='    '+str(zn.index)
				if zn.name:
					prntStr+=': '+zn.name
				prntStr+='\n'
		if self.variable:
			prntStr+='\nState......................\n'
			keys = self.variable.keys()
			if 'T' in keys: prntStr+='temperature = '+str(self.variable['T'])+'\n'
			if 'P' in keys: prntStr+='   pressure = '+str(self.variable['P'])+'\n'
			if 'S' in keys: prntStr+='  sat water = '+str(self.variable['S'])+'\n'
			if 'S_co2g' in keys: prntStr+='sat co2 gas = '+str(self.variable['S_co2g'])+'\n' 
			if 'S_co2l' in keys: prntStr+='sat co2 liq = '+str(self.variable['S_co2l'])+'\n' 
			if 'co2aq' in keys: prntStr+='dissolved co2 = '+str(self.variable['co2aq'])+'\n' 
			if 'eos' in keys: prntStr+='water phase = '+str(self.variable['eos'])+'\n' 
			if 'co2_eos' in keys: prntStr+='  co2 phase = '+str(self.variable['co2_eos'])+'\n' 
		if self.material:
			prntStr+='\nMaterial properties........\n'
			ks = self.material.keys()
			for k in node_props:
				if k in ks: prntStr += '    '+k + ' = ' +str(self.material[k])+'\n'
		if self.generator:
			prntStr+='\nGenerator properties.......\n'
			ks = self.generator.keys()
			for k in ks:
				prntStr += '    '+k + ' = ' +str(self.generator[k])+'\n'
		print prntStr
	what = property(_get_info)							#: Print to screen information about the node.
	def _get_connected_nodes(self):
		ndlist = []
		for con in self.connections:
			for nd in con.nodes:
				if nd != self: ndlist.append(nd)
		return ndlist
	connected_nodes = property(_get_connected_nodes)		#: (*lst[fnode]*) List of node objects connected to this node.
	def _get_connections(self): return self._connections
	connections = property(_get_connections)#: (*lst[fconn]*) List of connection objects of which the node is a member.
	def _get_elements(self): return self._elements
	elements = property(_get_elements)#: (*lst[felem]*) List of element objects of which the node is a member.
	def _get_zonelist(self): return [self.zone[k] for k in self.zone.keys()]
	zonelist = property(_get_zonelist)	#: (*lst[fzone]*) List of zones of which the node is a member
	def _get_control_volume(self):
		return None
	ctrl_vol = property(_get_control_volume)				#: (*fl64*) Control volume associated with the node *** NOT DONE ***.
	def _get_permeability(self): return self._permeability
	permeability = property(_get_permeability) #: (**)
	def _get_conductivity(self): return self._conductivity
	conductivity = property(_get_conductivity) #: (**)
	def _get_density(self): return self._density
	density = property(_get_density) #: (**)
	def _get_specific_heat(self): return self._specific_heat
	specific_heat = property(_get_specific_heat) #: (**)
	def _get_porosity(self): return self._porosity
	porosity = property(_get_porosity) #: (**)
	def _get_youngs_modulus(self): return self._youngs_modulus
	youngs_modulus = property(_get_youngs_modulus) #: (**)
	def _get_poissons_ratio(self): return self._poissons_ratio
	poissons_ratio = property(_get_poissons_ratio) #: (**)
	def _get_thermal_expansion(self): return self._thermal_expansion
	thermal_expansion = property(_get_thermal_expansion) #: (**)
	def _get_pressure_coupling(self): return self._pressure_coupling
	pressure_coupling = property(_get_pressure_coupling) #: (**)
	def _get_Pi(self): return self._Pi
	Pi = property(_get_Pi) #: (**)
	def _get_Ti(self): return self._Ti
	Ti = property(_get_Ti) #: (**)
	def _get_Si(self): return self._Si
	Si = property(_get_Si) #: (**)
	def _get_S_co2gi(self): return self._S_co2gi
	S_co2gi = property(_get_S_co2gi) #: (**)
	def _get_S_co2li(self): return self._S_co2li
	S_co2li = property(_get_S_co2li) #: (**)
	def _get_co2_aqi(self): return self._co2_aqi
	co2_aqi = property(_get_co2_aqi) #: (**)
	def _get_strsi(self): return self._strsi
	strsi = property(_get_strsi) #: (**)
	def _get_dispi(self): return self._dispi
	dispi = property(_get_dispi) #: (**)
	def _get_P(self): return self._P
	P = property(_get_P) #: (**)
	def _get_T(self): return self._T
	T = property(_get_T) #: (**)
	def _get_S(self): return self._S
	S = property(_get_S) #: (**)
	def _get_S_co2g(self): return self._S_co2g
	S_co2g = property(_get_S_co2g) #: (**)
	def _get_S_co2l(self): return self._S_co2l
	S_co2l = property(_get_S_co2l) #: (**)
	def _get_co2_aq(self): return self._co2_aq
	co2_aq = property(_get_co2_aq) #: (**)
	def _get_strs(self): return self._strs
	strs = property(_get_strs) #: (**)
	def _get_disp(self): return self._disp
	disp = property(_get_disp) #: (**)
class fconn(object):				#Connection object.
	"""Connection object, comprising two connected nodes, separated by some distance.

	A connection is associated with a distance between the two nodes.
	"""
	def __init__(self,nodes=[fnode(),fnode()]):
		self._nodes = nodes		
	def __repr__(self):	return 'n'+str(self.nodes[0].index)+':n'+str(self.nodes[1].index)	
	def _get_distance(self): 
		pos1 = self.nodes[0].position; pos2 = self.nodes[1].position
		return np.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)
	distance = property(_get_distance)	#: (*fl64*) Distance between the two connected nodes.
	def _get_nodes(self): return self._nodes
	nodes = property(_get_nodes)#: (*lst[fnode]*) List of node objects (fnode()) that define the connection.
class felem(object):				#Element object.
	"""Finite element object, comprising a set of connected nodes.
	
	A finite element is associated with an element centre and an element volume.
	"""
	def __init__(self,index=None,nodes = []):
		self._index = index			
		self._nodes = nodes			
	def __repr__(self): 
		retStr = 'e'+str(self.index)+': '
		for nd in self.nodes: retStr+='n'+str(nd.index)+', '
		return retStr
	def _get_index(self): return self._index
	index = property(_get_index)#: (*int*) Integer number denoting the element.		
	def _get_centre(self):
		c = np.array([0,0,0])
		for nd in self.nodes: c += np.array(nd.position)
		c = c/len(self.nodes)
		return c
	centre = property(_get_centre)	#: (*ndarray*) Coordinates of the element centroid.
	def get_info(self):
		print 'Element number: '+str(self.index)
		print '      Centroid: '
		print '            x = '+str(self.centre[0])
		print '            y = '+str(self.centre[1])
		print '            z = '+str(self.centre[2])
		prntStr =  'Contains nodes: ['
		for nd in self.nodes: prntStr += str(nd.index) +','
		prntStr = prntStr[:-1]+']'
		print prntStr
	what = property(get_info)				#: Print to screen information about the element.
	def _get_nodes(self): return self._nodes
	nodes = property(_get_nodes)#: (*lst[fnode]*) List of node objects that define the element.		
	def _get_volume(self):
		return None
	vol = property(_get_volume)		#: (*fl64*) Volume of the element *** NOT DONE ***.		
class octree(object):				#Octree object.
	"""Octree for spatial searching in 3D grids.
	
	The octree is automatically constructed when a new grid file is parsed.
	
	Note: direct interaction with this method is not generally required (or advised!).
	
	The octree is generated recursively, by sub-dividing the model domain into smaller and smaller compartments.
	If a compartment contains more than one point then it continues to sub-divide until there is no more than one point
	per compartment.
	"""
	def __init__(self,bounds,elements,parent=None):
		self.parent=parent		#: Parent compartment to which the compartment belongs.
		self.bounds=bounds		#: Bounding box defining the compartment.
		self.elements=elements	#: (*fnode*) Objects contained in the compartment.
		self.child=[] 			#: Children compartments contained within the compartment.
		self.neighbour = [None,None,None,None,None,None]	#: Six element list containing, in order, compartments +x, -x, +y, -y, +z, -z
		if self.parent: 		# if has a parent, is a child, iterate generation, inherit all elements
			self.generation=self.parent.generation+1	#: Generation of the compartment (a measure of how many sub-divisions have occurred.)
			self.all_elements=self.parent.all_elements
		else: 					# if not a parent, generation 0, elements only those existing
			self.generation=0
			self.all_elements=set(elements)
		if self.num_elements>1: 	# if more than one element in a cube, sub-divide
			cubes=sub_cubes(self.bounds)
			# update neighbour locations
			cube_elements=[[],[],[],[],[],[],[],[]]
			for elt in self.elements:
				for icube,cube in enumerate(cubes):
					if in_cube(elt.position,cube):
						cube_elements[icube].append(elt)
						break
			nbrs=[]
			for cube,elts in zip(cubes,cube_elements):
				if len(elts)>0: 
					self.child.append(octree(cube,elts,self))
					nbrs.append(1)
				else: nbrs.append(0)
			self.assign_neighbours(nbrs)
	def __repr__(self): return self.bounds.__repr__()
	def assign_neighbours(self,neighbours):
		cube_inds = []
		for cube_ind,exist in enumerate(neighbours):
			if exist: cube_inds.append(cube_ind)
		for child,update in zip(self.child,cube_inds):
			child.neighbour = deepcopy(self.neighbour)
			posI,nbrs = nbr_update[update]
			for pos,nbr in zip(posI,nbrs):
				if nbr in cube_inds:
					child.neighbour[pos] = self.child[cube_inds.index(nbr)]
			
	def get_num_elements(self): return len(self.elements)
	num_elements=property(get_num_elements)	#: (*int*) Number of objects in the compartment.
	def get_num_children(self): return len(self.child)
	num_children=property(get_num_children) #: (*int*) Number of sub-compartments in the compartment.
	def leaf(self,pos): 							# find leaf corresponding to position
		if in_cube(pos,self.bounds):
			for child in self.child:
				childleaf=child.leaf(pos)
				if childleaf: return childleaf
			return self
		else: return None
	def min_dist(self,pos):
		min_dist = 1.e10
		if self.child:
			for child in self.child:
				dist,cube = child.min_dist(pos)
				if dist<min_dist: min_dist = dist; min_cube = cube
		else:
			min_dist = (np.sqrt((self.elements[0].position[0]-pos[0])**2+
				(self.elements[0].position[1]-pos[1])**2+
				(self.elements[0].position[2]-pos[2])**2))
			min_cube = self
		#print min_cube.elements[0].index, ' ', min_dist
		return min_dist,min_cube
class fgrid(object):				#Grid object.
	""" FEHM grid object.
	
	"""
	def __init__(self,full_connectivity=False):
		self._nodelist=[]			
		self._node={}				
		self._connlist=[]			
		self._conn={}				
		self._elemlist=[]			
		self._elem={}				
		self._node_tree=None		
		self._filename=''			
		self._dimensions=3
		self._parent = None
		self._full_connectivity = full_connectivity
	def __repr__(self): return self.filename
	def read(self,meshfilename,full_connectivity=False): 
		"""Read data from an FEHM grid file.

		:param meshfilename: name of FEHM grid file, including path specification.
		:type meshfilename: str
		:param full_connectivity: read element and conection data and construct corresponding objects. Consumes siginifcant time and memory and is rarely used.
		:type full_connectivity: bool
		"""
		self._full_connectivity = full_connectivity
		self._filename = meshfilename 
		if meshfilename.endswith('.stor'):
			self.read_stor(meshfilename)
		else:
			self._read_inp(meshfilename)
			self.add_nodetree()
		if self._parent: self._parent._add_boundary_zones()
	def _read_inp(self,meshfilename): 		#Read in fehm meshfile for node,element data .
		self._nodelist = []
		self._connlist = []
		self._elemlist = []
		infile = open(meshfilename)
		if self._parent:
			if slash in meshfilename: meshfilename = meshfilename.split(slash)[-1]
			self._parent.meshfilename = meshfilename
			self._filename = meshfilename
			self._parent.files.grid = meshfilename
			
		ln = infile.readline()
		N = int(infile.readline())
		for i in range(N):
			nd = infile.readline().strip().split()
			new_node = fnode(index=int(nd[0]),position=np.array([float(nd[1]),float(nd[2]),float(nd[3])]))
			self.add_node(new_node)	
		
		infile.readline()
		infile.readline()
		N = infile.readline()
		connectivity = int(N.strip().split()[0])
		N = int(N.strip().split()[1])
		for i in range(N):
			el = [int(eli) for eli in infile.readline().strip().split()]
			
			if self._full_connectivity:
				new_elem = felem(index = el[0], nodes = [self.node[eli] for eli in el[1:]])
				self.add_elem(new_elem)
				if connectivity == 8:
					nds1 = [el[1],el[1],el[1],el[7],el[7],el[7],el[4],el[4],el[6],el[6],el[2],el[5]]
					nds2 = [el[5],el[2],el[4],el[6],el[3],el[8],el[3],el[8],el[5],el[2],el[3],el[8]]
				elif connectivity == 4:
					nds1 = [el[1],el[2],el[3],el[4]]
					nds2 = [el[2],el[3],el[4],el[1]]
				elif connectivity == 3:
					nds1 = [el[1],el[2],el[3]]
					nds2 = [el[2],el[3],el[1]]
				else:
					print 'ERROR: unrecognized connectivity'; return
				for nd1,nd2 in zip(nds1,nds2):
					if nd1>nd2: ndi = nd2; nd2 = nd1; nd1 = ndi
					nd1 = self.node[nd1]; nd2 = self.node[nd2]
					new_conn = fconn(nodes = [nd1,nd2])
					
					nd1inds = []
					for con in nd1.connections:
						for nd in con.nodes: nd1inds.append(nd.index)
					
					if nd2.index in nd1inds: continue
					
					self.add_conn(new_conn)
					nd1.connections.append(new_conn)
					nd2.connections.append(new_conn)
				# associate element nodes with element
				for nd in self.elemlist[-1].nodes: 
					self._node[nd.index].elements.append(self.elemlist[-1])		
			else:
				self._elemlist.append(el[1:])
				self._elem[el[0]] = self.elemlist[-1]
		infile.close()		

		if ((len(np.unique([nd.position[0] for nd in self.nodelist])) == 1) or
			 (len(np.unique([nd.position[1] for nd in self.nodelist])) == 1) or
			 (len(np.unique([nd.position[2] for nd in self.nodelist])) == 1)):
			self._dimensions = 2
		else: self._dimensions = 3
	def write(self,filename='',format='fehm'):
		"""Write grid object to an FEHM grid file.

		:param filename: name of FEHM grid file to write to, including path specification, e.g. 'c:\\path\\file_out.inp'
		:type filename: str
		:param format: FEHM grid file format - currently options are 'fehm' (default) and 'avs'. 'avs' is the required format for reading into LaGrit.
		:type format: str

		"""
		if filename: self._filename=filename
		if self.filename=='': self._filename='fdata.grid'
		if slash in self.filename:
			make_directory(self.filename)
		if format == 'fehm': self._write_fehm(filename)
		elif format == 'avs': self._write_avs(filename)
	def _write_fehm(self,filename):
		if self._parent:
			if self._parent.work_dir and not os.path.isdir(self._parent.work_dir): os.system('mkdir '+self._parent.work_dir)
			outfile = open(self._parent.work_dir+self.filename,'w')
		else:
			outfile = open(self.filename,'w')
		outfile.write('coor\n')
		outfile.write('   '+str(len(self.nodelist))+'\n')		
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')			
		outfile.write('\t0\n')
		outfile.write('elem\n')
		if self._full_connectivity:
			outfile.write(str(len(self.elemlist[0].nodes))+' '+str(len(self.elemlist))+'\n')	
			for el in self.elemlist:
				outfile.write(str(int(el.index))+'   ')
				for nd in el.nodes:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')	
		else:
			outfile.write(str(len(self.elemlist[0]))+' '+str(len(self.elemlist))+'\n')	
			for i,el in enumerate(self.elemlist):
				outfile.write(str(i+1)+'   ')
				for nd in el:
					outfile.write(str(nd)+'   ')
				outfile.write('\n')		
			
		outfile.write('\nstop\n')
		outfile.close()
	def _write_avs(self,filename):
		if self._parent:
			if self._parent.work_dir and not os.path.isdir(self._parent.work_dir): os.system('mkdir '+self._parent.work_dir)
			outfile = open(self._parent.work_dir+self.filename,'w')
		else:
			outfile = open(self.filename,'w')
		outfile.write(str(self.number_nodes)+' '+str(self.number_elems)+' 0 0 0\n')
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')
			
		if self._full_connectivity:			
			for el in self.elemlist:
				outfile.write(str(int(el.index))+'   1  ')
				if len(el.nodes) == 8: outfile.write('hex  ')
				for nd in el.nodes:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')
		else:		
			for i,el in enumerate(self.elemlist):
				outfile.write(str(i+1)+'   1  ')
				if len(el) == 8: outfile.write('hex  ')
				for nd in el:
					outfile.write(str(nd)+'   ')
				outfile.write('\n')
		outfile.close()
	def make(self,meshfilename,x,y,z):
		""" Generates an orthogonal mesh for input node positions. 
		
		The mesh is constructed using the ``fgrid.``\ **fmake** object and an FEHM grid file is written for the mesh.
		
		:param meshfilename: Name to which to save the grid file.
		:type meshfilename: str
		:param x: Unique set of x-coordinates.
		:type x: list[fl64]
		:param y: Unique set of y-coordinates.
		:type y: list[fl64]
		:param z: Unique set of z-coordinates.
		:type z: list[fl64]
		"""
		if self._parent:
			if self._parent.work_dir and not os.path.isdir(self._parent.work_dir): os.system('mkdir '+self._parent.work_dir)
			fm = fmake(self._parent.work_dir+meshfilename,x,y,z,self._full_connectivity)
			fm.write()		
			self.read(self._parent.work_dir+meshfilename)
			self._parent._add_boundary_zones()
		else:
			fm = fmake(meshfilename,x,y,z)
			fm.write()		
			self.read(meshfilename)
	def add_node(self,node=fnode()):		#Add a node object.
		self._nodelist.append(node)
		self._node[node.index] = self._nodelist[-1]
	def add_conn(self,conn=fconn()):		#Add a connection object.
		self._connlist.append(conn)
		self._conn[(conn.nodes[0].index,conn.nodes[1].index)] = self.connlist[-1]
	def add_elem(self,elem=felem()):		#Add an element object.
		self._elemlist.append(elem)
		self._elem[elem.index] = self.elemlist[-1]
	def node_nearest_point(self,pos = []):
		"""Return node object nearest to position in space.

		:param pos: Coordinates, e.g. [2300., -134.8, 0.].
		:type pos: list
		:returns:  fnode() -- node object closest to position.

		"""
		min_dist = 1.e10
		nd = None
		if self.octree.leaf(pos) == None:
			print 'Error: point outside domain'
			return None
		lf = self.octree.leaf(pos)
		min_dist= 1.e10
		
		for leaf in ([lf,]+lf.neighbour):
			if leaf == None: continue
			dist,cube = leaf.min_dist(pos)
			if dist<min_dist: min_dist = dist; min_cube = cube
		return min_cube.elements[0]
	def plot(self,save='',angle=[45,45],color='k',connections=False,equal_axes=True,
		xlabel='x / m',ylabel='y / m',zlabel='z / m',title='',font_size='small',cutaway=[]): 		#generates a 3-D plot of the zone.
		"""Generates and saves a 3-D plot of the grid.
		
		:param save: Name of saved zone image.
		:type save: str
		:param angle: 	View angle of zone. First number is tilt angle in degrees, second number is azimuth. Alternatively, if angle is 'x', 'y', 'z', view is aligned along the corresponding axis.
		:type angle: [fl64,fl64], str
		:param color: Color of zone.
		:type color: str, [fl64,fl64,fl64]
		:param connections: Plot connections. If ``True`` all connections plotted. If between 0 and 1, random proportion plotted. If greater than 1, specified number plotted.
		:type connections: bool
		:param equal_axes: Force plotting with equal aspect ratios for all axes.
		:type equal_axes: bool
		
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param zlabel: Label on z-axis.
		:type zlabel: str
		:param title: Title of plot.
		:type title: str
		
		:param font_size: Size of text on plot.
		:type font_size: str, int
		
		:param cutaway: Coordinate from which cutaway begins. Alternatively, specifying 'middle','centre' with choose the centre of the grid as the cutaway point.
		:type cutaway: [fl64,fl64,fl64], str
		
		"""
		if not save: save = self.filename.split('.')[0]+'.png'
		if cutaway in ['middle','center','centre','mid']:			
			cutaway = [(self.xmin+self.xmax)/2,(self.ymin+self.ymax)/2,(self.zmin+self.zmax)/2]
		if isinstance(angle,str):
			if angle == 'x': angle = [0,0]
			elif angle == 'y': angle = [0,90]
			elif angle == 'z': angle = [90,90]
			else: return
			face1 = True; face2 = True; face3 = True; face4 = True; face5 = True; face6 = True
		else:
			while angle[0]<-90: angle[0]+=180
			while angle[0]>90: angle[0]-=180
			while angle[1]<0: angle[1]+=180
			while angle[1]>360: angle[1]-=180
			if angle[0]>0: face1 = True; face2 = False
			else: face1 = False; face2 = True
			if angle[1]>270 or angle[1]<=90: face3 = True; face4 = False
			else: face3 = False; face4 = True
			if angle[1]>0 and angle[1]<=180: face5 = True; face6 = False
			else: face5 = False; face6 = True
		# plot bounding box
		plt.clf()
		fig = plt.figure(figsize=[10.5,8.275])
		ax = plt.axes(projection='3d')
		ax.set_aspect('equal', 'datalim')
		
		ax.set_xlabel(xlabel,size=font_size)
		ax.set_ylabel(ylabel,size=font_size)
		ax.set_zlabel(zlabel,size=font_size)
		ax.set_title(title,size=font_size)
		
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_zticklabels():
			t.set_fontsize(font_size)

		xmin,xmax = self.xmin, self.xmax
		ymin,ymax = self.ymin, self.ymax
		zmin,zmax = self.zmin, self.zmax
		
		if equal_axes:
			MAX = np.max([xmax-xmin,ymax-ymin,zmax-zmin])/2
			C = np.array([xmin+xmax,ymin+ymax,zmin+zmax])/2
			for direction in (-1, 1):
				for point in np.diag(direction * MAX * np.array([1,1,1])):
					ax.plot([point[0]+C[0]], [point[1]+C[1]], [point[2]+C[2]], 'w')
		ax.view_init(angle[0],angle[1])
		
		if cutaway: xmid,ymid,zmid = cutaway
		else: 
			if face1:
				if face5:
					if face3: xmid,ymid,zmid = xmax,ymax,zmax
					else: xmid,ymid,zmid = xmin,ymax,zmax
				else:
					if face3:xmid,ymid,zmid = xmax,ymin,zmax
					else:xmid,ymid,zmid = xmin,ymin,zmax
			else:
				if face5:
					if face3: xmid,ymid,zmid = xmax,ymax,zmin
					else: xmid,ymid,zmid = xmin,ymax,zmin
				else:
					if face3:xmid,ymid,zmid = xmax,ymin,zmin
					else:xmid,ymid,zmid = xmin,ymin,zmin
					
		p13 = [xmid,ymid,zmid]
		if face1:
			if face5:
				if face3:					
					p1=[xmin,ymin,zmax]; p2=[xmin,ymax,zmax]; p3=[xmin,ymax,zmin]; p4=[xmax,ymax,zmin]
					p5=[xmax,ymin,zmin]; p6=[xmax,ymin,zmax]; p7=[xmax,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymax,zmax]; p10=[xmid,ymax,zmid]; p11=[xmax,ymax,zmid];p12=[xmax,ymid,zmid]
				else:						
					p1=[xmax,ymin,zmax]; p2=[xmax,ymax,zmax]; p3=[xmax,ymax,zmin]; p4=[xmin,ymax,zmin]
					p5=[xmin,ymin,zmin]; p6=[xmin,ymin,zmax]; p7=[xmin,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymax,zmax]; p10=[xmid,ymax,zmid]; p11=[xmin,ymax,zmid];p12=[xmin,ymid,zmid]
			else:
				if face3:						
					p1=[xmin,ymax,zmax]; p2=[xmin,ymin,zmax]; p3=[xmin,ymin,zmin]; p4=[xmax,ymin,zmin]
					p5=[xmax,ymax,zmin]; p6=[xmax,ymax,zmax]; p7=[xmax,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymin,zmax]; p10=[xmid,ymin,zmid]; p11=[xmax,ymin,zmid];p12=[xmax,ymid,zmid]
				else:		
					p1=[xmax,ymax,zmax]; p2=[xmax,ymin,zmax]; p3=[xmax,ymin,zmin]; p4=[xmin,ymin,zmin]
					p5=[xmin,ymax,zmin]; p6=[xmin,ymax,zmax]; p7=[xmin,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymin,zmax]; p10=[xmid,ymin,zmid]; p11=[xmin,ymin,zmid];p12=[xmin,ymid,zmid]
		else:
			if face5:
				if face3:					
					p1=[xmin,ymin,zmin]; p2=[xmin,ymax,zmin]; p3=[xmin,ymax,zmax]; p4=[xmax,ymax,zmax]
					p5=[xmax,ymin,zmax]; p6=[xmax,ymin,zmin]; p7=[xmax,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymax,zmin]; p10=[xmid,ymax,zmid]; p11=[xmax,ymax,zmid];p12=[xmax,ymid,zmid]
				else:						
					p1=[xmax,ymin,zmin]; p2=[xmax,ymax,zmin]; p3=[xmax,ymax,zmax]; p4=[xmin,ymax,zmax]
					p5=[xmin,ymin,zmax]; p6=[xmin,ymin,zmin]; p7=[xmin,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymax,zmin]; p10=[xmid,ymax,zmid]; p11=[xmin,ymax,zmid];p12=[xmin,ymid,zmid]
			else:
				if face3:						
					p1=[xmin,ymax,zmin]; p2=[xmin,ymin,zmin]; p3=[xmin,ymin,zmax]; p4=[xmax,ymin,zmax]
					p5=[xmax,ymax,zmax]; p6=[xmax,ymax,zmin]; p7=[xmax,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymin,zmin]; p10=[xmid,ymin,zmid]; p11=[xmax,ymin,zmid];p12=[xmax,ymid,zmid]
				else:		
					p1=[xmax,ymax,zmin]; p2=[xmax,ymin,zmin]; p3=[xmax,ymin,zmax]; p4=[xmin,ymin,zmax]
					p5=[xmin,ymax,zmax]; p6=[xmin,ymax,zmin]; p7=[xmin,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymin,zmin]; p10=[xmid,ymin,zmid]; p11=[xmin,ymin,zmid];p12=[xmin,ymid,zmid]
		pt1=p1;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p3;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p3;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p5;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p5;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p1;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
					
		# minor lines
		xs = np.unique([nd.position[0] for nd in self.nodelist])
		ys = np.unique([nd.position[1] for nd in self.nodelist])
		zs = np.unique([nd.position[2] for nd in self.nodelist])
		
		for x in xs:
			if x>=np.min([p2[0],p9[0]]) and x<=np.max([p2[0],p9[0]]):
				ax.plot([x,x],[p1[1],p2[1]],[p2[2],p2[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p2[1],p2[1]],[p2[2],p3[2]],color=color,linewidth=0.5)
			else:
				ax.plot([x,x],[p12[1],p2[1]],[p10[2],p10[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p12[1],p6[1]],[p2[2],p2[2]],color=color,linewidth=0.5)				
				ax.plot([x,x],[p7[1],p7[1]],[p7[2],p11[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p11[1],p11[1]],[p11[2],p3[2]],color=color,linewidth=0.5)
		for y in ys:
			if y>=np.min([p6[1],p7[1]]) and y<=np.max([p6[1],p7[1]]):
				ax.plot([p6[0],p1[0]],[y,y],[p2[2],p2[2]],color=color,linewidth=0.5)				
				ax.plot([p6[0],p6[0]],[y,y],[p2[2],p3[2]],color=color,linewidth=0.5)
			else:             
				ax.plot([p10[0],p11[0]],[y,y],[p10[2],p10[2]],color=color,linewidth=0.5)
				ax.plot([p9[0],p2[0]],[y,y],[p2[2],p2[2]],color=color,linewidth=0.5)							
				ax.plot([p4[0],p4[0]],[y,y],[p4[2],p11[2]],color=color,linewidth=0.5)
				ax.plot([p10[0],p10[0]],[y,y],[p10[2],p9[2]],color=color,linewidth=0.5)
		for z in zs:
			if z>=np.min([p4[2],p11[2]]) and z<=np.max([p4[2],p11[2]]):
				ax.plot([p4[0],p4[0]],[p5[1],p4[1]],[z,z],color=color,linewidth=0.5)				
				ax.plot([p4[0],p3[0]],[p4[1],p4[1]],[z,z],color=color,linewidth=0.5)
			else:
				ax.plot([p4[0],p4[0]],[p6[1],p7[1]],[z,z],color=color,linewidth=0.5)			
				ax.plot([p10[0],p10[0]],[p7[1],p11[1]],[z,z],color=color,linewidth=0.5)						
				ax.plot([p2[0],p8[0]],[p4[1],p4[1]],[z,z],color=color,linewidth=0.5)				
				ax.plot([p7[0],p8[0]],[p7[1],p7[1]],[z,z],color=color,linewidth=0.5)				
		
		extension, save_fname, pdf = save_name(save,variable='grid',time=1)
		if self._parent:
			if self._parent.work_dir and not os.path.isdir(self._parent.work_dir): os.system('mkdir '+self._parent.work_dir)
			plt.savefig(self._parent.work_dir+save_fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
			format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		else:
			plt.savefig(save_fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
			format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: os.system('epstopdf ' + save_fname); os.system(delStr+' ' + save_fname)			
	def add_nodetree(self):
		""" Constuct octree for node positions. Call to update if changes made to grid.
		"""
		self.octree=octree(self.bounding_box,self.nodelist)	
	def _summary(self):		
		L = 62
		print ''
		print ' ####---------------------------------------------------------####'
		line = ' #### FEHM grid file \''+self.filename+'\' summary.'
		for i in range(L-len(line)): line += ' '
		print line+'####'
		print ' ####---------------------------------------------------------####'
		
		lines = []
		lines.append(' #### Domain extent:')
		lines.append('x = ['+str(self.xmin) + ', ' + str(self.xmax) + ']')		
		lines.append('y = ['+str(self.ymin) + ', ' + str(self.ymax) + ']')
		lines.append('z = ['+str(self.zmin) + ', ' + str(self.zmax) + ']')
		lines.append(' #### Statistics:')
		lines.append('nodes = ' +str(self.number_nodes))
		lines.append('elements = ' +str(self.number_elems))
		
		for line in lines:
			if line.startswith(' ##'):
				for i in range(L-len(line)): line += ' '
				print line+'####'
			else:
				prntStr = ' #### -'
				keepGoing = True
				line = line.split()
				while keepGoing:
					if not line: 
						for i in range(L-len(prntStr)): prntStr += ' '
						print prntStr+'####'
						prntStr = ' #### '
						break
					if len(prntStr)<(L-len(line[0])): 
						prntStr += ' '+line[0]
						line = line[1:]
					else:
						for i in range(L-len(prntStr)): prntStr += ' '
						print prntStr+'####'
						prntStr = ' ####   '
		print ' ####---------------------------------------------------------####'
		print ''
	def rotate(self,angle=0.,centre=[0.,0.]):
		'''Rotates the grid by some angle about a specified vertical axis.
		
		:param angle: Clockwise angle by which to rotate grid.
		:type angle: fl64
		:param centre: x and y coordinates of vertical axis about which to rotate. Alternatively, the centre of the computational domain can be specified by passing 'mid','middle','centre', or 'center'.
		:type centre: [fl64,fl64], str
		'''
		if centre in ['middle','mid','centre','center']:
			centre = [(self.xmin+self.xmax)/2.,(self.ymin+self.ymax)/2.]
		for nd in self.nodelist:
			old_pos = np.array(nd.position[0:2]) - np.array(centre) 	# position relative to centre of rotation
			theta_f = math.atan2(old_pos[1],old_pos[0]) + angle/180.*math.pi
			dist = np.sqrt(np.dot(old_pos,old_pos))
			new_pos = [dist*math.cos(theta_f),dist*math.sin(theta_f)]
			nd.position[0] = new_pos[0]+centre[0]
			nd.position[1] = new_pos[1]+centre[1]
	def get_bounding_box(self):
		minx = self.nodelist[0].position[0]
		maxx = self.nodelist[0].position[0]
		miny = self.nodelist[0].position[1]
		maxy = self.nodelist[0].position[1]
		minz = self.nodelist[0].position[2]
		maxz = self.nodelist[0].position[2]
		for node in self.nodelist:
			pos = node.position
			if pos[0]<minx: minx = pos[0]
			if pos[1]<miny: miny = pos[1]
			if pos[2]<minz: minz = pos[2]
			if pos[0]>maxx: maxx = pos[0]
			if pos[1]>maxy: maxy = pos[1]
			if pos[2]>maxz: maxz = pos[2]
		xr = maxx-minx+1.; yr = maxy-miny+1.; zr = maxz-minz+1.;
		
		c = np.array([minx+maxx,miny+maxy,minz+maxz])/2.
		
		r = np.max([xr,yr,zr])
		
		if minx == maxx:
			c[0] = minx; r = np.array([1,r,r])
			return [c-0.51*r,c+0.51*r]
		elif miny == maxy:
			c[1] = miny; r = np.array([r,1,r])
			return [c-0.51*r,c+0.51*r]
		elif minz == maxz:
			c[2] = minz; r = np.array([r,r,1])
			return [c-0.51*r,c+0.51*r]
		else:			
			return [c-0.51*r,c+0.51*r]
	bounding_box = property(get_bounding_box)	
	def _get_filename(self): return self._filename
	filename = property(_get_filename)#: (*str*) Name of FEHM grid file.
	def _get_node(self): return self._node
	node = property(_get_node)#: (*dict[fnode]*) Dictionary of grid nodes, indexed by node integer.
	def _get_nodelist(self): return self._nodelist
	nodelist = property(_get_nodelist)#: (*lst[fnode]*) List of all node objects in the grid.
	def _get_elem(self): return self._elem
	elem = property(_get_elem)#: (*dict[felem]*) Dictionary of elements, indexed by element integer.
	def _get_elemlist(self): return self._elemlist
	elemlist = property(_get_elemlist)#: (*lst[felem]*) List of all element objects in the grid.
	def _get_conn(self): return self._conn
	conn = property(_get_conn)#: (*dict[fconn]*) Dictionary of connections, indexed by a two element tuple of the member node integers.
	def _get_connlist(self): return self._connlist
	connlist = property(_get_connlist)#: (*lst[fconn]*) List of all connection objects in the grid.
	def _get_dimensions(self): return self._dimensions
	dimensions = property(_get_dimensions) #: (*int*) Dimensions of the grid.
	def get_xmin(self): return np.min([nd.position[0] for nd in self.nodelist])
	xmin = property(get_xmin) 				#: Minimum x-coordinate for all nodes.
	def get_xmax(self): return np.max([nd.position[0] for nd in self.nodelist])
	xmax = property(get_xmax)				#: Maximum x-coordinate for all nodes.
	def get_ymin(self): return np.min([nd.position[1] for nd in self.nodelist])
	ymin = property(get_ymin)				#: Minimum y-coordinate for all nodes.
	def get_ymax(self): return np.max([nd.position[1] for nd in self.nodelist])
	ymax = property(get_ymax)				#: Maximum y-coordinate for all nodes.
	def get_zmin(self): return np.min([nd.position[2] for nd in self.nodelist])
	zmin = property(get_zmin)				#: Minimum z-coordinate for all nodes.
	def get_zmax(self): return np.max([nd.position[2] for nd in self.nodelist])
	zmax = property(get_zmax)				#: Maximum z-coordinate for all nodes.
	def get_node_number(self): return len(self.nodelist)
	number_nodes = property(get_node_number)#: Number of nodes in grid.
	def get_element_number(self): return len(self.elemlist)
	number_elems = property(get_element_number)#: Number of elements in grid.
	def get_info(self):
		print 'FEHM grid file \''+self.filename+'\' summary.'
		print 'Model domain: x = ['+str(self.xmin) + ', ' + str(self.xmax) + ']'
		print '              y = ['+str(self.ymin) + ', ' + str(self.ymax) + ']'
		print '              z = ['+str(self.zmin) + ', ' + str(self.zmax) + ']'
		print '          nodes = ' +str(self.number_nodes)
		print '       elements = ' +str(self.number_elems)
		print ' '
	what = property(get_info) 				#: Print to screen information about the grid.
class fmake(object): 				#Rectilinear grid constructor.
	"""Generate an orthogonal mesh corresponding to vectors of nodal positions.
	"""
	def __init__(self,meshname,x=None,y=None,z=None,full_connectivity=False):
		self._x = list(np.unique(x))
		self._y = list(np.unique(y))
		self._z = list(np.unique(z))
		self._dimension = None
		self._meshname = ''
		self._full_connectivity = full_connectivity
		if meshname: self._meshname = meshname
	def seed(self,edge='x',method='equal', number = None, size = None, bias = None):
		"""Allocate mesh seeds to edges. NOT FINISHED (not started...)
		"""
		
	def write(self,meshname=''):
		"""Write out the grid file.
		
		:param meshname: Name of the grid file.
		:type meshname: str
		"""
		self.refresh()
		if meshname: self._meshname = meshname
		if self.meshname=='': self._meshname='fdata.grid'
		outfile = open(self.meshname,'w')
		outfile.write('coor\n')
		outfile.write('   '+str(len(self.nodelist))+'\n')
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')
		outfile.write('\t0\n')
		outfile.write('elem\n')
		if self._full_connectivity:
			outfile.write(str(len(self.elemlist[0].nodes))+' '+str(len(self.elemlist))+'\n')		
			for el in self.elemlist:
				outfile.write(str(int(el.index))+'   ')
				for nd in el.nodes:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')
		else:
			outfile.write(str(len(self.elemlist[0]))+' '+str(len(self.elemlist))+'\n')		
			for i,el in enumerate(self.elemlist):
				outfile.write(str(i+1)+'   ')
				for nd in el:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')
		outfile.write('\nstop\n')
		outfile.close()
	def refresh(self):
		"""Generate grid node and element objects corresponding to x y and z seeds.
		"""
		# first determine dimension of grid
		xF = self.x != None; yF = self.y != None; zF = self.z != None
		if xF and yF and zF: self.dimension = 3
		elif (xF and yF and not zF) or (xF and zF and not yF) or (zF and yF and not xF): self.dimension = 2
		else: print 'ERROR: not enough dimensions specified'; return
		if self.dimension == 2: print 'ERROR: two dimensional grids not supported'; return
		if self.dimension == 3:
			# create nodes
			self._nodelist = []
			ind = 1
			self._z = list(np.sort(self._z))
			self._x = list(np.sort(self._x))
			self._y = list(np.sort(self._y))
			for zi in self._z:
				for yi in self._y:
					for xi in self._x:
						self._nodelist.append(fnode(index=ind,position=[xi,yi,zi]))
						ind +=1
			
			# create elements
			self._elemlist = []
			ind = 1
			xL = len(self._x); yL = len(self._y); zL = len(self._z)
			for i in range(1,len(self._z)):
				for j in range(1,len(self._y)):
					for k in range(1,len(self._x)):
						nodes=[
							self._nodelist[i*xL*yL + (j-1)*xL + k-1],
							self._nodelist[i*xL*yL + (j-1)*xL + k],
							self._nodelist[i*xL*yL + j*xL + k],
							self._nodelist[i*xL*yL + j*xL + k-1],
							self._nodelist[(i-1)*xL*yL + (j-1)*xL + k-1],
							self._nodelist[(i-1)*xL*yL + (j-1)*xL + k],
							self._nodelist[(i-1)*xL*yL + j*xL + k],
							self._nodelist[(i-1)*xL*yL + j*xL + k-1],
							]
						if self._full_connectivity:
							self._elemlist.append(felem(index=ind,nodes=nodes))
							ind +=1
						else:
							self._elemlist.append(nodes)
	def _get_x(self): return self._x
	def _set_x(self,value): self._x = value
	x = property(_get_x, _set_x) #: (*lst[fl64]*) x coordinates of nodes.
	def _get_y(self): return self._y
	def _set_y(self,value): self._y = value
	y = property(_get_y, _set_y) #: (*lst[fl64]*) y coordinates of nodes.
	def _get_z(self): return self._z
	def _set_z(self,value): self._z = value
	z = property(_get_z, _set_z) #: (*lst[fl64]*) z coordinates of nodes.
	def _get_meshname(self): return self._meshname
	def _set_meshname(self,value): self._meshname = value
	meshname = property(_get_meshname, _set_meshname) #: (*str*) Name of grid file to write out.
	def _get_nodelist(self): return self._nodelist
	def _set_nodelist(self,value): self._nodelist = value
	nodelist = property(_get_nodelist, _set_nodelist) #: (*lst[fnode]*) List of node objects in the grid.
	def _get_elemlist(self): return self._elemlist
	def _set_elemlist(self,value): self._elemlist = value
	elemlist = property(_get_elemlist, _set_elemlist) #: (*lst[felem]*) List of element objects in the grid.
