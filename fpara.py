"""Classes for PyFEHM VTK output"""

"""
Copyright 2013.
Los Alamos National Security, LLC. 
This material was produced under U.S. Government contract DE-AC52-06NA25396 for 
Los Alamos National Laboratory (LANL), which is operated by Los Alamos National 
Security, LLC for the U.S. Department of Energy. The U.S. Government has rights 
to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES 
ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce 
derivative works, such modified software should be clearly marked, so as not to 
confuse it with the version available from LANL.

Additionally, this library is free software; you can redistribute it and/or modify 
it under the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation; either version 2.1 of the License, or (at your option) 
any later version. Accordingly, this library is distributed in the hope that it 
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General 
Public License for more details.
"""

import numpy as np
from ftool import*
import pyvtk as pv
from copy import copy

class fVtkData(pv.VtkData):
	def __init__(self,*args,**kws):
		pv.VtkData.__init__(self,*args,**kws)
		self.times = []
		self.material = pv.PointData()
		self.contour = {}
	def to_string(self, time=None, format = 'ascii',material=False):
		ret = ['# vtk DataFile Version 2.0',
			   self.header,
			   format.upper(),
			   self.structure.to_string(format)
			   ]
		if self.cell_data.data:
			ret.append(self.cell_data.to_string(format))
		if material:
			ret.append(self.material.to_string(format))
		else:
			if self.contour[time].data:
				ret.append(self.contour[time].to_string(format))
		return '\n'.join(ret)
	def tofile(self, filename, format = 'ascii'):
		"""Save VTK data to file.
		"""
		written_files = []
		if not pv.common.is_string(filename):
			raise TypeError,'argument filename must be string but got %s'%(type(filename))
		if format not in ['ascii','binary']:
			raise TypeError,'argument format must be ascii | binary'
		filename = filename.strip()
		if not filename:
			raise ValueError,'filename must be non-empty string'
		if filename[-4:]!='.vtk':
			filename += '.vtk'
		
		# first write material properties file
		filename_int = ''.join(filename[:-4]+'_mat.vtk')
		f = open(filename_int,'wb')
		f.write(self.to_string(format,material=True))
		f.close()
		written_files.append(filename_int)
		
		# write contour output file
		times = np.sort(self.contour.keys())
		for i,time in enumerate(times):
			if len(times)>1:
				filename_int = ''.join(filename[:-4]+'.%04i'%i+'.vtk')
			else:
				filename_int = filename
			#print 'Creating file',`filename`
			f = open(filename_int,'wb')
			f.write(self.to_string(time,format))
			f.close()
			written_files.append(filename_int)
		return written_files
class fvtk(object):
	def __init__(self,parent,filename,contour,show_zones,diff):
		self.parent = parent
		self.path = fpath(parent = self)
		self.path.filename = filename
		self.data = None
		self.contour = contour
		self.variables = []
		self.materials = []
		self.zones = []
		self.show_zones = show_zones
		self.diff = diff
	def assemble(self):		
		"""Assemble all information in pyvtk objects."""			
		self.assemble_grid()		# add grid information
		self.assemble_zones()		# add zone information
		self.assemble_properties()	# add permeability data
		if self.contour != None:	# add contour data
			self.assemble_contour()
	def assemble_grid(self):
		"""Assemble grid information in pyvtk objects."""
		# node positions, connectivity information
		nds = [nd.position for nd in self.parent.grid.nodelist]		
		cns = [[nd.index-1 for nd in el.nodes] for el in self.parent.grid.elemlist]
		
		# make grid
		self.data = fVtkData(pv.UnstructuredGrid(nds,hexahedron=cns),'PyFEHM VTK model output')
		
		# grid information
		dat = np.array([nd.position for nd in self.parent.grid.nodelist])
		nds = np.array([nd.index for nd in self.parent.grid.nodelist])
		self.data.material.append(pv.Scalars(nds,name='n',lookup_table='default'))
		self.data.material.append(pv.Scalars(dat[:,0],name='x',lookup_table='default'))
		self.data.material.append(pv.Scalars(dat[:,1],name='y',lookup_table='default'))
		self.data.material.append(pv.Scalars(dat[:,2],name='z',lookup_table='default'))
		
		self.x_lim = [np.min(dat[:,0]),np.max(dat[:,0])]
		self.y_lim = [np.min(dat[:,1]),np.max(dat[:,1])]
		self.z_lim = [np.min(dat[:,2]),np.max(dat[:,2])]
		self.n_lim = [1,len(self.parent.grid.nodelist)]
	def assemble_zones(self):
		"""Assemble zone information in pyvtk objects."""
		# zones will be considered material properties as they only need to appear once
		N = len(self.parent.grid.nodelist)
		nds = np.zeros((1,N))[0]
		self.parent.zonelist.sort(key=lambda x: x.index)
		for zn in self.parent.zonelist:
			if zn.index == 0: continue
			name = 'zone%04i'%zn.index
			if zn.name: name += '_'+zn.name.replace(' ','_')
			self.zones.append(name)
			zn_nds = copy(nds)
			for nd in zn.nodelist: zn_nds[nd.index-1] = 1
			self.data.material.append(
				pv.Scalars(zn_nds,
				name=name,
				lookup_table='default'))
	def assemble_properties(self):
		"""Assemble material properties in pyvtk objects."""
		# permeabilities
		perms = np.array([nd.permeability for nd in self.parent.grid.nodelist])
		if np.mean(perms)>0.: perms = np.log10(perms)
		
		self.add_material('perm_x',perms[:,0])
		self.add_material('perm_y',perms[:,1])
		self.add_material('perm_z',perms[:,2])

		props = np.array([[nd.density, nd.porosity, nd.specific_heat, nd.youngs_modulus,nd.poissons_ratio,nd.thermal_expansion,nd.pressure_coupling,nd.Ti,nd.Pi,nd.Si] 	for nd in self.parent.grid.nodelist])
		names = ['density','porosity','specific_heat','youngs_modulus','poissons_ratio','thermal_expansion','pressure_coupling','Pi','Ti','Si']
		for name, column in zip(names,props.T):
			self.add_material(name,column)
	def add_material(self,name,data):
		if all(v is None for v in data): return 		# if all None, no data to include
		data = np.array([dt if dt != None else -1.e30 for dt in data]) 		# check for None, replace with -1.e30
		self.data.material.append(pv.Scalars(data,name=name,lookup_table='default'))
		self.materials.append(name)
		self.__setattr__(name+'_lim',[np.min(data),np.max(data)])
	def assemble_contour(self):
		"""Assemble contour output in pyvtk objects."""
		self.data.contour = dict([(time,pv.PointData()) for time in self.contour.times])
		if self.diff: time0 = self.contour.times[0]
		for time in self.contour.times:
			do_lims = (time == self.contour.times[-1])
			for var in self.contour.variables:
				if time != self.contour.times[0] and var in ['x','y','z','n']: continue
				if var not in self.variables: self.variables.append(var)
				self.data.contour[time].append(pv.Scalars(self.contour[time][var],name=var,lookup_table='default'))
				if self.diff:
					self.data.contour[time].append(pv.Scalars(self.contour[time][var]-self.contour[time0][var],name='diff_'+var,lookup_table='default'))
				if do_lims: self.__setattr__(var+'_lim',[np.min(self.contour[time][var]),np.max(self.contour[time][var])])
	def write(self):	
		"""Call to write out vtk files."""
		if self.parent.work_dir: wd = self.parent.work_dir
		else: wd = self.parent._path.absolute_to_file
		fls = self.data.tofile(wd+slash+self.path.filename)
		# save file names for later use
		self.material_file = fls[0]
		self.contour_files = []
		if len(fls)>1:
			self.contour_files = fls[1:]
		return fls
	def initial_display(self,show):
		"""Determines what variable should be initially displayed."""
		mat_vars = ['n','x','y','z','perm_x','perm_y','perm_z','porosity','density','cond_x','cond_y','cond_z']
		if self.contour:
			cont_vars = self.contour.variables
		
		# convert k* format to perm_*
		if show == 'kx': show = 'perm_x'
		elif show == 'ky': show = 'perm_y'
		elif show == 'kz': show = 'perm_z'
		
		# check for unspecified coordinate in potentially anisotropic properties
		if show in ['permeability','perm']:
			print 'NOTE: plotting z-component of permeability, for other components specify show=\'perm_x\', etc.'
			show = 'perm_z'
		if show in ['conducitivity','cond']:
			print 'NOTE: plotting z-component of conductivity, for other components specify show=\'cond_x\', etc.'
			show = 'cond_z'
		
		# check if material property or contour output requested for display
		if show in mat_vars:
			self.initial_show = 'material'
			self.default_material_property = show
			self.default_material_lims = self.__getattribute__(show+'_lim')
			if self.contour:
				# get default contour variable to display
				for var in self.contour.variables:
					if var not in ['x','y','z','n']: break
				self.default_contour_variable = var
				self.default_contour_lims = self.__getattribute__(var+'_lim')
		elif show in cont_vars:
			self.initial_show = 'contour'
			self.default_contour_variable = show
			self.default_contour_lims = self.__getattribute__(show+'_lim')
			self.default_material_property = 'perm_x' 		# default 
			self.default_material_lims = self.__getattribute__('perm_x_lim')
		else:
			print 'ERROR: requested property or variable does not exist, available options are...'
			print 'Material properties:'
			for mat in mat_vars:
				print '  - '+mat
			print 'Contour output variables:'
			for var in cont_vars:
				print '  - '+var
			print ''
	def startup_script(self):
		x0,x1 = self.parent.grid.xmin, self.parent.grid.xmax
		y0,y1 = self.parent.grid.ymin, self.parent.grid.ymax
		z0,z1 = self.parent.grid.zmin, self.parent.grid.zmax
		xm,ym,zm = (x0+x1)/2., (y0+y1)/2., (z0+z1)/2.
		xr,yr,zr = (x1-x0), (y1-y0), (z1-z0)
		
		dflt_mat = '\''+self.default_material_property+'\''		
		mat_lim = self.default_material_lims
		
		f = open('pyfehm_paraview_startup.py','w')
		
		contour_files=[file for file in self.contour_files]
		
		################################### load paraview modules ######################################
		lns = [
			'try: paraview.simple',
			'except: from paraview.simple import *',
			'paraview.simple._DisableFirstRenderCameraReset()',
			'',
			]
			
		################################### load material properties ###################################
		lns += ['mat_prop = LegacyVTKReader( FileNames=[']
		file = self.material_file.replace('\\','/')
		lns += ['\''+file+'\',']
		lns += ['] )']
		lns += ['RenameSource("model", mat_prop)']
			
		################################### initial property display ###################################
		lns += [
			'rv = GetRenderView()',
			'dr = Show()',
			'dr.ScalarOpacityUnitDistance = 1.7320508075688779',
			'dr.EdgeColor = [0.0, 0.0, 0.5]',
			'',
			'rv.CenterOfRotation = [%10.5f, %10.5f, %10.5f]'%(xm,ym,zm),
			'',
			'rv.CameraViewUp = [-0.4, -0.11, 0.92]',
			'rv.CameraPosition = [%10.5f, %10.5f, %10.5f]'%(xm+2.5*xr,ym+1.5*yr,zm+1.5*zr),
			'rv.CameraFocalPoint = [%10.5f, %10.5f, %10.5f]'%(xm,ym,zm),
			'',
			'mr = GetDisplayProperties(mat_prop)',
			'mr.Representation = \'Surface With Edges\'',
			'',
			'lt = GetLookupTableForArray( '+dflt_mat+', 1, RGBPoints=[%4.2f, 0.23, 0.299, 0.754, %4.2f, 0.706, 0.016, 0.15], VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )'%tuple(mat_lim),
			'',
			'pf = CreatePiecewiseFunction( Points=[%4.2f, 0.0, 0.5, 0.0, %4.2f, 1.0, 0.5, 0.0] )'%tuple(mat_lim),
			'',
			'mr.ScalarOpacityFunction = pf',
			'mr.ColorArrayName = (\'POINT_DATA\', '+dflt_mat+')',
			'mr.LookupTable = lt',
			'',
			'lt.ScalarOpacityFunction = pf',
			'',
			'ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='+dflt_mat+', LabelFontSize=12, Enabled=1, TitleFontSize=12 )',
			'GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)',
			'',
			'lt = GetLookupTableForArray('+dflt_mat+', 1 )',
			'',
			'ScalarBarWidgetRepresentation1.LookupTable = lt',
			'',
			]
			
		################################### load in nodes as glyphs ###################################	
		
		ndRadius = np.min([con.distance for con in self.parent.grid.connlist])/10.
		lns += [
			'AnimationScene1 = GetAnimationScene()',
			'AnimationScene1.AnimationTime = 0.0',
			'rv.ViewTime = 0.0',
			'source = FindSource("model")',
			'SetActiveSource(source)',
			'',
			'G = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )',
			'G.GlyphTransform = "Transform2"',
			'G.GlyphType = "Sphere"',
			'G.RandomMode = 0',
			'G.ScaleMode = \'off\'',
			'G.MaskPoints = 0',
			'G.GlyphType.Radius = %10.5f'%ndRadius,
			'',
			'RenameSource("nodes", G)',
			'',
			'rv = GetRenderView()',
			'mr = GetDisplayProperties(source)',
			'dr = Show()',
			'dr.ColorArrayName = (\'POINT_DATA\', \'n\')',
			'dr.ScaleFactor = 1.1',
			'dr.SelectionPointFieldDataArrayName = "nodes"',
			'dr.EdgeColor = [0.0, 0.0, 0.5000076295109483]',
			'dr.ColorArrayName = (\'POINT_DATA\', \'\')',
			'dr.DiffuseColor = [0.,0.,0.]',
			'dr.Visibility = 0',
			]
		################################### load in zones as glyphs ###################################		
		colors = [
			[1.,1.,0.],
			[1.,0.,1.],
			[0.,1.,1.],
			[1.,1.,0.5],
			[1.,0.5,1.],
			[0.5,1.,1.],
			[1.,1.,0.25],
			[1.,0.25,1.],
			[0.25,1.,1.],
			[1.,1.,0.75],
			[1.,0.75,1.],
			[0.75,1.,1.],
			
			[0.5,1.,0.5],
			[1.,0.5,0.5],
			[0.5,0.5,1.],			
			[0.5,0.75,0.5],
			[0.75,0.5,0.5],
			[0.5,0.5,0.75],		
			[0.5,0.25,0.5],
			[0.25,0.5,0.5],
			[0.5,0.5,0.25],
			
			[0.75,1.,0.75],
			[1.,0.75,0.75],
			[0.75,0.75,1.],			
			[0.75,0.5,0.75],
			[0.5,0.75,0.75],
			[0.75,0.75,0.5],			
			[0.75,0.25,0.75],
			[0.25,0.75,0.75],
			[0.75,0.75,0.25],			
			
			[0.25,1.,0.25],
			[1.,0.25,0.25],
			[0.25,0.25,1.],			
			[0.25,0.75,0.25],
			[0.75,0.25,0.25],
			[0.25,0.25,0.75],			
			[0.25,0.5,0.25],
			[0.5,0.25,0.25],
			[0.25,0.25,0.5],
			]
		
		zones = []; cols = []
		for zone,color in zip(self.zones,colors):
			if self.show_zones == 'user':			
				if ('XMIN' in zone) or ('XMAX' in zone) or ('YMIN' in zone) or ('YMAX' in zone) or ('ZMIN' in zone) or ('ZMAX' in zone): continue
			zones.append(zone)
			cols.append(color)
		
		lns += ['cols = [']
		for col in cols:
			lns += ['[%3.2f,%3.2f,%3.2f],'%tuple(col)]
		lns += [']']
		lns += ['zones = [']
		for zone in zones:
			lns += ['\''+zone+'\',']
		lns += [']']
		lns += ['for zone,col in zip(zones,cols):',
			'\tAnimationScene1 = GetAnimationScene()',
			'\tAnimationScene1.AnimationTime = 0.0',
			'\trv.ViewTime = 0.0',
			'\tsource = FindSource("model")',
			'\tSetActiveSource(source)',
			'\t',
			'\tG = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )',
			'\tG.GlyphTransform = "Transform2"',
			'\tG.Scalars = [\'POINTS\', zone]',
			'\tG.ScaleMode = \'scalar\'',
			'\tG.GlyphType = "Sphere"',
			'\tG.RandomMode = 0',
			'\tG.MaskPoints = 0',
			'\t',
			'\tG.GlyphType.Radius = %10.5f'%(2*ndRadius),
			'\t',
			'\tRenameSource(zone, G)',
			'\t',
			'\trv = GetRenderView()',
			'\tmr = GetDisplayProperties(source)',
			'\tdr = Show()',
			'\tdr.ColorArrayName = (\'POINT_DATA\', \'n\')',
			'\tdr.ScaleFactor = 1.1',
			'\tdr.SelectionPointFieldDataArrayName = zone',
			'\tdr.EdgeColor = [0.0, 0.0, 0.5000076295109483]',
			'\tdr.Opacity = 0.5',
			'\t',
			'\tlt = GetLookupTableForArray(zone, 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.5, 0.865, 0.865, 0.865, 1.0]+col, VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )',
			'\t',
			'\tpf = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )',
			'\t',
			'\tdr.ColorArrayName = (\'POINT_DATA\', zone)',
			'\tdr.LookupTable = lt',
			'\tdr.Visibility = 0',
			'\t',
			'\tlt.ScalarOpacityFunction = pf',
			]
		
		
		################################### load in contour output ###################################	
		if len(contour_files)>0:
			lns += ['contour_output = LegacyVTKReader( FileNames=[']
			for file in contour_files:
				file = file.replace('\\','/')
				lns += ['\''+file+'\',']
			lns += ['] )']
			lns += ['RenameSource("contour_output", contour_output)']
		
		################################### set up initial visualisation ###################################	
			dflt_cont = '\''+self.default_contour_variable+'\''		
			cont_lim = self.default_contour_lims
		
			viewTime = len(self.contour.times)-1
		
			lns += [
				'lt = GetLookupTableForArray('+dflt_cont+', 1, RGBPoints=[%10.5f, 0.23, 0.299, 0.754, %10.5f, 0.706, 0.016, 0.15], VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )'%tuple(cont_lim),
				'',
				'pf = CreatePiecewiseFunction( Points=[%10.5f, 0.0, 0.5, 0.0, %10.5f, 1.0, 0.5, 0.0] )'%tuple(cont_lim),
				'',
				'dr = Show() #dr = DataRepresentation1',
				'dr.Representation = \'Surface With Edges\'',
				'dr.EdgeColor = [0.15, 0.15, 0.15]',
				'dr.ScalarOpacityFunction = pf',
				'dr.ColorArrayName = (\'POINT_DATA\', '+dflt_cont+')',
				'dr.ScalarOpacityUnitDistance = 1.7320508075688779',
				'dr.LookupTable = lt',
				'',
				'rv.ViewTime = %4i'%viewTime,
				'',
				'ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='+dflt_cont+', LabelFontSize=12, Enabled=1, LookupTable=lt, TitleFontSize=12 )',
				'GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)',
				'',
				]
			
		if len(contour_files)>0:
			lns+= [
				'model = FindSource("model")',
				'model_rep = GetDisplayProperties(model)',
				'contour_output = FindSource("contour_output")',
				'cont_rep = GetDisplayProperties(contour_output)',
				]	
			if self.initial_show == 'material':
				lns+=[
					'model_rep.Visibility = 1',
					'cont_rep.Visibility = 0	',
					]					
			elif self.initial_show == 'contour':
				lns += [
					'model_rep.Visibility = 0',
					'cont_rep.Visibility = 1',
					]			
		f.writelines('\n'.join(lns))
		f.close()
	def _get_filename(self): return self.path.absolute_to_file+slash+self.path.filename
	filename = property(_get_filename) #: (**)