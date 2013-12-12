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

class fVtkData(pv.VtkData):
	def __init__(self,*args,**kws):
		pv.VtkData.__init__(self,*args,**kws)
		self.times = []
		self.material = pv.PointData()
		self.contour = {}
		
	def to_string(self, time, format = 'ascii',material=False):
		ret = ['# vtk DataFile Version 2.0',
			   self.header,
			   format.upper(),
			   self.structure.to_string(format)
			   ]
		if self.cell_data.data:
			ret.append(self.cell_data.to_string(format))
		if material:
			for data in self.contour[time].data:
				self.material.append(data)
			ret.append(self.material.to_string(format))
		else:
			if self.contour[time].data:
				ret.append(self.contour[time].to_string(format))
		return '\n'.join(ret)

	def tofile(self, filename, format = 'ascii'):
		"""Save VTK data to file.
		"""
		written_files = []
		times = np.sort(self.contour.keys())
		for i,time in enumerate(times):
			if not pv.common.is_string(filename):
				raise TypeError,'argument filename must be string but got %s'%(type(filename))
			if format not in ['ascii','binary']:
				raise TypeError,'argument format must be ascii | binary'
			filename = filename.strip()
			if not filename:
				raise ValueError,'filename must be non-empty string'
			if filename[-4:]!='.vtk':
				filename += '.vtk'
			if len(times)>1:
				filename_int = ''.join(filename[:-4]+'.%04i'%i+'.vtk')
			else:
				filename_int = filename
			#print 'Creating file',`filename`
			f = open(filename_int,'wb')
			if i == 0:
				f.write(self.to_string(time,format,material=True))
			else:
				f.write(self.to_string(time,format))
			f.close()
			written_files.append(filename_int)
		return written_files
		
class fvtk(object):
	def __init__(self,parent,filename,contour):
		self.parent = parent
		self.path = fpath(parent = self)
		self.path.filename = filename
		self.data = None
		self.contour = contour
		self.variables = []
	def assemble(self):			
		# node positions, connectivity information
		nds = [nd.position for nd in self.parent.grid.nodelist]		
		cns = [[nd.index-1 for nd in el.nodes] for el in self.parent.grid.elemlist]
		
		# make grid
		self.data = fVtkData(pv.UnstructuredGrid(nds,hexahedron=cns))
		
		# add permeability data
		self.assemble_properties()
		
		# add contour data
		if self.contour != None:
			self.assemble_contour()
	
	def assemble_properties(self):
		nan = float('NaN')
		
		perms = np.array([nd.permeability for nd in self.parent.grid.nodelist])
		if np.mean(perms)>0.: perms = np.log10(perms)
		
		self.kx_lim = [np.min(perms[:,0]),np.max(perms[:,0])]
		self.ky_lim = [np.min(perms[:,1]),np.max(perms[:,1])]
		self.kz_lim = [np.min(perms[:,2]),np.max(perms[:,2])]
		
		self.data.material.append(pv.Scalars(perms[:,0],name='kx',lookup_table='default'))
		self.data.material.append(pv.Scalars(perms[:,1],name='ky',lookup_table='default'))
		self.data.material.append(pv.Scalars(perms[:,2],name='kz',lookup_table='default'))
		self.variables.append('kx')
		self.variables.append('ky')
		self.variables.append('kz')
		
		dens = np.array([nd.density if nd.density != nan else nan for nd in self.parent.grid.nodelist])
		self.data.material.append(pv.Scalars(dens,name='density',lookup_table='default'))
		self.variables.append('density')
		
	def assemble_contour(self):
		self.data.contour = dict([(time,pv.PointData()) for time in self.contour.times])
		for time in self.contour.times:
			do_lims = (time == self.contour.times[-1])
			for var in self.contour.variables:
				if time != self.contour.times[0] and var in ['x','y','z','n']: continue
				if var not in self.variables: self.variables.append(var)
				self.data.contour[time].append(pv.Scalars(self.contour[time][var],name=var,lookup_table='default'))
				
				if do_lims: self.__setattr__(var+'_lim',[np.min(self.contour[time][var]),np.max(self.contour[time][var])])
	def write(self):	
		if self.parent.work_dir: wd = self.parent.work_dir
		else: wd = self.parent._path.absolute_to_file
		return self.data.tofile(wd+slash+self.path.filename)
	def startup_script(self,show):
		x0,x1 = self.parent.grid.xmin, self.parent.grid.xmax
		y0,y1 = self.parent.grid.ymin, self.parent.grid.ymax
		z0,z1 = self.parent.grid.zmin, self.parent.grid.zmax
		xm,ym,zm = (x0+x1)/2., (y0+y1)/2., (z0+z1)/2.
		xr,yr,zr = (x1-x0), (y1-y0), (z1-z0)
		initial_display = '\''+show+'\''
		
		lim = self.__getattribute__(show+'_lim')
		
		f = open('pyfehm_paraview_startup.py','w')
		viewTime = len(self.contour.times)-1
		lns = [
			'try: paraview.simple',
			'except: from paraview.simple import *',
			'paraview.simple._DisableFirstRenderCameraReset()',
			'',
			'temp_vtk = GetActiveSource()',
			'RenderView1 = GetRenderView()',
			'DataRepresentation1 = Show()',
			'DataRepresentation1.ScalarOpacityUnitDistance = 1.7320508075688779',
			'DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5]',
			'',
			'Render()',
			'RenderView1.CenterOfRotation = [%10.5f, %10.5f, %10.5f]'%(xm,ym,zm),
			'',
			'RenderView1.CameraViewUp = [-0.4, -0.11, 0.92]',
			'RenderView1.CameraPosition = [%10.5f, %10.5f, %10.5f]'%(xm+2.5*xr,ym+1.5*yr,zm+1.5*zr),
			'RenderView1.CameraFocalPoint = [%10.5f, %10.5f, %10.5f]'%(xm,ym,zm),
			'',
			'RenderView1.ViewTime = %4i'%viewTime,
			'',
			'my_representation0 = GetDisplayProperties(temp_vtk)',
			'my_representation0.Representation = \'Surface With Edges\'',
			'',
			'a1_PVLookupTable = GetLookupTableForArray( '+initial_display+', 1, RGBPoints=[%4.2f, 0.23, 0.299, 0.754, %4.2f, 0.706, 0.016, 0.15], VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )'%tuple(lim),
			'',
			'a1_PiecewiseFunction = CreatePiecewiseFunction( Points=[%4.2f, 0.0, 0.5, 0.0, %4.2f, 1.0, 0.5, 0.0] )'%tuple(lim),
			'',
			'my_representation0.ScalarOpacityFunction = a1_PiecewiseFunction',
			'my_representation0.ColorArrayName = (\'POINT_DATA\', '+initial_display+')',
			'my_representation0.LookupTable = a1_PVLookupTable',
			'',
			'a1_PVLookupTable.ScalarOpacityFunction = a1_PiecewiseFunction',
			'',
			'ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='+initial_display+', LabelFontSize=12, Enabled=1, TitleFontSize=12 )',
			'GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)',
			'',
			'a1_PVLookupTable = GetLookupTableForArray('+initial_display+', 1 )',
			'',
			'ScalarBarWidgetRepresentation1.LookupTable = a1_PVLookupTable',
			'',
			'',
		]
		f.writelines('\n'.join(lns))
		f.close()
		
	def _get_filename(self): return self.path.absolute_to_file+slash+self.path.filename
	filename = property(_get_filename) #: (**)