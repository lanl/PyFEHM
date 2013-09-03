"""Environment file for PyFEHM. Set default attribute values."""

class fdflt(object):
	def __init__(self):
		# material properties - these values will be assigned as defaults if not otherwise set
		self.permeability			=	1.e-15
		self.conductivity			=	2.2
		self.density				=	2500.
		self.porosity				=	0.1	
		self.specific_heat			=	1.e3	
		self.youngs_modulus			=	1.e4 	# MPa
		self.poissons_ratio			=	0.25
		self.pressure_coupling		=	1.
		self.thermal_expansion 		=	3.e-5	# / K

		# output data formats
		self.hist_format 			= 	'tec'
		self.cont_format 			= 	'surf'

		# set this to the fehm executable to be used if no default assigned
		self.fehm_path 				=	'c:\\users\\264485\\fehm\\source\\src\\fehm.exe'
		self.files 					=	['outp','hist','check']
		self.co2_interp_path 		= 	'c:\\users\\264485\\python\\pyfehm\\co2_interp_table.txt'

		# fdata booleans
		self.associate 				= 	True		# associate macro, zone information with nodes
		self.sticky_zones 			= 	True		# print zone definitions immediately before use in input file
		
