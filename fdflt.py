"""Environment file for PyFEHM. Set default attribute values."""

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

import os

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
		
		# initial conditions
		self.Pi						=	1. 		# pressure
		self.Ti						=	30. 	# temperature 
		self.Si						=	1.		# saturation

		# output data formats
		self.hist_format 			= 	'tec'
		self.cont_format 			= 	'surf'
		self.parental_cont 			= 	True

		# set this to the fehm executable to be used if no default assigned
		self.fehm_path 				=	'c:\\users\\264485\\fehm\\source\\src\\fehm.exe'
		self.files 					=	['outp','hist','check']
		self.co2_interp_path 		= 	'c:\\users\\264485\\python\\pyfehm\\co2_interp_table.txt'
		self.co2_interp_path_2 		= 	'/home/ddempsey/python/pyfehm/co2_interp_table.txt'
		if not os.path.isfile(self.co2_interp_path):
			self.co2_interp_path = self.co2_interp_path_2

		# fdata booleans
		self.associate 				= 	True		# associate macro, zone information with nodes
		self.sticky_zones 			= 	True		# print zone definitions immediately before use in input file
		self.full_connectivity 		=	True	
		self.sleep_time 			= 	1.
		self.keep_unknown 			= 	True 		# set true if PyFEHM should preserve unknown macros in future output files
		
		# default values for mactro ITER (parameters controlling solver)
		self.iter = {
			'linear_converge_NRmult_G1':1.e-5, 			# convergence criteria
			'quadratic_converge_NRmult_G2':1.e-5,
			'stop_criteria_NRmult_G3':1.e-3,
			'machine_tolerance_TMCH':-1.e-5,
			'overrelaxation_factor_OVERF':1.1,
			
			'reduced_dof_IRDOF':0,
			'reordering_param_ISLORD':0,
			'IRDOF_param_IBACK':0,
			'number_SOR_iterations_ICOUPL':0,
			'max_machine_time_RNMAX':3600, 				# number of minutes at which FEHM will cut a simulation
			}
		# default values for macro CTRL (parameters controlling simulation)
		self.ctrl = {
			'max_newton_iterations_MAXIT':10,			# solver parameters
			'newton_cycle_tolerance_EPM':1.e-5,         # solver parameters
			'number_orthogonalizations_NORTH':8,        # solver parameters
			
			'max_solver_iterations_MAXSOLVE':24,
			'acceleration_method_ACCM':'gmre',
			
			'JA':1,'JB':0,'JC':0,
			'order_gauss_elim_NAR':2,
			
			'implicitness_factor_AAW':1,
			'gravity_direction_AGRAV':3, 				# direction of gravity
			'upstream_weighting_UPWGT':1.0,
			
			'max_multiply_iterations_IAMM':7,
			'timestep_multiplier_AIAA':1.5, 			# acceleration, time step multiplier
			'min_timestep_DAYMIN':1.e-5, 				# minimum allowable time step (days)
			'max_timestep_DAYMAX':30.,					# maximum allowable time step (days)
			
			'geometry_ICNL':0, 							# problem geometry (0 = 3-D)
			'stor_file_LDA':0 							# flag to use stor file
			}
		# default values for macro TIME
		self.time = {
			'initial_timestep_DAY':1., 					# initial time step size (days)
			'max_time_TIMS':365., 						# maximum simulation time (days)
			'max_timestep_NSTEP':200, 					# maximum number of time steps 
			'print_interval_IPRTOUT':1,					# for printing information to screen
			'initial_year_YEAR':None, 					# initial simulation time (years)	
			'initial_month_MONTH':None,					# (months)
			'initial_day_INITTIME':None					# (years)
			}
		# default values for macro SOL
		self.sol = {
			'coupling_NTT':1,
			'element_integration_INTG':-1
			}
		# default values for macro TRAC
		self.trac = {
			'init_solute_conc_ANO':0.,
			'implicit_factor_AWC':1.,
			'tolerance_EPC':1.e-7,
			'upstream_weight_UPWGTA':0.5,
			
			'solute_start_DAYCS':1.,
			'solute_end_DAYCF':2.,
			'flow_end_DAYHF':1.,
			'flow_start_DAYHS':2.,
			
			'max_iterations_IACCMX':50,
			'timestep_multiplier_DAYCM':1.2,
			'initial_timestep_DAYCMM':1.,
			'max_timestep_DAYCMX':1000.,
			'print_interval_NPRTTRC':1.
			}
			
		self.adsorption = {
			'type_IADSF':None,
			'alpha1_A1ADSF':None,
			'alpha2_A2ADSF':None,
			'beta_BETADF':None
			}