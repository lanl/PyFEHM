# 7.1 PyFEHM tutorial 1
# Cube with fixed pressures/temperatures
# 7.1.1 Getting started

#import os,sys
#sys.path.append('c:\\python\\pyfehm')

from fdata import*
from fpost import*

root = 'tut1'
dat = fdata()

# 7.1.2 Grid generation
x = np.linspace(0,10,11)
dat.grid.make(root+'_GRID.inp',x=x,y=x,z=x)
dat.grid.plot(root+'_GRID.png',color='r',angle=[45,45])

# 7.1.3 Zone creation
zn = fzone(index=1,name='lower')
zn.rect([-0.,-0.,-0.1],[10.1,10.1,3.1])
dat.add(zn)

zn = fzone(index=2,name='middle')
zn.rect([-0.,-0.,3.1],[10.1,10.1,6.1])
dat.add(zn)

zn = fzone(index=3,name='upper')
zn.rect([-0.,-0.,6.1],[10.1,10.1,10.1])
dat.add(zn)

# 7.1.4 Adding macros
rm = fmacro('rock',param=(('density',2500),('specific_heat',1000),('porosity',0.1)))
dat.add(rm)

pm = fmacro('perm',zone=1,param=(('kx',1.e-15),('ky',1.e-15),('kz',1.e-16)))
dat.add(pm)

kmid = 1.e-20
pm = fmacro('perm',zone=2,param=(('kx',kmid),('ky',kmid),('kz',kmid)))
dat.add(pm)

pm = fmacro('perm',zone='upper',param=(('kx',1.e-14),('ky',1.e-14),('kz',1.e-14)))
dat.add(pm)

# 7.1.5 Initial and boundary conditions
pres = fmacro('pres',param=(('pressure',5.),('temperature',60.),('saturation',1)))
dat.add(pres)

hflx = fmacro('hflx',zone='ZMAX',param=(('heat_flow',30),('multiplier',1.e10)))
dat.add(hflx)

dat.zone['ZMIN'].fix_temperature(80)

flow = fmacro('flow',zone='YMIN',param=(('rate',6.),('energy',-60.),('impedance',1.e6)))
dat.add(flow)

dat.zone['YMAX'].fix_pressure(P=4.,T=60.)

# 7.1.6 Running the simulation
dat.cont.variables.append(['xyz','temperature','pressure'])

dat.tf=10.
#dat.time['max_time_TIMS']=10.

dat.files.root = root
dat.run(root+'_INPUT.dat',exe='fehm.exe',files=['outp'])

# 7.1.7 Visualisation
c = fcontour('*.csv',latest=True)
c.slice_plot_fill(save='Tslice.png',cbar=True,levels=11,slice=['x',5],variable='T',method='linear',title='temperature / degC',
xlabel='y / m', ylabel = 'z / m')
c.slice_plot_fill(save='Pslice.png',cbar=True,levels=np.linspace(4,6,9),slice=['x',5],variable='P',method='linear',title='pressure / MPa',
xlabel='y / m', ylabel = 'z / m')
