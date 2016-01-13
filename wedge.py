# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 22:49:35 2015

@author: mohit
"""

import dsmc_particles as dm_p
import dsmc_cells as dm_c
import dsmc_solver as dm_sol
from time import time

"""
'_col_' denotes collision while '_f_' denotes free
"""


def main():
    inlet = [((0.0, 0.0), (0.0, 1.0))]
    outlet = [((1.0, 1.0), (1.0, 0.8))]
    zero_grad  = [((0.2, 0.0), (0.0, 0.0)), ((0.0, 1.0), (1.0, 1.0))]
    surface = ((0.2, 0.0), (1.0, 0.8))
#    surface_temp = 300.0
    
    center = (0.5,0.5)
    length = 1.0
    width = 1.0
    volume = 0.68
    domain = dm_p.Domain(inlet, zero_grad, outlet)
    domain.set_surface(surface)
#    domain.set_surface(surface, surface_temp)
    domain.set_volume(volume)
    
#    ensemble_sample = 10
    time_av_sample = 10
    dof = 3.0
    mass = 66.3e-27
    viscosity_coeff = 2.117
    viscosity_index = 0.81
    mole_fraction = [0.5, 0.5]
    dia = 4.17e-10
    mach = [2.5, 0.0, 0.0]
    temperature = 500.0
    ref_temperature = 273.0
    number_density = 1.699e19
    gamma = 5.0 / 3.0
    n_particles_in_cell = 10
    ref_point = (0.1, 0.5)
    argon = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 0,
                ref_temperature, gamma, volume, number_density)
    argon1 = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                ref_temperature, gamma, volume, number_density)
    gas = dm_p.Gas([argon, argon1], mole_fraction, mach, temperature)
    gas.setup()
#    dl = min(gas.mean_f_path)
#    dt = min(gas.mean_col_time)
    dt = 1.0e-5
#    print dt
#    cell_x, cell_y = np.ceil(length / dl), np.ceil(width / dl)
    cell_x, cell_y = 10, 10
    cells = dm_c.RectCells(cell_x, cell_y, length, width, center, 2)
    
    start_time = time()
    solver = dm_sol.DsmcSolver(cells, gas, domain, n_particles_in_cell,
                               ref_point, dt, time_av_sample)
    solver.run(1, 1, 1)
    end_time = time()
    simulation_time = end_time - start_time
    print "simulation time in mins (upper bound) = ", int(simulation_time / 60) + 1
    print "simulation time in sec = ", simulation_time
    number_density = solver.get_2d_num_den(cell_x, cell_y)
    print "number_density = ", number_density
    temperature = solver.get_2d_temperature(cell_x, cell_y)
    print "temperature = ", temperature
    
    num_den_msg = 'this file contains number density( x1e18 ) of each cell'
    temp_msg = 'this file contains temperature of each cell'
    
    dump_output('wedge_super_temperature.txt', temperature,  temp_msg)
    dump_output('wedge_super_number_density.txt', number_density / 1e18, num_den_msg)



def dump_output(filename, data, data2=None, msg=None):
    f = open(filename, 'w')
    if msg != None:
        print >> f, msg
    if data2 != None:
        print >> f, data2
    print >> f, data
    f.close



if __name__ == '__main__':
    main()
