# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 18:46:14 2016

@author: mohit
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 22:49:35 2015

@author: mohit
"""

import dsmc_particles as dm_p
import dsmc_cells as dm_c
import dsmc_solver as dm_sol
import dsmc_geometry as dm_g
from time import time

"""
'_col_' denotes collision while '_f_' denotes free
"""


def main():

    surface = [((0.0, 1.0), (0.0, 0.0)), ((0.0, 0.0), (1.0, 0.0)),
            ((1.0, 0.0), (1.0, 1.0)), ((1.0, 1.0), (0.0, 1.0))]
    
    center = (0.5,0.5)
    length = 1.0
    width = 1.0
    volume = 1.0
    
    domain = dm_g.Domain(volume)
    surf_group = dm_g.SurfaceGroup()
    surf_group.add_new_group(surface, center)
    
#    ensemble_sample = 10
    sample_size = 10
    dof = 3.0
    mass = 66.3e-27
    viscosity_coeff = 2.117
    viscosity_index = 0.81
    mole_fraction = [0.5, 0.5]
    dia = 4.17e-10
    mach = [0.0, 0.0, 0.0]
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
    solver = dm_sol.DsmcSolver(cells, gas, domain, surf_group,
                               n_particles_in_cell, ref_point, dt, sample_size)
    solver.run(1, 1, 1)
    end_time = time()
    simulation_time = end_time - start_time
    print "simulation time in mins (upper bound) = ", int(simulation_time / 60) + 1
    print "simulation time in sec = ", simulation_time
    number_density = solver.get_2d_num_den(10, 10)
    print "number_density = ", number_density
    temperature = solver.get_2d_temperature(10, 10)
    print "temperature = ", temperature
    
    f = open('box_temperature.txt', 'w')
    np.savetxt(f, temperature)
    f.close()
    
    f = open('box_number_density.txt', 'w')
    np.savetxt(f, number_density)
    f.close()


if __name__ == '__main__':
    main()