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

import numpy as np
from time import time
from dsmc.dsmc.cells import RectCells
from dsmc.dsmc.solver import DsmcSolver
from dsmc.dsmc.particles import Molecules, Gas
from dsmc.dsmc.reflection_models import Diffuse
from dsmc.dsmc.collider_models import VhsCollider
from dsmc.dsmc.collision import CollisionDetector
from dsmc.dsmc.geometry import SurfaceGroup, Domain


"""
'_col_' denotes collision while '_f_' denotes free
"""


def main():

    surface1 = [((0.0, 1.0), (0.0, 0.0)), ((0.0, 0.0), (1.0, 0.0)), ((1.0, 0.0), (1.0, 1.0))]
    surface1_temp = 300.0
    
    surface2 = [((1.0, 1.0), (0.0, 1.0))]
    surface2_temp = 500.0
    
    center = (0.5, 0.5)
    length = 1.0
    width = 1.0
    volume = 1.0
    
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
    
    argon = Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 0,
                ref_temperature, gamma, volume, number_density)
    argon1 = Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                ref_temperature, gamma, volume, number_density)
    gas = Gas([argon, argon1], mole_fraction, mach, temperature)
    gas.setup()
    
    domain = Domain(volume)
    surf_group = SurfaceGroup()
    surf_group.add_new_group(surface1, center, surface1_temp)
    surf_group.add_new_group(surface2, center, surface2_temp)
    dt = 1.0e-5
    cell_x, cell_y = 5, 5
    cells = RectCells(cell_x, cell_y, length, width, center, 2)
    
    start_time = time()
    solver = DsmcSolver(cells, gas, domain, surf_group, n_particles_in_cell,
                        ref_point, dt, sample_size, CollisionDetector,
                        VhsCollider, Diffuse)
    solver.run()
    end_time = time()
    
    simulation_time = end_time - start_time
    print "simulation time in mins (upper bound) = ", int(simulation_time / 60) + 1
    print "simulation time in sec = ", simulation_time
    number_density = solver.get_2d_num_den(10, 10)
    print "number_density = ", number_density
    temperature = solver.get_2d_temperature(10, 10)
    print "temperature = ", temperature
    
    f = open('box_temp_grad_temperature.txt', 'w')
    np.savetxt(f, temperature)
    f.close()
    
    f = open('box_temp_grad_number_density.txt', 'w')
    np.savetxt(f, number_density)
    f.close()


if __name__ == '__main__':
    main()
