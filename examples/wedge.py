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
from dsmc.dsmc.reflection_models import Specular, Diffuse
from dsmc.dsmc.collider_models import VhsCollider
from dsmc.dsmc.collision import CollisionDetector
from dsmc.dsmc.geometry import SurfaceGroup, Domain

"""
'_col_' denotes collision while '_f_' denotes free
"""


def main():
    
    wedge = wedge = [((0.0, 0.6268), (0.1, 1.0))]
    ref_point = (0.0, 1.0)
    surf_temp = 1000.0
    # creating surface
    surf_group = SurfaceGroup()
    surf_group.add_new_group(wedge, ref_point, surf_temp)
    
    inlet = [((0.0, 0.0), (1.0, 0.0))]
    outlet = [((1.0, 1.0), (0.0, 1.0))]
    zero_grad  = [((0.0, 0.0), (0.0, 0.6268)),
		  ((1.0, 0.0), (1.0, 1.0))]
    center = (0.5, 0.5)
    length = 1.0
    width = 1.0
    volume = 1.0 - 0.01866
    domain = Domain(volume, inlet, zero_grad, outlet)
    
    
#    ensemble_sample = 10
    time_av_sample = 1
    dof = 3.0
    mass = 66.3e-27
    viscosity_coeff = 2.117
    viscosity_index = 0.81
    mole_fraction = [0.5, 0.5]
    dia = 4.17e-10
    mach = [0.0, 10.0, 0.0]
    temperature = 200.0
    ref_temperature = 273.0
    number_density = 1.699e20
    gamma = 5.0 / 3.0
    n_particles_in_cell = 30
    argon = Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 0,
                ref_temperature, gamma, volume, number_density)
    argon1 = Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                ref_temperature, gamma, volume, number_density)
    gas = Gas([argon, argon1], mole_fraction, mach, temperature)
    gas.setup()
#    dl = min(gas.mean_f_path)
#    dt = min(gas.mean_col_time)
    dt = 1.0e-7
#    print dt
#    cell_x, cell_y = np.ceil(length / dl), np.ceil(width / dl)
    cell_x, cell_y = 250, 250
    cells = RectCells(cell_x, cell_y, length, width, center, 2)
    
    start_time = time()
    solver = DsmcSolver(cells, gas, domain, surf_group, n_particles_in_cell,
                         dt, time_av_sample, CollisionDetector,
                        VhsCollider, Specular, ignore_frac=0.25)
    solver.run()
    end_time = time()
    
    simulation_time = end_time - start_time
    print "simulation time in mins (upper bound) = ", round(simulation_time / 60)
    print "simulation time in sec = ", simulation_time
    number_density_0 = solver.get_2d_num_den(cell_x, cell_y, 0)
    number_density_1 = solver.get_2d_num_den(cell_x, cell_y, 1)
    temperature = solver.get_2d_temperature(cell_x, cell_y)
    speed = solver.get_2d_speed(cell_x, cell_y)
    u = solver.get_2d_u(cell_x, cell_y)
    v = solver.get_2d_v(cell_x, cell_y)
    w = solver.get_2d_w(cell_x, cell_y)
    pressure = (number_density_0 + number_density_1) * 1.3806488e-23 * temperature
    
    f = open('wedge_super_temperature.txt', 'w')
    np.savetxt(f, temperature)
    f.close()
    
    f = open('wedge_super_number_density_0.txt', 'w')
    np.savetxt(f, number_density_0)
    f.close()
    
    f = open('wedge_super_number_density_1.txt', 'w')
    np.savetxt(f, number_density_1)
    f.close()
    
    f = open('wedge_super_speed.txt', 'w')
    np.savetxt(f, speed)
    f.close()
    
    f = open('wedge_super_u.txt', 'w')
    np.savetxt(f, u)
    f.close()
    
    f = open('wedge_super_v.txt', 'w')
    np.savetxt(f, v)
    f.close()
    
    f = open('wedge_super_w.txt', 'w')
    np.savetxt(f, w)
    f.close()
    
    f = open('wedge_super_pressure.txt', 'w')
    np.savetxt(f, pressure)
    f.close()


if __name__ == '__main__':
    main()
