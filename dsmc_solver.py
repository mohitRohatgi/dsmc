# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 04:19:43 2015

@author: mohit
"""

import numpy as np
import dsmc_cells as dm_c
import dsmc_collision as dm_col
import dsmc_initialiser as dm_i
import dsmc_particles as dm_p
import dsmc_reflector as dm_r
import dsmc_sampler as dm_s
import dsmc_boundary as dm_b


class DsmcSolver:
    def __init__(self, cells, gas, domain, n_particles_in_cell, ref_point, dt,
                 n_steps):
        
        self.n_steps = n_steps
        self.particles, self.cells_in = dm_i.Initialiser.run(cells, gas, domain,
                                                n_particles_in_cell, ref_point)
        
        self.distributor = dm_c.Distributor(cells, self.particles)
        
        self.collision_manager = dm_col.CollisionManager(cells, self.particles,
                                                          gas, dt)
        
        self.sampling_manager = dm_s.Sampling_manager(cells, self.cells_in, gas,
                                                      self.particles, n_steps)
        self.flag = False
        if (len(domain.inlet) > 1):
            self.flag = True
            self.boundary = dm_b.Boundary(cells, gas, domain, n_particles_in_cell,
                                      self.particles, ref_point)
        self.movement_manager =dm_r.MovementManager(self.particles, domain.surface,
                                                    domain.surface_temperature)
        self.dt = dt
        self.temperature = np.zeros(len(cells.temperature))
        self.number_density = np.zeros_like(self.temperature)
    
    
    def run(self, detector_key, collider_key, reflection_key):
        for i in range(self.n_steps):
            self.movement_manager.move_all(reflection_key, self.dt)
            if (self.flag):
                self.boundary.run(self.dt)
            self.distributor.distribute_all_particles()
            self.collision_manager.run(collider_key, detector_key)
            self.sampling_manager.run()
        print self.sampling_manager.get_temperature()
#        print self.sampling_manager.cells.temperature
    
    
    def _find_moving_particles(self, reflected_particles):
        index_list = range(len(self.particles.x))
        for index in reflected_particles:
            np.delete(index_list, index)
        return index_list
    
    
    def extract(self):
        return (self.temperature, self.number_density)
        



def main():
    inlet = ((0.2, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.8))
    outlet = ((0.2, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.8))
    surface = ((0.2, 0.0), (1.0, 0.8))
#    inlet = ()
#    outlet = ()
#    surface = ((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0))
    center = (0.5,0.5)
    length = 1.0
    width = 1.0
    volume = 0.36
    ensemble_sample = 100
    time_av_sample = 100
    time = 2.0
    dof = 3.0
    mass = 66.3e-27
    viscosity_coeff = 2.117
    viscosity_index = 0.81
    mole_fraction = [0.5, 0.5]
    dia = 4.17e-10
    mach = [0.0, 0.0, 0.0]
    temperature = 500.0
    surface_temperature = 300.0
    ref_temperature = 273.0
    number_density = 1.699e19
    gamma = 5.0 / 3.0
    n_particles_in_cell = 8
    ref_point = (0.1, 0.5)
    domain = dm_p.Domain(inlet, outlet, surface, volume, surface_temperature)
    argon = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 0,
                ref_temperature, gamma, volume, number_density)
    argon1 = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                ref_temperature, gamma, volume, number_density)
    gas = dm_p.Gas([argon, argon1], mole_fraction, mach, temperature)
    gas.setup()
    dl = min(gas.mean_f_path)
#    dt = min(gas.mean_col_time)
    dt = 1.0e-5
#    print dt
#    cell_x, cell_y = np.ceil(length / dl), np.ceil(width / dl)
    cell_x, cell_y = 10, 10
    cells = dm_c.RectCells(cell_x, cell_y, length, width, center, gas)
    solver = DsmcSolver(cells, gas, domain, n_particles_in_cell, ref_point, dt,
                        time_av_sample)
    solver.run(1, 1, 1)
#    output = solver.extract()
#    output_manager = dm_out.Output_manager()
#    output_manager.show(output)



if __name__ == '__main__':
    main()