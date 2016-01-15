# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 04:19:43 2015

@author: mohit
"""

import numpy as np
import dsmc_cells as dm_c
import dsmc_collision as dm_col
import dsmc_initialiser as dm_i
import dsmc_reflector as dm_r
import dsmc_sampler as dm_s
import dsmc_boundary as dm_b


class DsmcSolver:
    def __init__(self, cells, gas, domain, n_particles_in_cell,
                 ref_point, dt, n_steps):
        
        self.n_steps = n_steps
        self.cells = cells
        self.particles, self.cells_in = dm_i.Initialiser.run(cells, gas, domain,
                                                n_particles_in_cell, ref_point)
        
        self.distributor = dm_c.Distributor(cells, self.particles)
        
        self.collision_manager = dm_col.CollisionManager(cells, self.particles,
                                                    gas.get_reduced_mass(), dt)
        
        self.sampling_manager = dm_s.Sampling_manager(cells, self.cells_in, 
                                gas.get_n_species(), self.particles, n_steps)
        
        self.flag = False
        self.movement_manager = dm_r.MovementManager(self.particles,
                                                     domain.get_surface())
        self.dt = dt
        self.temperature = np.zeros(len(cells.get_temperature()))
        self.number_density = np.zeros((gas.get_n_species(), 
                                        len(cells.get_temperature())))
        if (domain.is_open()):
            self.flag = True
            self.boundary_manager = dm_b.BoundaryManager(cells, domain, gas, 
                                        self.particles, n_particles_in_cell)
            self.boundary_manager.run([], self.dt)
    
    
    def run(self, detector_key, collider_key, reflection_key):
        
        for i in range(self.n_steps):
            self.movement_manager.move_all(reflection_key, self.dt)
            self.distributor.distribute_all_particles()
            self.collision_manager.run(collider_key, detector_key)
            self.sampling_manager.run()
            if (self.flag):
                particles_out = self.distributor.get_particles_out()
                self.boundary_manager.run(particles_out, self.dt)
        self.temperature = self.sampling_manager.get_temperature()
        self.number_density = self.sampling_manager.get_number_density()
        constt = 0.0
        for temp in self.temperature:
            constt += temp
        
        print "system temperature = ", temp
    
    
    def _find_moving_particles(self, reflected_particles):
        index_list = range(self.particles.get_particles_count())
        for index in reflected_particles:
            np.delete(index_list, index)
        return index_list
    
    
    # this function should return number density in 1d numpy array
    def get_1d_num_den(self):
        return self.number_density
    
    
    # this function should return temperature in 1d numpy array
    def get_1d_temperature(self):
        return self.temperature
    
    
    # this function should return temperature in 2d numpy array
    def get_2d_temperature(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, self.temperature)
    
    
    # this function should return number density in 2d numpy array
    def get_2d_num_den(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, self.number_density)
    
    
    # this function should convert 1d numpy array to 2d numpy array given
    # x dimension and y dimension.
    def _convert_to_2d(self, cell_x, cell_y, array):
        cell_x, cell_y = int(cell_x), int(cell_y)
        new_array = np.zeros((cell_x, cell_y), dtype=float)
        index = cell_x * cell_y
        
        for i in range(cell_x):
            index -= cell_x
            for j in range(cell_y):
                new_array[i][j] = array[index + j]
        
        return new_array