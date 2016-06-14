# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 04:19:43 2015

@author: mohit
"""

import numpy as np
from collision import CollisionManager
from initialiser import Initialiser
from movement import MovementManager
from sampler import SamplingManager
import boundary as dm_b
import time as time


class DsmcSolver:
    def __init__(self, cells, gas, domain, surf_group,n_par_in_cell,
                 dt, n_steps, Detector, Collider, Reflector,
                 ignore_frac=0.1, multiphase=False):
        
        self.n_steps = n_steps
        self.cells = cells
        self.open = domain.is_open()
        self.particles, self.cells_in = Initialiser.run(cells, gas, surf_group,
                                                        domain, n_par_in_cell,
                                                        not self.open, multiphase)
        
        self.collision_manager = CollisionManager(cells, self.particles, gas,
                                                  dt, Detector, Collider)
        
        self.sampling_manager = SamplingManager(cells, self.cells_in, 
                                            gas.get_n_species(), self.particles,
                                            n_steps, ignore_frac)
        
        self.movement_manager = MovementManager(self.particles, surf_group, 
                                                Reflector)
        self.dt = dt
        self.temperature = np.zeros(len(cells.get_temperature()))
        self.number_density = np.zeros((gas.get_n_species(), 
                                        len(cells.get_temperature())))
        if (self.open):
            self.boundary_manager = dm_b.BoundaryManager(cells, domain, gas, 
                                        self.particles, n_par_in_cell)
            self.boundary_manager.run([], self.dt)
    
    
    def run(self):
        move_time = 0.0
        dist_time = 0.0
        col_time = 0.0
        samp_time = 0.0
        bound_time = 0.0
        for i in range(self.n_steps):
            print "time_sample = ", i, " ",
            start = time.time()
            self.movement_manager.move_all(self.dt)
            move_time += time.time() - start
            
            start = time.time()
            self.cells.distribute_all_particles(self.particles)
            dist_time += time.time() - start
            
            start = time.time()
            self.collision_manager.collide()
            col_time += time.time() - start
            
            start = time.time()
            self.sampling_manager.sample()
            samp_time += time.time() - start
            
            if (self.open):
                start = time.time()
                particles_out = self.cells.get_particles_out()
#                print "length of particles out = ", len(particles_out)
                self.boundary_manager.run(particles_out, self.dt)
                bound_time += time.time() - start
        
        print "movement time = ", move_time
        print "distribution time = ", dist_time
        print "collisio time = ", col_time
        print "sampling time = ", samp_time
        print "bound_time = ", bound_time
    
    
    def _find_moving_particles(self, reflected_particles):
        index_list = range(self.particles.get_particles_count())
        for index in reflected_particles:
            np.delete(index_list, index)
        return index_list
    
    
    # this function should return temperature in 2d numpy array
    def get_2d_temperature(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, 
                                   self.sampling_manager.get_temperature())

    
    # this function should return mach in 2d numpy array
    def get_2d_speed(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, 
                                   self.sampling_manager.get_speed())

    
    # this function should return mach in 2d numpy array
    def get_2d_u(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, 
                                   self.sampling_manager.get_u())

    
    # this function should return mach in 2d numpy array
    def get_2d_v(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, 
                                   self.sampling_manager.get_v())

    
    # this function should return mach in 2d numpy array
    def get_2d_w(self, cell_x, cell_y):
        return self._convert_to_2d(cell_x, cell_y, 
                                   self.sampling_manager.get_w())
    
    
    # this function should return number density in 2d numpy array
    def get_2d_num_den(self, cell_x, cell_y, tag=0):
        return self._convert_to_2d(cell_x, cell_y,
                                   self.sampling_manager.get_number_density(tag))
    
    
    # this function should convert 1d numpy array to 2d numpy array given
    # x dimension and y dimension.
    def _convert_to_2d(self, cell_x, cell_y, array):
        cell_x, cell_y = int(cell_x), int(cell_y)
        new_array = np.zeros((cell_y, cell_x), dtype=float)
        
        for i in range(cell_x):
            for j in range(cell_y):
                new_array[i][j] = array[j * cell_x + i]
        
        return new_array
