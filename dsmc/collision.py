# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:13:46 2015

@author: mohit
"""

import numpy as np
import math as math
import time as time


class CollisionManager:
    def __init__(self, cells, particles, gas, dt, Detector, Collider):
        self.detector = Detector(cells, particles, gas, dt)
        self.collider = Collider(particles)
    
    
    def set_dt(self, dt):
        self.dt = dt
    
    
    def collide(self):
        collision_pair, rel_speed = self.detector.run()
        self.collider.run(collision_pair, rel_speed)



# this class is a binary collision detector. It doesn't takes into account
# the trace particle collision problem.
class CollisionDetector:    
    # setup needs to be called after each time step.
    def __init__(self, cells, particles, gas, dt):
        self.particles = particles
        self.cells = cells
        self.n_species = gas.get_n_species()
        self.reduced_mass = gas.get_reduced_mass()
        self.dt = dt
        length = len(cells.get_temperature())
        self.ref_max_area = np.ones(length)
        self.ref_max_area *= (self.particles.get_dia(0) * self.particles.get_dia(0)
                             * np.pi / 2.0)
        
        u_max = max(self.particles.get_velx())
        v_max = max(self.particles.get_vely())
        w_max = max(self.particles.get_velz())
        self.ref_max_area *= (u_max * u_max + v_max * v_max + w_max * w_max)
        
        self.n_particle = np.zeros(length, dtype=int)
        self.int_collisions = np.zeros(length, dtype = int)
        self.remaining_collisions = np.zeros(length, dtype = float)
        self.dsmc_collisions = np.zeros(length, dtype = int)
        self.uncol_particles = []
    
    
    # Returns a tuple of 2 elements. The first element is list of particle pairs
    # to collide and the second element would be their respective relative speed.
    # Both pairs and relative speeds are contained in a list .
    # data structure can be represented as (for 3 cells each having 1 collision):
    # ([(x, y), (p, q), (m, n)], [a, b, c])
    def run(self):
        collision_pair = []
        rel_speed = []
        self.uncol_particles = [list(self.cells.get_particles_inside(i))
                            for i in range(len(self.cells.get_temperature()))]
        
        self._find_collisions()
        for cell_index in range(len(self.cells.get_temperature())):
            total_collided = len(rel_speed)
            self._detect_cell(collision_pair, rel_speed, cell_index)
            self.dsmc_collisions[cell_index] = len(rel_speed) - total_collided
        
        return (collision_pair, rel_speed)
    
    
    # assuming all the particle represents same number of molecules, n_eff.
    # assuming the cells are rectangular.
    def _find_collisions(self):
        probability = self.particles.get_n_eff() * self.ref_max_area * self.dt
        probability /= self.cells.get_cell_volume()
        #n_particle = len(self.cells.get_particles_inside(cell_index))
        self._find_n_particle()
        collisions = 0.5 * self.n_particle * (self.n_particle - 1) * probability
        self.int_collisions = (collisions + self.remaining_collisions).astype(int)
        self.remaining_collisions += collisions - self.int_collisions
    
    
    # helper function for _find_collision to find the array of no. of particles
    # each cell has.
    # introducing this function as it is not easy to vectorize...
    def _find_n_particle(self):
        for cell_index in range(len(self.n_particle)):
            self.n_particle[cell_index] = len(self.cells.get_particles_inside(cell_index))
    
    
    # this function returns list of particles detected for collision in a cell
    # identified by cell_index as the first element and their relative speed
    # as the second element.
    def _detect_cell(self, detected, rel_speed, cell_index):
        for i in range (self.int_collisions[cell_index]):
            self._detect_pair(detected, rel_speed, cell_index)
    
    
    # this function appends particle pair and their relative speed to their 
    # respective list
    # most time consuming
    def _detect_pair(self, detected, rel_speed, cell_index):
        length =  len(self.uncol_particles[cell_index])
        if length > 1:
            pair = self._select_pair(length, cell_index)
            index1 = int(pair[0])
            index2 = int(pair[1])
            relative_speed = self._find_relative_speed(index1, index2)
            if (self._check((index1, index2), cell_index, relative_speed)):
                detected.append((index1, index2))
                rel_speed.append(relative_speed)
                self.uncol_particles[cell_index].remove(index1)
    
    
    # index are indices of cells[cell_index.particles_inside and 
    # represents the particle inside the cell that would checked for collision.
    # this function returns the index of particle in particle array.
    def _select_pair(self, length, cell_index):
        index1 = int(np.random.randint(length))
        index2 = int(np.random.randint(length))
        while (index1 == index2):
            index2 = int(np.random.rand() * length)
        index = (self.uncol_particles[cell_index][index1], 
                 self.uncol_particles[cell_index][index2])
        return (index)
    
    
    # pair represents the pair of particles to check for collision.
    # it constains the indices of particles in particle array while cell_index
    # contains the cell_index in cell array.
    def _check(self, pair, cell_index, relative_speed):
        col_area = self._find_col_area(pair, cell_index, relative_speed)
        probability = col_area / self.ref_max_area[cell_index]
        threshold = np.random.rand()
        if (probability > threshold):
            return True
        else:
            return False
    
    
    def _find_relative_speed(self, index1, index2):
        u_rel = self.particles.get_velx(index1) - self.particles.get_velx(index2)
        v_rel = self.particles.get_vely(index1) - self.particles.get_vely(index2)
        w_rel = self.particles.get_velz(index1) - self.particles.get_velz(index2)
        relative_speed =  np.sqrt(u_rel ** 2.0  + v_rel ** 2.0 + w_rel ** 2.0)
        return (relative_speed)
    
    
    # this function finds the collision area and updates the max collision area.
    # most time consuming...
    def _find_col_area(self, pair, cell_index, relative_speed):
        if np.abs(relative_speed) < 1.0e-8:
            return 0.0
        
        d_ref = self.particles.get_dia(pair[0]) +  self.particles.get_dia(pair[1])
        d_ref *= 0.5
        
        T_ref = self.particles.get_ref_temp(pair[0]) + self.particles.get_ref_temp(pair[1])
        T_ref *= 0.5
        
        omega_ref = self.particles.get_visc_index(pair[0]) + self.particles.get_visc_index(pair[1])
        omega_ref *= 0.5
        
        k = 1.3806488e-23
        constt = (2.0 * k * T_ref / self.reduced_mass / (relative_speed * 
                  relative_speed)) ** (omega_ref - 0.5)
        constt /= math.gamma(2.5 - omega_ref)
        col_area = relative_speed * (0.5 * np.pi * d_ref * d_ref) * constt
        
        
        if (self.ref_max_area[cell_index] < col_area):
            self.ref_max_area[cell_index] = col_area
        
        return col_area
