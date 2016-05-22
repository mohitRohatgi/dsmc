# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:13:46 2015

@author: mohit
"""

import numpy as np
import math as math
#from numba import jit

class CollisionManager:
    def __init__(self, cells, particles, gas, dt, Detector, Collider):
        self.detector = Detector(cells, particles, gas, dt)
        self.ref_max_area = np.ones(len(cells.particles_inside))
        self.ref_max_area *= (particles.get_dia(0) * particles.get_dia(0)
                             * np.pi / 2.0) * 300.0
        self.rem_cols = np.zeros(len(cells.particles_inside), dtype = float)
        self.cells = cells
        self.particles = particles
        self.reduced_mass = gas.reduced_mass
        self.n_species = gas.get_n_species()
        self.collider = Collider(particles)
        self.dt = dt
        self.tag = particles.tag
        self.ref_temp = np.zeros(len(self.particles.species))
        self.visc_index = np.zeros_like(self.ref_temp)
        self.dia = np.zeros_like(self.ref_temp)
        self.n_particles = np.zeros(len(cells.particles_inside), dtype=int)
        
        for i, mol in enumerate(self.particles.species):
            self.ref_temp[i] = mol.ref_temp
            self.visc_index[i] = mol.visc_index
            self.dia[i] = mol.dia
    
    
    def collide(self):
        col_pair, rel_speed = self.detector.run()
#        col_pair, rel_speed = detect_collision(self.cells.particles_inside, 
#                    self.cells.cell_volume, self.dt, self.dia, 
#                    self.particles.n_eff, self.reduced_mass, self.n_species,
#                    self.particles.u, self.particles.v, self.particles.w,
#                    self.ref_temp, self.visc_index, self.ref_max_area, 
#                    self.rem_cols, self.tag, self.n_particles)
        
        self.collider.run(col_pair, rel_speed)

#@jit
def detect_collision(particles_inside, cell_volume, dt, dia, n_eff, reduced_mass,
                     n_species, u, v, w, ref_temp, visc_index, ref_max_area,
                     rem_cols, tag, n_particles):
    
    
    # finding number of collisions ...
    probability = n_eff * ref_max_area  * dt / cell_volume
    for cell_index in range(len(n_particles)):
        n_particles[cell_index] = len(particles_inside[cell_index])
    collisions = 0.5 * n_particles * (n_particles - 1) * probability
    collisions += rem_cols
    int_cols = collisions.astype(int)
    rem_cols = collisions - int_cols
    
    # finding random numbers ...
    threshold = np.random.random(sum(int_cols))
    next_index = 0
    detected = []
    rel_speed = []
    
    # detecting the collisions for each cell...
    # iterating for every cell...
    for cell_index, n_particle in enumerate(n_particles):        
        # restricting collisions to not exceed the max possible no. of collision,
        # max possible no. of collisions = no. of particles in cell / 2.
        max_poss_cols = n_particle / 2
        if int_cols[cell_index] > max_poss_cols:
            rem_cols[cell_index] += int_cols[cell_index] - max_poss_cols
            int_cols[cell_index] = max_poss_cols
        
        if n_particle > 1:
#             generating random pairs for checking collision ...
            col_pair = np.random.permutation(n_particles[cell_index])

            for col_num in range(int_cols[cell_index]):
                index1 = col_pair[2 * col_num]
                index2 = col_pair[2 * col_num + 1]
                index1 = particles_inside[cell_index][index1]
                index2 = particles_inside[cell_index][index2]
                
                # finding relative speed ...
                u_rel = u[index1] - u[index2]
                v_rel = v[index1] - v[index2]
                w_rel = w[index1] - w[index2]
                relative_speed = np.sqrt(u_rel * u_rel  + v_rel * v_rel + 
                                         w_rel * w_rel)
                
                # finding collision area ...
                tag1, tag2 = tag[index1], tag[index2]
                d_ref = (dia[tag1] +  dia[tag2]) * 0.5
                T_ref = (ref_temp[tag1] + ref_temp[tag2]) * 0.5
                omega_ref = (visc_index[tag1] + visc_index[tag2]) * 0.5

                k = 1.3806488e-23
                if relative_speed < 1.0e-6:
                    col_area = 0.0
                
                else:
                    constt = (2.0 * k * T_ref / reduced_mass / (relative_speed * 
                              relative_speed)) ** (omega_ref - 0.5)
                    constt /= math.gamma(2.5 - omega_ref)
                    col_area = relative_speed * (0.5 * np.pi * d_ref * d_ref) * constt
                
                # updating reference maximum area ...
                if (ref_max_area[cell_index] < col_area):
                    ref_max_area[cell_index] = col_area
                
                # finding the probability of collision and performing collision
                # accordingly...
                probability = col_area / ref_max_area[cell_index]
                if (probability > threshold[next_index]):
                    detected.append((index1, index2))
                    rel_speed.append(relative_speed)
                
                # increasing the index for use in the next iteration
                next_index += 1
                
        else:
            # adding up collision to remaining collision if there is no particle
            # to collide in the cell
            rem_cols[cell_index] += int_cols[cell_index]
    return (detected, rel_speed)



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
                             * np.pi / 2.0) * 300.0
        
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
        # n_particle = len(self.cells.get_particles_inside(cell_index))
        for cell_index in range(len(self.n_particle)):
            self.n_particle[cell_index] = len(self.cells.get_particles_inside(cell_index))
        collisions = 0.5 * self.n_particle * (self.n_particle - 1) * probability
        collisions += self.remaining_collisions
        self.int_collisions = collisions.astype(int)
        self.remaining_collisions = collisions - self.int_collisions
        
    
    
    # this function returns list of particles detected for collision in a cell
    # identified by cell_index as the first element and their relative speed
    # as the second element.
    def _detect_cell(self, detected, rel_speed, cell_index):
        for col_num in range (self.int_collisions[cell_index]):
            self._detect_pair(detected, rel_speed, cell_index)
    
    
    # this function appends particle pair and their relative speed to their 
    # respective list
    # most time consuming ...
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
            index2 = int(np.random.randint(length))
        index = (self.uncol_particles[cell_index][index1], 
                 self.uncol_particles[cell_index][index2])
        return (index)
    
    
    def _find_relative_speed(self, index1, index2):
        u_rel = self.particles.get_velx(index1) - self.particles.get_velx(index2)
        v_rel = self.particles.get_vely(index1) - self.particles.get_vely(index2)
        w_rel = self.particles.get_velz(index1) - self.particles.get_velz(index2)
        relative_speed =  np.sqrt(u_rel ** 2.0  + v_rel ** 2.0 + w_rel ** 2.0)
        return (relative_speed)
    
    
    # pair represents the pair of particles to check for collision.
    # it constains the indices of particles in particle array while cell_index
    # contains the cell_index in cell array.
    def _check(self, pair, cell_index, relative_speed):
        col_area = self._find_col_area(pair, cell_index, relative_speed)
        probability = col_area / self.ref_max_area[cell_index]
        return probability > np.random.rand()
    
    
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
        constt = (2.0 * k * T_ref / (self.reduced_mass * relative_speed * 
                  relative_speed)) ** (omega_ref - 0.5)
        constt /= math.gamma(2.5 - omega_ref)
        col_area = relative_speed * (0.5 * np.pi * d_ref * d_ref) * constt
        
        if (self.ref_max_area[cell_index] < col_area):
            self.ref_max_area[cell_index] = col_area
        
        return col_area
