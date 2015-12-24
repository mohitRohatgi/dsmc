# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:13:46 2015

@author: mohit
"""

import numpy as np
import math

def mean(a, b):
    return ((a + b) / 2.0)



class CollisionManager:
    def __init__(self, cells, particles, reduced_mass, dt):
        self.cells = cells
        self.particles = particles
        self.dt = dt
        self.detector = Detector(cells, particles, reduced_mass, dt)
        self.collider = Collider(self.particles)
    
    
    def set_dt(self, dt):
        self.dt = dt
    
    
    def add_collider_key(self, key, object_ref):
        self.collider.model[key] = object_ref
    
    
    def add_detector_key(self, key, object_ref):
        self.detector.model[key] = object_ref
    
    
    def show_detector_key(self):
        for key in self.detector.model:
            print key
    
    
    def show_collider_key(self):
        for key in self.collider.model:
            print key
    
    
    # collision_pairs and rel_speed are not stored as they dont have much use.
    def run(self, collider_key, detector_key):
        collision_pair, rel_speed = self.detector.model[detector_key].run()
        self.collider.run(collision_pair, rel_speed, collider_key)


class Detector:
    def __init__(self, cells, particles, reduced_mass, dt):
        self.cells = cells
        self.particles = particles
        self.dt = dt
        self.collision_pair = []
        self.rel_speed = np.zeros(cells.n_x * cells.n_y)
        self.model = {'binary detector' : CollisionDetector(cells, particles, 
                                                            reduced_mass, dt), 
                    1 : CollisionDetector(cells, particles, reduced_mass, dt)}
    
    
    def add_key(self, key, object_ref):
        self.model[key] = object_ref
    
    
    def run(self, model_key=1):
        detector = self.model[model_key]
        self.collision_pair, self.rel_speed = detector.run()



# this class is a binary collision detector. It doesn't takes into account
# the trace particle collision problem.
class CollisionDetector:    
    # setup needs to be called after each time step.
    def __init__(self, cells, particles, reduced_mass, dt):
        self.particles = particles
        self.cells = cells
        self.reduced_mass = reduced_mass
        self.dt = dt
        length = cells.n_x * cells.n_y
        self.ref_max_area = np.ones(length)
        self.ref_max_area *= self.particles.dia[0] ** 2.0 * np.pi / 2.0 
        self.ref_max_area *= (max(self.particles.u) ** 2.0 + max(self.particles.v)
                                ** 2.0 + max(self.particles.w) ** 2.0)
        
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
        stack1 = []
        stack2 = []
        self.uncol_particles = [list(self.cells.particles_inside[i])
                                for i in range(self.cells.n_x * self.cells.n_y)]
        self._detect_subset(stack1, stack2, 0, self.cells.n_x * self.cells.n_y - 1)
        return (stack1, stack2)
    
    
    def _detect_subset(self, stack1, stack2, start, end):
        if (start == end):
            self._detect_cell(stack1, stack2, start)
            self.dsmc_collisions[start] += len(stack2)
        else:
            mid = int((start + end) / 2)
            self._detect_subset(stack1, stack2, start, mid)
            self._detect_subset(stack1, stack2, mid + 1, end)
    
    
    # this function returns list of particles detected for collision in a cell
    # identified by cell_index as the first element and their relative speed
    # as the second element.
    def _detect_cell(self, detected, rel_speed, cell_index):
        self._find_collisions(cell_index)
        for i in range (self.int_collisions[cell_index]):
            self._detect_pair(detected, rel_speed, cell_index)
    
    
    # assuming all the particle represents same number of molecules, n_eff.
    # assuming the cells are rectangular.
    def _find_collisions(self, cell_index):
        probability = self.particles.n_eff * self.ref_max_area[cell_index]
        probability /= self.cells.volume[cell_index] 
        probability *= self.dt
        n_particle = len(self.cells.particles_inside[cell_index])
        collisions = 0.5 * n_particle * (n_particle - 1) * probability
        self.int_collisions[cell_index] = int(collisions + 
                                        self.remaining_collisions[cell_index])
        self.remaining_collisions[cell_index] += (collisions - 
                                                self.int_collisions[cell_index])
    
    
    # this function appends particle pair and their relative speed to their 
    # respective list
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
#        print col_area
        probability = col_area / self.ref_max_area[cell_index]
        threshold = np.random.rand()
        if (probability > threshold):
            return True
        else:
            return False
    
    
    def _find_relative_speed(self, index1, index2):
        u_rel = self.particles.u[index1] - self.particles.u[index2]
        v_rel = self.particles.v[index1] - self.particles.v[index2]
        w_rel = self.particles.w[index1] - self.particles.w[index2]
        relative_speed =  np.sqrt(u_rel ** 2  + v_rel ** 2 + w_rel ** 2)
        return (relative_speed)
    
    
    # this function finds the collision area and updates the max collision area.
    def _find_col_area(self, pair, cell_index, relative_speed):
        d_ref = mean(self.particles.dia[pair[0]], self.particles.dia[pair[1]])
        T_ref = mean(self.particles.ref_temp[pair[0]], 
                     self.particles.ref_temp[pair[1]])
        omega_ref = mean(self.particles.visc_index[pair[0]], 
                         self.particles.visc_index[pair[1]])
        
        k = 1.3806488e-23
        constt1 = (2.0 * k * T_ref / self.reduced_mass
                    / (relative_speed ** 2.0)) ** (omega_ref - 0.5)
        constt1 /= math.gamma(2.5 - omega_ref)
        col_area = relative_speed*(0.5*np.pi*d_ref ** 2.0)*constt1			
        
        
        if (self.ref_max_area[cell_index] < col_area):
            self.ref_max_area[cell_index] = col_area
        return col_area



# this class uses different models to collide the given number of particles.
# In other words, this function control which model to use for collision 
# according to the argument. 
class Collider:    
    def __init__(self, particles):
        self.particles = particles
        self.model = {'vhs': VhsCollider(self.particles),
                      1: VhsCollider(self.particles)}

    
    def add_key(self, key, object_ref):
        self.model[key] = object_ref
    
    
    def run(self, collision_pairs, rel_speed, collider_key=None):
        if collider_key == None:
            collider = VhsCollider(self.particles)
            collider.run(collision_pairs, rel_speed)
        else:
            collider = self.model[collider_key]
            collider.run(collision_pairs, rel_speed)
            



class VhsCollider():
    def __init__(self, particles):
        self.particles = particles
        self.collision_pairs = []
        self.rel_speed = []

    
    def run(self, collision_pairs, rel_speed):
        self.collision_pairs = collision_pairs
        self.rel_speed = rel_speed
        if len(self.collision_pairs) > 0:
            self._vhs(0, len(self.collision_pairs) - 1)
            print "collided in the time step"
        else: 
            print "no collision in the time step"
    
    
    def _vhs(self, start, end):
        if (start == end):
            self._vhs_pair(start)
        else:
            mid = int((start + end) / 2)
            self._vhs(start, mid)
            self._vhs(mid + 1, end)
    
    
    def _vhs_pair(self, pair_index):
        pair_index = int(pair_index)
        avg_vel = self._find_pair_avg_vel(pair_index)
        rel_vel = self._find_vhs_post(self.rel_speed[pair_index])
        vel11, vel21 = self._find_velocity(avg_vel[0], rel_vel[0])
        vel12, vel22 = self._find_velocity(avg_vel[1], rel_vel[1])
        vel13, vel23 = self._find_velocity(avg_vel[2], rel_vel[2])
        index1 = self.collision_pairs[pair_index][0]
        index2 = self.collision_pairs[pair_index][1]
        particle = self.particles
        particle.u[index1], particle.v[index1], particle.w[index1] = vel11, vel12, vel13
        particle.u[index2], particle.v[index2], particle.w[index2] = vel21, vel22, vel23
    
    
    def _find_vhs_post(self, rel_speed):
        random_multiplication_factor = 2.0 * np.random.random() - 1.0
        post_u_rel = random_multiplication_factor * rel_speed
        Random_yz_angle = 2.0 * np.pi * np.random.random()
        yz_speed = (1.0 - random_multiplication_factor ** 2)
        yz_speed = np.sqrt(yz_speed) * rel_speed
        post_v_rel = yz_speed * np.sin(Random_yz_angle)
        post_w_rel = yz_speed * np.cos(Random_yz_angle)
        return (post_u_rel, post_v_rel, post_w_rel)
    
    
    def _find_pair_avg_vel(self, index):
        particle = self.particles
        index1 = self.collision_pairs[index][0]
        index2 = self.collision_pairs[index][1]
        u_avg = mean(particle.u[index1], particle.u[index2])
        v_avg = mean(particle.v[index1], particle.v[index2])
        w_avg = mean(particle.w[index1], particle.w[index2])
        return (u_avg, v_avg, w_avg)
    
    
    def _find_velocity(self, avg_vel, rel_vel):
        vel1 = (2.0 * avg_vel + rel_vel) / 2.0
        vel2 = (2.0 * avg_vel - rel_vel) / 2.0
        return (vel1, vel2)