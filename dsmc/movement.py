# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:12:08 2015

@author: mohit
"""

import detector as dm_d

class MovementManager:
    def __init__(self, particles, surf_group, Reflector):
        self.particles = particles
        self.detector = dm_d.IntersectionDetector(surf_group)
        self.reflector = Reflector(surf_group, particles)
    
    def move_all(self, dt):
        for i in range(self.particles.get_particles_count()):
            self._reflect_n_move(i, dt)
    
    
    # this function checks if the particle would reflect and correspondingly
    # would either reflect or move the particle.
    def _reflect_n_move(self, particle_index, dt):
        particle_index = int(particle_index)
        point = (self.particles.get_x(particle_index), 
                 self.particles.get_y(particle_index))
        u = self.particles.get_velx(particle_index)
        v = self.particles.get_vely(particle_index)
        if self.detector.detect_point(point, u, v, dt):
            group_index = self.detector.get_surf_group_index()
            surf_index = self.detector.get_surface_index()
            intersect_time = self.detector.get_intersect_time()
            por = self.detector.get_por()
            self.reflector.reflect(particle_index, group_index, surf_index, por)
            self.particles.move(particle_index, dt * 0.01)
            remaining_time = dt * 0.99 - intersect_time
            if remaining_time < dt * 0.01:
                self.particles.move(particle_index, dt * 0.01)
            else:
                self._reflect_n_move(particle_index, remaining_time)
        else:
            self.particles.move(particle_index, dt)