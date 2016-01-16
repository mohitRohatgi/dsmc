# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:12:08 2015

@author: mohit
"""

import dsmc_detector as dm_d
import numpy as np



class MovementManager:
    def __init__(self, particles, surf_group, surface_temperature=0):
        self.particles = particles
        self.surf_group = surf_group
        self.detector = dm_d.IntersectionDetector(surf_group)
        self.reflector = Reflector(surf_group, particles)
    
    def move_all(self, model_key, dt):
        self._move_subset(0, self.particles.get_particles_count() - 1, dt, model_key)
    
    def _move_subset(self, start, end, dt, model_key):
        start, end = int(start), int(end)
        if start == end:
            self._reflect_n_move(start, dt, model_key)
        else:
            mid = int((start + end) / 2)
            self._move_subset(start, mid, dt, model_key)
            self._move_subset(mid + 1, end, dt, model_key)
    
    
    # this function checks if the particle would reflect and correspondingly
    # would either reflect or move the particle.
    def _reflect_n_move(self, particle_index, dt, model_key):
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
            self.reflector.reflect(particle_index, model_key, group_index,
                                   surf_index, por)
            self.particles.move(particle_index, dt * 0.01)
            remaining_time = dt * 0.99 - intersect_time
            if remaining_time < dt * 0.01:
                self.particles.move(particle_index, dt * 0.01)
            else:
                self._reflect_n_move(particle_index, remaining_time, model_key)
        else:
            self.particles.move(particle_index, dt)


# por stands for the point of reflection of particle from the surface and is 
# stored according to the index of the particle.
class Reflector:
    def __init__(self, surf_group, particles):
        self.surface = surf_group
        self.particles = particles
        self.model = {'specular' : Specular(surf_group, particles), 
                      1 : Specular(surf_group, particles), 
                      'diffuse' : Diffuse(surf_group, particles),
                        2 : Diffuse(surf_group, particles)}
    
    
    def reflect(self, particle_out, model_key, group_index, surf_index, por):
        self.model[model_key].run(particle_out, group_index, surf_index)
        self._modify_location(particle_out, por)
    
    
    def _modify_location(self, particle_out, por):        
        self.particles.set_x(por[0], particle_out)
        self.particles.set_y(por[1], particle_out)



# this class only changes the velocity and puts the particle at the por.
# from por movement class would take care of movement part.
class Specular:
    def __init__(self, surf_group, particles):
        self.surf_group = surf_group
        self.particles = particles    
    
    def run(self, particle_index, group_index, surf_index):
        tangent = self.surf_group.get_surf_tangent(group_index, surf_index)
        dx = tangent[0]
        dy = tangent[1]
        v = self.particles.get_vely(particle_index) * dy * dy
        v -= self.particles.get_vely(particle_index) * dx * dx
        v += self.particles.get_velx(particle_index) * dx * dy * 2.0
        v /= (dx * dx + dy * dy)
        u = self.particles.get_velx(particle_index) * dx * dx
        u -= self.particles.get_velx(particle_index) * dy * dy
        u += self.particles.get_vely(particle_index) * dx * dy * 2.0
        u /= (dx * dx + dy * dy)
        self.particles.set_velx(u, particle_index)
        self.particles.set_vely(v, particle_index)



class Diffuse:
    def __init__(self, surf_group, particles):
        self.surface = surf_group
        self.particles = particles
    
    
    def run(self, particle_index, group_index, surface_index):
        mpv = np.sqrt(1.0 / self.particles.get_mass(particle_index))
        normal = self._find_normal(surface_index)
        c1 = np.random.normal(0.0, 0.5) * mpv
        c2 = np.random.normal(0.0, 0.5) * mpv
        theta = np.random.random() * 2.0 * np.pi
        self.particles.set_velz(c2 * np.sin(theta))
        c2 *= np.cos(theta)
        s = np.sqrt(normal[0] ** 2.0 + normal[1] ** 2.0)
        u = c1 * self.surface_tangent[surface_index][0] + c2 * normal[0]
        v = c1 * self.surface_tangent[surface_index][1] + c2 * normal[1]
        self.particles.set_velx(u / s)
        self.particles.set_vely(v / s)
    
    
    def _find_normal(self, surface_index):
        pass