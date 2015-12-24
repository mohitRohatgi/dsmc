# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:12:08 2015

@author: mohit
"""

import numpy as np
import dsmc_detector as dm_d



class MovementManager:
    def __init__(self, particles, surface, surface_temperature=None):
        self.particles = particles
        self.surface = surface
        self.detector = dm_d.IntersectionDetector(surface)
        self.reflector = Reflector(surface, particles, surface_temperature)
    
    def move_all(self, model_key, dt):
        self._move_subset(0, len(self.particles.x) - 1, dt, model_key)
    
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
        point = (self.particles.x[particle_index], self.particles.y[particle_index])
        u = self.particles.u[particle_index]
        v = self.particles.v[particle_index]
        if self.detector.detect_point(point, u, v, dt):
            refl_surface = self.detector.get_surface_index()
            intersect_time = self.detector.get_intersect_time()
            por = self.detector.get_por()
            self.reflector.reflect(particle_index, model_key, refl_surface,
                                       por)
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
    def __init__(self, surface, particles, surface_temperature):
        self.surface = surface
        self.particles = particles
        self.surface_temperature = surface_temperature
        self.model = {'specular' : Specular(surface, particles), 
                      1 : Specular(surface, particles), 
                      'diffuse' : Diffuse(surface, particles,surface_temperature),
                        2 : Diffuse(surface, particles, surface_temperature)}
    
    
    def reflect(self, particle_out, model_key, refl_surface, por):
        self.model[model_key].run(particle_out, refl_surface, por)
    
    
    def add_key(self, key, object_ref):
        self.model[key] = object_ref
    
    
    def show_key(self):
        for key in self.model:
            print key



# this class only changes the velocity and puts the particle at the por.
# from por movement class would take care of movement part.
class Specular:
    def __init__(self, surface, particles):
        self.surface = surface
        self.n_surface_vertices = len(surface)
        self.particles = particles
        self.tangent_finder = dm_d.TangentFinder(surface)
        self.surface_tangent = self.tangent_finder.get__tangent()
        self.refl_surface = 0
    
    
    def run(self, particle_index, surface_index, por):
        self._modify_vel(particle_index, surface_index)
        self._modify_location(particle_index, por, surface_index)
    
    
    def _modify_vel(self, particle_index, surface_index):
        dx = self.surface_tangent[surface_index][0]
        dy = self.surface_tangent[surface_index][1]
        v = self.particles.v[particle_index] * dy * dy
        v -= self.particles.v[particle_index] * dx * dx
        v += self.particles.u[particle_index] * dx * dy * 2.0
        v /= (dx * dx + dy * dy)
        u = self.particles.u[particle_index] * dx * dx
        u -= self.particles.u[particle_index] * dy * dy
        u += self.particles.v[particle_index] * dx * dy * 2.0
        u /= (dx * dx + dy * dy)
        self.particles.v[particle_index] = v
        self.particles.u[particle_index] = u
    
    
    def _modify_location(self, index, por, surface_index):        
        self.particles.x[index] = por[0]
        self.particles.y[index] = por[1]



class Diffuse:
    def __init__(self, surface, particles, surface_temperature):
        self.surface = surface
        self.particles = particles
        self.surface_temperature = surface_temperature
        self.dt = 0.0
        self.particles_out = []
        self.intersect_time = []
        self.refl_surface = []
        self.por = []
    
    
    def run(self, particles_out, dt, model_key, intersect_time, refl_surface,
            por):
        self.dt = dt
        self.particles_out = particles_out
        self.intersect_time = intersect_time
        self.refl_surface = refl_surface
        self.por = por
        self._diffuse_subset(0, len(self.particles_out) - 1)
    
    
    def _diffuse_subset(self, start, end):
        start, end = int(start), int(end)
        if start == end:
            self._specular(start)
        else:
            mid = int((start + end) / 2)
            self._specular_subset(start, mid)
            self._specular_subset(mid + 1, end)
    
    
    def _diffuse(self, index):
        pass