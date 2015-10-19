# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:12:08 2015

@author: mohit
"""

import numpy as np
import dsmc_detector as dm_d



class MovementManager:
    def __init__(self, particles, surface, surface_temperature):
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
    
    
    def _reflect_n_move(self, particle_index, dt, model_key):
        particle_index = int(particle_index)
        point = (self.particles.x[particle_index], self.particles.y[particle_index])
        u = self.particles.u[particle_index]
        v = self.particles.v[particle_index]
        if self.detector.detect_point(point, u, v, dt):
            refl_surface = self.detector.get_surface_index()
            intersect_time = self.detector.get_intersect_time()
            por = self.detector.get_por()
            self.reflector.reflect(particle_index, model_key, refl_surface, por)
#            remaining_time = dt - intersect_time
#            self._reflect_n_move(particle_index, remaining_time, model_key)
        else:
            self.particles.move(particle_index, dt)
            
            
#            if (self.particles.x[particle_index] > 1.0 or self.particles.y[particle_index] > 1.0
#            or self.particles.x[particle_index] < 0.0 or self.particles.y[particle_index] < 0.0):
#                print "no reflection = ", particle_index, self.particles.x[particle_index], self.particles.y[particle_index]



# por stands for the point of reflection of particle from the surface and is 
# stored according to the index of the particle.
class Reflector:
    def __init__(self, surface, particles, surface_temperature):
        self.surface = surface
        self.particles = particles
        self.surface_temperature = surface_temperature
        self.model = {'specular' : Specular(surface, particles), 1 : Specular(
                    surface, particles), 'diffuse' : Diffuse(surface, particles, 
                    surface_temperature), 2 : Diffuse(surface, particles, 
                    surface_temperature)}
    
    
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
        self.particles = particles
        self.normal_finder = dm_d.NormalFinder(surface)
        self.surface_tangent = self.normal_finder.get__tangent()
        self.refl_surface = 0
    
    
    def run(self, particle_index, surface_index, por):
        self._modify_vel(particle_index, surface_index)
        self._modify_location(particle_index, por)
    
    
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
        u /= dx * dx + dy * dy
        self.particles.v[particle_index] = v
        self.particles.u[particle_index] = u
    
    
    def _modify_location(self, index, por):
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


# Reflector and Reflection_detector classes should be visible to this class.
# for adding new model
class ReflectionManager:
    def __init__(self, surface, particles, surface_temperature):
        self.detector = ReflectionDetector(surface, particles)
        self.particles = particles
        self.reflector = Reflector(surface, particles, surface_temperature)
    
    
    def run(self, model, dt):
        self.detector.detect_all(dt)
        particles_out = self.detector.get_particles_out()
        por = self.detector.get_por()
        intersect_time = self.detector.get_intersect_time()
        refl_surface = self.detector.get_surface()
        if len(particles_out) > 0:
            self.reflector.run(particles_out, dt, model, intersect_time, 
                               refl_surface, por)
#            print "re-checking for particles_out"
#            self.run(model, dt)
        return particles_out
    
    
    # after adding models in Reflector this function should be updated.
    def add_model(self, key, obj_ref):
        self.reflector.add_key(key, obj_ref)



# this class detects if reflect takes place, and if it does, which particle
# and their respective surface.
class ReflectionDetector:
    # considering surface to be a continuous union of different linear segment
    # It is also assumed that it contain vertices in an ascending order stating 
    # from the leftmost vertex as its first element.
    # ref_point is an arbitary point inside domain and not on surface or out of
    # domain.
    # ref_sign is the sign of ref_point w.r.t. surface segments. For particles
    # out of domain, the product of ref_sign and its sign should be negative.
    # 1 stands for point is in increasing y-direction and -1 for decreasing 
    # y-direction.
    # setup needs to be called after each time step.
    def __init__(self, surface, particles):
        self.particles = particles
        self.surface = surface
        self.detector = dm_d.IntersectionDetector(surface)
        self.intersect_time = []
        self.por = []
        self.refl_surface = []
        self.particles_out = []
    
    
    def get_surface(self):
        return self.refl_surface
    
    
    def get_por(self):
        return self.por
    
    
    def get_intersect_time(self):
        return self.intersect_time
    
    
    def get_particles_out(self):
        return self.particles_out
    
    
    # this function would find particles crossing the surface.
    def detect_all(self, dt):
        points = []
        slopes = []
        vels = []
        for index in range(len(self.particles.x)):
            point = (self.particles.x[index], self.particles.y[index])
            slope = self.particles.v[index] / self.particles.u[index]
            vel = np.sqrt(self.particles.u[index] ** 2.0 + 
                            self.particles.v[index] ** 2.0)
            points.append(point)
            slopes.append(slope)
            vels.append(vel)
        self.detector.detect_all(points, slopes, vels, dt)
        self._extract()
        return True
    
    
    def _extract(self):
        self.intersect_time = self.detector.get_intersect_time()
        self.por = self.detector.get_por()
        self.refl_surface = self.detector.get_surface()
        self.particles_out = self.detector.get_points()