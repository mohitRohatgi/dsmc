# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:12:08 2015

@author: mohit
"""

import numpy as np
import dsmc_detector as dm_d


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
    
    
    def run(self, particles_out, dt, model_key, intersect_time, refl_surface,
            por):
        self.model[model_key].run(particles_out, dt, intersect_time, 
                                  refl_surface, por)
    
    
    def add_key(self, key, object_ref):
        self.model[key] = object_ref
    
    
    def show_key(self):
        for key in self.model:
            print key



class Specular:
    def __init__(self, surface, particles):
        self.surface = surface
        self.particles = particles
        self.dt = 0.0
        self.locator = dm_d.SurfaceLocator(surface)
        self.surface_slope = self.locator.get_slope()
        self.particles_out = []
        self.intersect_time = []
        self.refl_surface = []
        self.por = []
    
    
    def run(self, particles_out, dt, intersect_time, refl_surface, por):
        self.dt = dt
        self.particles_out = particles_out
        self.intersect_time = intersect_time
        self.refl_surface = refl_surface
        self.por = por
        self._specular_subset(0, len(self.particles_out) - 1)
    
    
    def _specular_subset(self, start, end):
        start, end = int(start), int(end)
        if start == end:
            self._specular(start)
        else:
            mid = int((start + end) / 2)
            self._specular_subset(start, mid)
            self._specular_subset(mid + 1, end)
    
    
    def _specular(self, index):
        self._modify_vel(index)
        self._modify_location(index)
    
    
    def _modify_vel(self, index):
        tan = self.surface_slope[self.refl_surface[index]]
        u = 2.0 * self.particles.v[index] * tan
        u += self.particles.u[index] * (1.0 - tan * tan)
        u /= (1 + tan * tan)
        v = 2.0 * self.particles.u[index] * tan
        v += self.particles.v[index] * (1.0 - tan * tan)
        v /= (1 + tan * tan)
        self.particles.u[index] = u
        self.particles.v[index] = v
    
    
    def _modify_location(self, index):
        dt = self.dt - self.intersect_time[index]
        self.particles.x[index] = self.por[index][0] + self.particles.u[index] * dt
        self.particles.y[index] = self.por[index][1] + self.particles.v[index] * dt



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