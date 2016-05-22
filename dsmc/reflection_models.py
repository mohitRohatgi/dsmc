# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 13:00:05 2016

@author: mohit
"""

import numpy as np

# this class only changes the velocity and puts the particle at the por.
# from por movement class would take care of movement part.
class Specular:
    def __init__(self, surf_group, particles):
        self.surf_group = surf_group
        self.particles = particles
    
    
    def reflect(self, particle_out, group_index, surf_index, por):
        self._run(particle_out, group_index, surf_index)
        self._modify_location(particle_out, por)
    
    
    def _run(self, particle_index, group_index, surf_index):
        tangent = self.surf_group.get_surf_tangent(group_index, surf_index)
        dx = tangent[0]
        dy = tangent[1]
        v = self.particles.get_vely(particle_index) * dy * dy
        v -= self.particles.get_vely(particle_index) * dx * dx
        v += self.particles.get_velx(particle_index) * dx * dy * 2.0
        u = self.particles.get_velx(particle_index) * dx * dx
        u -= self.particles.get_velx(particle_index) * dy * dy
        u += self.particles.get_vely(particle_index) * dx * dy * 2.0
        self.particles.set_velx(u, particle_index)
        self.particles.set_vely(v, particle_index)
    
    
    def _modify_location(self, particle_out, por):        
        self.particles.set_x(por[0], particle_out)
        self.particles.set_y(por[1], particle_out)



class Diffuse:
    def __init__(self, surf_group, particles):
        self.surf_group = surf_group
        self.particles = particles
        self.k = 1.3806488e-23
    
    
    def reflect(self, particle_out, group_index, surf_index, por):
        self._run(particle_out, group_index, surf_index)
        self._modify_location(particle_out, por)
    
    
    def _run(self, particle_index, group_index, surf_index):
        T = self.surf_group.get_surf_temp(group_index, surf_index)
        mpv = np.sqrt(2.0 * self.k * T / self.particles.get_mass(particle_index))
        normal = self.surf_group.get_surf_normal(group_index, surf_index)
        tangent = self.surf_group.get_surf_tangent(group_index, surf_index)
        c1 = np.random.normal(0.0, 1.0) * mpv
        c2 = np.random.normal(0.0, 1.0) * mpv
        theta = np.random.random() * 2.0 * np.pi
        self.particles.set_velz(c2 * np.sin(theta), particle_index)
        c2 *= np.cos(theta)
        u = (c1 * tangent[0] + c2 * normal[0])
        v = (c1 * tangent[1] + c2 * normal[1])
        self.particles.set_velx(u, particle_index)
        self.particles.set_vely(v, particle_index)
    
    
    def _modify_location(self, particle_out, por):        
        self.particles.set_x(por[0], particle_out)
        self.particles.set_y(por[1], particle_out)
