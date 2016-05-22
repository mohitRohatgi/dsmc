# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:03:57 2015

@author: mohit
"""

import numpy as np


# p_time is the present time of the simulation while time0 is the initial
# f_time is the final time.
class SamplingManager:
    def __init__(self, cells, cells_in, n_species, 
                 particles, n_steps, ignore_frac):
        self.cells = cells
        self.particles = particles
        self.n_steps = int(n_steps)
        self.instant_sampler = Instant_sampler(cells, cells_in, particles,
                                               n_species)
        self.time_sampler = TimeSampler(cells, cells_in, particles, n_steps,
                                         n_species, ignore_frac)
    
    
    def run(self):
        self.instant_sampler.run()
        self.time_sampler.sample_domain()
    
    
    def get_temperature(self):
        return self.time_sampler.get_temperature()
    
    
    def get_mach(self):
        return self.time_sampler.get_mach()
    
    
    def get_all_number_density(self):
        return self.time_sampler.get_all_number_density()
    
    
    def get_number_density(self, tag):
        return self.time_sampler.get_number_density(tag)



class TimeSampler:
    def __init__(self, cells, cells_in, particles, 
                 n_steps, n_species, ignore_frac):
        self.cells = cells
        self.cells_in = cells_in
        self.particles = particles
        self.number_density = np.zeros((n_species, len(cells.get_temperature())))
        self.temperature = np.zeros(len(cells.get_temperature()))
        self.u = np.zeros_like(self.temperature)
        self.v = np.zeros_like(self.temperature)
        self.w = np.zeros_like(self.temperature)
        self.gamma = np.zeros_like(self.temperature)
        self.mass = np.zeros_like(self.temperature)
        self.mach = np.zeros_like(self.temperature)
        self.n_steps = n_steps
        self.step_counter = 0
        self.n_species = n_species
        self.ignore_frac = ignore_frac
    
    
    def get_all_number_density(self):
        return self.number_density
    
    
    def get_number_density(self, species_index=0):
        return self.number_density[species_index]
    
    
    def get_temperature(self):
        return self.temperature
    
    
    def get_mach(self):
        return self.u
    
    
    def sample_domain(self):
        self.step_counter += 1
        if (self.step_counter >= self.n_steps * self.ignore_frac):
            self._sample_cell()
            
            if (self.step_counter >= self.n_steps):
                self._output()
    
    
    def _sample_cell(self):
        self.u += self.cells.get_velx()
        self.v += self.cells.get_vely()
        self.w += self.cells.get_velz()
        self.mass += self.cells.get_mass()
        self.gamma += self.cells.get_gamma()
        self.temperature += self.cells.get_temperature()
        self.number_density += self.cells.get_number_density()
    
    
    def _output(self):
        k = 1.3806488e-23
        sample_size = self.n_steps * (1.0 - self.ignore_frac)
        self.gamma /= sample_size
        self.mass /= sample_size
        self.u /= sample_size
        self.v /= sample_size
        self.w /= sample_size
        self.temperature /= sample_size
        self.number_density /= sample_size
        sound_speed_sq = self.gamma * k * self.temperature / self.mass
        self.mach = self.u * self.u + self.v * self.v + self.w * self.w
        self.mach /= sound_speed_sq
        self.mach = np.sqrt(self.mach)
        
        self.mach[self.temperature < 1.0e-6] = 0.0
        
        



# this class finds temperature and number density of function at that time
# step, i.e., instataneous sampling is done.
class Instant_sampler:
    def __init__ (self, cells, cells_in, particles, n_species):
        self.cells = cells
        self.particles = particles
        self.cells_in = cells_in
        # n_species represents the no. of species in the gas.
        self.n_species = n_species
    
    
    # this function samples all the cell at a paricular time step.
    def run(self):
        self.particles.compute_energy()
        self._find_properties()
        self._sample_domain()
    
    
    def _find_properties(self):
        for index in self.cells_in:
            self._set_property(index)
    
    
    def _set_property(self, cell_index):
        u, v, w, mass, gamma = 0.0, 0.0, 0.0, 0.0, 0.0
        
        if len(self.cells.get_particles_inside(cell_index)) > 0:
            for index in self.cells.get_particles_inside(cell_index):
                u += self.particles.get_velx(index) * self.particles.get_mass(index)
                v += self.particles.get_vely(index) * self.particles.get_mass(index)
                w += self.particles.get_velz(index) * self.particles.get_mass(index)
                mass += self.particles.get_mass(index)
                gamma += self.particles.get_gamma(index)
            
            u /= mass
            v /= mass
            w /= mass
            mass /= len(self.cells.get_particles_inside(cell_index))
            gamma /= len(self.cells.get_particles_inside(cell_index))
            
        self.cells.set_velx(u, cell_index)
        self.cells.set_vely(v, cell_index)
        self.cells.set_velz(w, cell_index)
        self.cells.set_mass(mass, cell_index)
        self.cells.set_gamma(gamma, cell_index)
    
    
    def _sample_domain(self):
        for cell_index in self.cells_in:
            self._find_temperature(cell_index)
            self._find_number_density(cell_index)
    
    
    def _find_temperature(self, cell_index):
        k = 1.3806488e-23
        if len(self.cells.get_particles_inside(cell_index)) > 1:
            energy = self._find_cell_energy(cell_index)
            vel_energy = self._find_vel_energy(cell_index)
            temperature = (energy - vel_energy) * 2.0 / 3.0 / k
            self.cells.set_temperature(temperature, cell_index)
        else:
            self.cells.set_temperature(0.0, cell_index)
    
    def _find_cell_energy(self, cell_index):
        energy = 0.0
        for index in self.cells.get_particles_inside(cell_index):
            energy += (self.particles.get_eu(index) + self.particles.get_ev(
                        index) + self.particles.get_ew(index))
        return energy / len(self.cells.get_particles_inside(cell_index))
    
    
    def _find_vel_energy(self, cell_index):
        u = self.cells.get_velx(cell_index)
        v = self.cells.get_vely(cell_index)
        w = self.cells.get_velz(cell_index)
        return (u * u + v * v + w * w) * self.cells.get_mass(cell_index)
    
    
    def _find_number_density(self, index):
        for tag in range(self.n_species):
            constt = (self.cells.get_n_particles(tag, index) * 
                self.particles.get_n_eff()) / self.cells.get_cell_volume()
            
            self.cells.set_number_density(constt, tag, index)
