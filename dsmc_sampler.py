# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:03:57 2015

@author: mohit
"""

import numpy as np


# p_time is the present time of the simulation while time0 is the initial
# f_time is the final time.
class Sampling_manager:
    def __init__(self, cells, cells_in, gas, particles, n_steps):
        self.cells = cells
        self.particles = particles
        self.p_time = 0.0
        self.n_steps = int(n_steps)
        self.p_step = 0
        self.instant_sampler = Instant_sampler(cells, cells_in, particles, gas)
        self.time_sampler = Time_sampler(cells, cells_in, particles, n_steps)
    
    
    def run(self):
        self.instant_sampler.run()
        self.time_sampler.sample_domain()
        self.p_step += 1
    
    
    def get_temperature(self):
        return self.time_sampler.temperature
    
    
    def get_number_density(self):
        return self.time_sampler.number_density



class Time_sampler:
    def __init__(self, cells, cells_in, particles, n_steps):
        self.cells = cells
        self.cells_in = cells_in
        self.particles = particles
        self.number_density = np.zeros((cells.n_x * cells.n_y, 
                                        len(cells.n_species[0])))
        self.temperature = np.zeros(cells.n_x * cells.n_y)
        self.n_steps = n_steps
        self.step_counter = 0
    
    
    def sample_domain(self):
        self.step_counter += 1
        for index in self.cells_in:
            self._sample_cell(index)
    
    
    def _sample_cell(self, cell_index):
        total = sum (self.cells.n_species[cell_index])
        if total > 0.0:
            self.temperature[cell_index] += self.cells.temperature[cell_index]
            self.number_density[cell_index] += (self.cells.number_density[cell_index])
        
        if (self.step_counter >= self.n_steps):
            self._output(cell_index)
    
    
    def _output(self, cell_index):
        self.temperature[cell_index] /= self.step_counter
        self.number_density[cell_index] /= self.step_counter



# this class finds temperature and number density of function at that time
# step, i.e., instataneous sampling is done.
class Instant_sampler:
    def __init__ (self, cells, cells_in, particles, gas):
        self.cells = cells
        self.particles = particles
        self.gas = gas
        self.cells_in = cells_in
        # n_species represents the no. of species in the gas.
        self.n_species = len(gas.species)
    
    
    # this function samples all the cell at a paricular time step.
    def run(self):
        self.particles.compute_energy()
        self._find_properties()
        self._sample_subset(0, len(self.cells.temperature) - 1)
    
    
    def _find_properties(self):
        for index in self.cells_in:
            self._find_property(index)
    
    
    def _find_property(self, cell_index):
        self.cells.u[cell_index] = 0.0
        self.cells.v[cell_index] = 0.0
        self.cells.w[cell_index] = 0.0
        self.cells.mass[cell_index] = 0.0
        count = 0
        self.cells.number_density[cell_index] = len(self.cells.particles_inside[cell_index])
        if len(self.cells.particles_inside[cell_index]):
            for index in self.cells.particles_inside[cell_index]:
                self.cells.u[cell_index] += self.particles.u[index] * self.particles.mass[index]
                self.cells.v[cell_index] += self.particles.v[index] * self.particles.mass[index]
                self.cells.w[cell_index] += self.particles.w[index] * self.particles.mass[index]
                self.cells.mass[cell_index] += self.particles.mass[index]
                count += 1
            self.cells.u[cell_index] /= self.cells.mass[cell_index]
            self.cells.v[cell_index] /= self.cells.mass[cell_index]
            self.cells.w[cell_index] /= self.cells.mass[cell_index]
            self.cells.mass[cell_index] /= count
    
    
    def _sample_subset(self, start, end):
        start = int(start)
        end = int(end)
        if start == end:
            self._sample_cell(start)
        else:
            mid = int((start + end) / 2)
            self._sample_subset(start, mid)
            self._sample_subset(mid + 1, end)
    
    
    def _sample_cell(self, cell_index):
        self._find_temperature(cell_index)
        self._find_number_density(cell_index)
    
    
    def _find_temperature(self, cell_index):
        k = 1.3806488e-23
        constt = len(self.cells.particles_inside[cell_index])
        self.cells.temperature[cell_index] = 0.0
        if constt > 0.0:
            energy = self._find_cell_energy(cell_index)
            vel_energy = self._find_vel_energy(cell_index)
            self.cells.temperature[cell_index] = ((energy - vel_energy) * 2.0
                                                    / 3.0 / k)
    
    def _find_cell_energy(self, cell_index):
        energy = 0.0
        for index in self.cells.particles_inside[cell_index]:
            energy += (self.particles.eu[index] + self.particles.ev[index] + 
                        self.particles.ew[index])
        return energy / len((self.cells.particles_inside[cell_index]))
    
    
    def _find_vel_energy(self, cell_index):
        u = self.cells.u[cell_index]
        v = self.cells.v[cell_index]
        w = self.cells.w[cell_index]
        mass = self.cells.mass[cell_index]
        return (u ** 2.0 + v ** 2.0 + w ** 2.0) * mass
    
    
    def _find_number_density(self, index):
        self.cells.number_density[index] = (sum(self.cells.n_species[index]) / 
                                            self.cells.volume[index])