# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:08:11 2015

@author: mohit
"""

import numpy as np


# assuming domain is rectangular and surface is a polygon.
# inlet, outlet and surface need to be wrapped in a list even if they are just 
# represented by a single point.
class Domain:
    def __init__(self, inlet, outlet, surface, volume, surface_temperature):
        self.inlet = inlet
        self.outlet = outlet
        self.surface = surface
        self.volume = volume
        self.surface_temperature = surface_temperature



class Molecules:
    def __init__(self, dia, viscosity_index, mass, viscosity_coeff, dof, tag,
              ref_temperature, gamma, volume, number_density):
        self.dia = dia
        self.number_density = number_density
        self.viscosity_index = viscosity_index
        self.mass = mass
        self.viscosity_coeff = viscosity_coeff
        self.dof = dof
        self.tag = tag
        self.ref_temperature = ref_temperature
        self.gamma = gamma
        self.volume = volume



# mole_fraction should be wraapped up in a list.
class Gas:
    def __init__(self, species, mole_fraction, mach, temperature):
        self.species = species
        self.mole_fraction = mole_fraction
        self.number_density = 0.0
        self.mach_x = mach[0]
        self.mach_y = mach[1]
        self.mach_z = mach[2]
        self.temperature = temperature
        #using tag in substitue of finding length of the array.
        constt = ((species[-1].tag + 2) * (species[-1].tag + 1)) / 2
        self.mean_f_path = np.zeros(species[-1].tag + 1)
        self.mean_col_time = np.zeros(species[-1].tag + 1)
        self.mean_col_rate = np.zeros(constt)
        self.col_density_rate = np.zeros(constt)
        self.dia = np.zeros(constt)
        self.ref_temperature = np.zeros(constt)
        self.mpv = np.zeros(len(species))
        self.reduced_mass = 0.0
        self.viscosity_index = np.zeros(constt)
        self.gamma = 0.0
        self.mass = 0.0
    
    
    def setup(self):
        self._find_number_density()
        self._find_gas_property()
        self._find_col_property()
    
    
    def _find_number_density(self):
        self.number_density = 0.0
        for index, species in enumerate(self.species):
            self.number_density += species.number_density
    
    
    def _find_gas_property(self):
        self.gamma = 0.0
        self.mass = 0.0
        prod = 1.0
        count = 0.0
        k = 1.3806488e-23
        for index,species in enumerate(self.species):
            self.gamma += species.gamma * self.mole_fraction[index]
            self.mass += species.mass * self.mole_fraction[index]
            self.mpv[index] = np.sqrt(2.0 * k * self.temperature 
                                            / species.mass)
            if len(self.species) > 1:
                prod *= species.mass
                count += species.mass
        if len(self.species) > 1:
            self.reduced_mass = prod / count
        else:
            self.reduced_mass = 0.5 * self.species[0].mass
#        print "reduced = ", self.reduced_mass
    
    
    # bug here.
    def _find_col_property(self):
        # boltzmann constant
        k = 1.3806488e-23
        for index1, molecule1 in enumerate(self.species):
            for index2 in range(index1, self.species[-1].tag + 1):
                molecule2 = self.species[index2]
                index = ((index1 * (index1 + 1)) / 2)
                index += index2
                self.ref_temperature[index] = self.species[index2].ref_temperature
                self.dia[index] = 0.5 * (molecule1.dia + molecule2.dia)
                self.viscosity_index[index] = 0.5 * (molecule1.viscosity_index 
                                            + molecule2.viscosity_index)
                self.ref_temperature[index] = 0.5 * (molecule1.ref_temperature 
                                            + molecule2.ref_temperature)
            constt = 0.0
            for index2, molecule2 in enumerate(self.species):
                if index1 < index2:
                    index = ((index1 * (index1 + 1)) / 2)
                    index += index2
                else:
                    index = ((index2 * (index2 + 1)) / 2)
                    index += index1
                self.mean_f_path[index1] += ((1 + molecule1.mass / molecule2.mass)
                        ** 0.5 * (self.ref_temperature[index] / self.temperature)
                        ** (self.viscosity_index[index] - 0.5) * self.dia[index]
                        ** 2.0 * self.mole_fraction[index2])
                
                constt = (self.dia[index] ** 2.0 * self.mole_fraction[index2]
                        * (self.temperature / self.ref_temperature[index]) ** 
                        (1.0 - self.viscosity_index[index]) * 
                        self.ref_temperature[index] ** 0.5)
                
                self.mean_col_rate[index1] += constt
                
                if index1 != index2:
                    self.col_density_rate[index] = self.mole_fraction[index1] * constt
                else:
                    self.col_density_rate[index] = self.mole_fraction[index1] * constt * 0.5
        
            self.mean_f_path[index1] = 1 / self.mean_f_path[index1]
#            print self.mean_f_path[index1]
            self.mean_f_path[index1] /= (np.pi * self.number_density)
            self.mean_col_rate[index1] *= 2.0 * self.number_density * np.sqrt(2.0 * 
                                        k * np.pi / self.reduced_mass)
            self.mean_col_time[index1] = 1 / self.mean_col_rate[index1]
#            print "col_rate = ", self.mean_col_rate[index1]


# tag represents the species type of particle.
# particles are stored in an ascending order of the species.
# for eg. species 1 would be stored in the first n1 particles the species2 for
# the next n2 particles and so on.
class Particles:
    def __init__(self, n_particles, gas):
        self.gas = gas
        self.x = np.zeros(n_particles)
        self.y = np.zeros(n_particles)
        self.u = np.zeros(n_particles)
        self.v = np.zeros(n_particles)
        self.w = np.zeros(n_particles)
        self.eu = np.zeros(n_particles)
        self.ev = np.zeros(n_particles)
        self.ew = np.zeros(n_particles)
        self.mass = np.ones(n_particles)
        self.dia = np.ones(n_particles)
        self.cross_area = np.ones(n_particles)
        self.viscosity_index = np.ones(n_particles)
        self.ref_temperature = np.ones(n_particles)
        self.viscosity_coeff = np.ones(n_particles)
        self.dof = np.ones(n_particles)
        self.gamma = np.ones(n_particles)
        self.mpv = np.ones(n_particles)
        self.num = n_particles
        self.n_eff = 0.0
        self.tag = np.zeros(n_particles, dtype = int)
    
    
    # this function should only be called once.
    # sets the particles according to the molecules it represents.
    def setup(self, domain_volume):
        n_species = self.num * np.asarray(self.gas.mole_fraction)
        count = 0
        self.n_eff = self.gas.number_density * domain_volume / self.num
        for index, molecule in enumerate(self.gas.species):
            try:
                end = n_species[index + 1]
            except:
                end = self.num
            self.mass[count:end] = molecule.mass
            self.dia[count:end] = molecule.dia
            self.cross_area[count:end] = 0.25 * np.pi * molecule.dia ** 2.0
            self.viscosity_index[count:end] = molecule.viscosity_index
            self.ref_temperature[count:end] = molecule.ref_temperature
            self.viscosity_coeff[count:end] = molecule.viscosity_coeff
            self.dof[count:end] = molecule.dof
            self.gamma[count:end] = molecule.gamma
            self.mpv[count:end] = self.gas.mpv[index]
            self.tag[count:end] = index
            count += n_species[index]
         
    
    def move(self, dt, index_list):
        for index in index_list:
            self.x[index] += self.u[index] * dt
            self.y[index] += self.v[index] * dt
    
    
    def randomize_location(self):
        self.x = np.random.rand(self.n_particles)
        self.y = np.random.rand(self.n_particles)
    
    
    def compute_energy(self):
        self.eu = self.u * self.u * self.mass
        self.ev = self.v * self.v * self.mass
        self.ew = self.w * self.w * self.mass