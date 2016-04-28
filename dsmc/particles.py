# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:08:11 2015

@author: mohit
"""

import numpy as np


class Molecules:
    def __init__(self, dia, visc_index, mass, visc_coeff, dof, tag,
              ref_temp, gamma, volume, number_density):
        self.dia = dia
        self.number_density = number_density
        self.visc_index = visc_index
        self.mass = mass
        self.visc_coeff = visc_coeff
        self.dof = dof
        self.tag = tag
        self.ref_temp = ref_temp
        self.gamma = gamma
        self.volume = volume
    
    
    def get_dia(self):
        return self.dia
    
    
    def get_visc_index(self):
        return self.visc_index
    
    
    def get_ref_temp(self):
        return self.ref_temp
    
    
    def get_tag(self):
        return self.tag
    
    
    def get_number_density(self):
        return self.number_density
    
    
    def get_gamma(self):
        return self.gamma
    
    
    def get_mass(self):
        return self.mass
    
    
    def get_dof(self):
        return self.dof
    
    
    def get_visc_coeff(self):
        return self.visc_coeff



# mol_frac should be wraapped up in a list.
# species should organised in ascending order w.r.t. tag.
# species represents the list of molecules.
class Gas:
    def __init__(self, species, mol_frac, mach, temperature):
        self.species = species
        self.mol_frac = mol_frac
        self.number_density = 0.0
        self.mach = mach
        self.temperature = temperature
        #using tag in substitue of finding length of the array.
        self.mpv = np.zeros(len(species))
        self.reduced_mass = 0.0
        self.gamma = 0.0
        self.mass = 0.0
    
    
    def setup(self):
        # finding number density ...
        for index, species in enumerate(self.species):
            self.number_density += species.get_number_density()
        
        # finding gas properties ...
        prod = 1.0
        count = 0.0
        k = 1.3806488e-23
        for index,species in enumerate(self.species):
            self.gamma += species.get_gamma() * self.mol_frac[index]
            self.mass += species.get_mass() * self.mol_frac[index]
            self.mpv[index] = np.sqrt(2.0 * k * self.temperature 
                                            / species.get_mass())
            if len(self.species) > 1:
                prod *= species.get_mass()
                count += species.get_mass()
        if len(self.species) > 1:
            self.reduced_mass = prod / count
        else:
            self.reduced_mass = 0.5 * self.species[0].get_mass()
    
    
    # tag1, tag2 are the tags of the molecules, i.e., it represents the species
    # of the molecule the particle represents.
#    def get_gamma_function(self):        
#        return self.gamma_func
    
    
    def get_mass(self):
        return self.mass
    
    
    def get_number_density(self):
        return self.number_density
    
    
    def get_temperature(self):
        return self.temperature
    
    
    def get_mol_frac(self):
        return self.mol_frac
    
    
    def get_species(self):
        return self.species
    
    
    def get_n_species(self):
        return len(self.species)
    
    
    def get_reduced_mass(self):
        return self.reduced_mass
    
    
    def get_mach_x(self):
        return self.mach[0]
    
    
    def get_mach_y(self):
        return self.mach[1]
    
    
    def get_mach_z(self):
        return self.mach[2]
    
    
    def get_gamma(self):
        return self.gamma
    
    
    def get_mpv(self):
        return self.mpv


# tag represents the species type of particle.
# particles are stored in an ascending order of the species.
# for eg. species 1 would be stored in the first n1 particles the species2 for
# the next n2 particles and so on.
# tag starts from zero.
class Particles:
    def __init__(self, n_particles, n_eff=0.0):
        self.num = n_particles
        self.x = np.zeros(n_particles)
        self.y = np.zeros(n_particles)
        self.u = np.zeros(n_particles)
        self.v = np.zeros(n_particles)
        self.w = np.zeros(n_particles)
        self.eu = np.zeros(n_particles)
        self.ev = np.zeros(n_particles)
        self.ew = np.zeros(n_particles)
        self.n_eff = n_eff
        self.tag = np.zeros(n_particles, dtype = int)
        self.mpv = np.zeros(n_particles)
        self.species = []


    # this function should only be called once.
    # sets the particles according to the molecules it represents.
    # species is a list of molecule instances.
    def setup(self, mole_fraction, species, mpv):
        self.species = species
        n_species = len(species)
        for index in range(self.num):
            tag = index % n_species
            self.tag[index] = tag
            self.mpv[index] = mpv[tag]


    def move_all(self, dt):
        self.x += self.u * dt
        self.y += self.v * dt
    
    
    def move_sub(self, dt, particle_list):
        for index in particle_list:
            self.move(index, dt)
         
    
    def move(self, index, dt):
        self.x[index] += self.u[index] * dt
        self.y[index] += self.v[index] * dt
    
    
    def randomize_location(self):
        self.x = np.random.rand(self.n_particles)
        self.y = np.random.rand(self.n_particles)
    
    
    def compute_energy(self):
        for index in range(len(self.x)):
            self.set_particle_energy(index)
    
    
    def set_particle_energy(self, index):
        mass = self.get_mass(index)
        self.eu[index] = self.u[index] * self.u[index] * mass
        self.ev[index] = self.v[index] * self.v[index] * mass
        self.ew[index] = self.w[index] * self.w[index] * mass
    
    
    def get_tag(self, index):
        return self.tag[index]
    
    
    def get_dia(self, index):
        return self.species[self.get_tag(index)].get_dia()
    
    
    def get_gamma(self, index):
        return self.species[self.get_tag(index)].get_gamma()
    
    
    def get_mass(self, index):
        return self.species[self.get_tag(index)].get_mass()
    
    
    def get_ref_temp(self, index):
        return self.species[self.get_tag(index)].get_ref_temp()
    
    
    def get_visc_index(self, index):
        return self.species[self.get_tag(index)].get_visc_index()
    
    
    def get_mpv(self, index=None):
        if index == None:
            return self.mpv
        return self.mpv
    
    
    def get_eu(self, index):
        return self.eu[index]
    
    
    def get_ev(self, index):
        return self.ev[index]
    
    
    def get_ew(self, index):
        return self.ew[index]
    
    
    def get_n_eff(self):
        return self.n_eff
    
    
    def get_particles_count(self):
        return self.num
    
    
    def get_x(self, index=None):
        if index == None:
            return self.x
        return self.x[index]
    
    
    def get_y(self, index):
        if index == None:
            return self.y
        return self.y[index]
    
    
    def get_velx(self, index=None):
        if index == None:
            return self.u
        else:
            return self.u[index]
    
    
    def get_vely(self, index=None):
        if index == None:
            return self.v
        else:
            return self.v[index]
    
    
    def get_velz(self, index=None):
        if index == None:
            return self.w
        else:
            return self.w[index]
    
    
    def set_x(self, x, index=None):
        if index == None:
            self.x = x
        else:
            self.x[index] = x
    
    
    def set_y(self, y, index=None):
        if index == None:
            self.y = y
        else:
            self.y[index] = y
    
    
    def set_velx(self, velx, index=None):
        if index == None:
            self.u = velx
        else:
            self.u[index] = velx
    
    
    def set_vely(self, vely, index=None):
        if index == None:
            self.v = vely
        else:
            self.v[index] = vely
    
    
    def set_velz(self, velz, index=None):
        if index == None:
            self.w = velz
        else:
            self.w[index] = velz