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
    def __init__(self, inlet=None, zero_grad=None, outlet=None):
        self.inlet = inlet
        self.zero_grad = zero_grad
        self.outlet = outlet
        self.surface = []
        self.volume = 0.0
        self.surface_temp = 0.0
        if (inlet == None or len(inlet) == 0) and (
            outlet == None or len(outlet == 0)):
            self.bool = False
        else:
            self.bool = True
    
    def get_inlet(self):
        return self.inlet
    
    def get_outlet(self):
        return self.outlet
    
    def get_zero_grad(self):
        return self.zero_grad
    
    def get_surface(self):
        return self.surface
    
    def get_volume(self):
        return self.volume
    
    def get_surface_temp(self):
        return self.surface_temp
    
    def get_count_zero_grad(self):
        return len(self.zero_grad)
    
    def set_outlet(self, outlet):
        self.outlet = outlet
    
    def set_surface(self, surface, surface_temp=0.0):
        self.surface = surface
        self.surface_temp = surface_temp
    
    def set_volume(self, volume):
        self.volume = volume
    
    def set_inlet(self, inlet):
        self.inlet = inlet
    
    def set_zero_grad(self, zero_grad):
        self.zero_grad = zero_grad
    
    def is_open(self):
        return self.bool
    
    def is_inlet_on(self):
        return self.inlet != None
    
    def is_outlet_on(self):
        return self.outlet != None
    
    def is_zero_grad_on(self):
        return self.zero_grad != None



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


# tag represents the species type of particle.
# particles are stored in an ascending order of the species.
# for eg. species 1 would be stored in the first n1 particles the species2 for
# the next n2 particles and so on.
class Particles:
    def __init__(self, n_particles):
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
        self.visc_index = np.ones(n_particles)
        self.ref_temp = np.ones(n_particles)
        self.viscosity_coeff = np.ones(n_particles)
        self.dof = np.ones(n_particles)
        self.gamma = np.ones(n_particles)
        self.mpv = np.ones(n_particles)
        self.num = n_particles
        self.n_eff = 0.0
        self.tag = np.zeros(n_particles, dtype = int)


    # this function should only be called once.
    # sets the particles according to the molecules it represents.
    def setup(self, domain_volume, mole_fraction, number_density, species, mpv):
        n_species = len(mole_fraction)
        self.n_eff = number_density * domain_volume / self.num
        for index in range(self.num):
            tag = index % n_species
            self.mass[index] = species[tag].get_mass()
            self.dia[index] = species[tag].get_dia()
            self.cross_area[index] = 0.25 * np.pi * species[tag].get_dia() ** 2.0
            self.visc_index[index] = species[tag].get_visc_index()
            self.ref_temp[index] = species[tag].get_ref_temp()
            self.viscosity_coeff[index] = species[tag].get_visc_coeff()
            self.dof[index] = species[tag].get_dof()
            self.gamma[index] = species[tag].get_gamma()
            self.mpv[index] = mpv[tag]
            self.tag[index] = tag


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
        self.eu = self.u * self.u * self.mass
        self.ev = self.v * self.v * self.mass
        self.ew = self.w * self.w * self.mass
    
    
    def get_tag(self, index):
        return self.tag[index]
    
    
    def get_dia(self, index=None):
        if index == None:
            return self.dia
        else:
            return self.dia[index]
    
    
    def get_mass(self, index):
        return self.mass[index]
    
    
    def get_eu(self, index):
        if index == None:
            return self.eu
        else:
            return self.eu[index]
    
    
    def get_ev(self, index):
        if index == None:
            return self.ev
        else:
            return self.ev[index]
    
    
    def get_ew(self, index):
        if index == None:
            return self.ew
        else:
            return self.ew[index]
    
    
    def get_n_eff(self):
        return self.n_eff
    
    
    def get_ref_temp(self, index):
        return self.ref_temp[index]
    
    
    def get_visc_index(self, index):
        return self.visc_index[index]
    
    
    def get_mpv(self, index=None):
        if index == None:
            return self.mpv
        else:
            return self.mpv[index]
    
    
    def get_particles_count(self):
        return self.num
    
    
    def get_x(self, index=None):
        if index == None:
            return self.x
        else:
            return self.x[index]
    
    
    def set_x(self, x, index=None):
        if index == None:
            self.x = x
        else:
            self.x[index] = x
    
    
    def get_y(self, index=None):
        if index == None:
            return self.y
        else:
            return self.y[index]
    
    
    def set_y(self, y, index=None):
        if index == None:
            self.y = y
        else:
            self.y[index] = y
    
    
    def get_velx(self, index=None):
        if index == None:
            return self.u
        else:
            return self.u[index]
    
    
    def set_velx(self, velx, index=None):
        if index == None:
            self.u = velx
        else:
            self.u[index] = velx
    
    
    def get_vely(self, index=None):
        if index == None:
            return self.v
        else:
            return self.v[index]
    
    
    def set_vely(self, vely, index=None):
        if index == None:
            self.v = vely
        else:
            self.v[index] = vely
    
    
    def get_velz(self, index=None):
        if index == None:
            return self.w
        else:
            return self.w[index]
    
    
    def set_velz(self, velz, index=None):
        if index == None:
            self.w = velz
        else:
            self.w[index] = velz



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
        constt = ((species[-1].get_tag() + 2) * (species[-1].get_tag() + 1)) / 2
        self.mean_f_path = np.zeros(species[-1].get_tag() + 1)
        self.mean_col_time = np.zeros(species[-1].get_tag() + 1)
        self.mean_col_rate = np.zeros(constt)
        self.col_density_rate = np.zeros(constt)
        self.dia = np.zeros(constt)
        self.ref_temp = np.zeros(constt)
        self.mpv = np.zeros(len(species))
        self.reduced_mass = 0.0
        self.visc_index = np.zeros(constt)
        self.gamma = 0.0
        self.mass = 0.0
    
    
    def setup(self):
        self._find_number_density()
        self._find_gas_property()
        self._find_col_property()
    
    
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
    
    
    def _find_number_density(self):
        for index, species in enumerate(self.species):
            self.number_density += species.get_number_density()
    
    
    def _find_gas_property(self):
        self.gamma = 0.0
        self.mass = 0.0
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
#        print "reduced = ", self.reduced_mass
    
    
    # bug here.
    def _find_col_property(self):
        # boltzmann constant
        k = 1.3806488e-23
        for index1, mol1 in enumerate(self.species):
            for index2 in range(index1, self.species[-1].get_tag() + 1):
                mol2 = self.species[index2]
                index = ((index1 * (index1 + 1)) / 2)
                index += index2
                self.ref_temp[index] = self.species[index2].get_ref_temp()
                self.dia[index] = 0.5 * (mol1.get_dia() + mol2.get_dia())
                self.visc_index[index] = 0.5 * (mol1.get_visc_index() 
                                            + mol2.get_visc_index())
                self.ref_temp[index] = 0.5 * (mol1.get_ref_temp() 
                                            + mol2.get_ref_temp())
            constt = 0.0
            for index2, mol2 in enumerate(self.species):
                if index1 < index2:
                    index = ((index1 * (index1 + 1)) / 2)
                    index += index2
                else:
                    index = ((index2 * (index2 + 1)) / 2)
                    index += index1
                self.mean_f_path[index1] +=((1+mol1.get_mass()/mol2.get_mass())
                        ** 0.5 * (self.ref_temp[index] / self.temperature)
                        ** (self.visc_index[index] - 0.5) * self.dia[index]
                        ** 2.0 * self.mol_frac[index2])
                
                constt = (self.dia[index] ** 2.0 * self.mol_frac[index2]
                        * (self.temperature / self.ref_temp[index]) ** 
                        (1.0 - self.visc_index[index]) * 
                        self.ref_temp[index] ** 0.5)
                
                self.mean_col_rate[index1] += constt
                
                if index1 != index2:
                    self.col_density_rate[index] = self.mol_frac[index1] * constt
                else:
                    self.col_density_rate[index]=self.mol_frac[index1]*constt*0.5
        
            self.mean_f_path[index1] = 1 / self.mean_f_path[index1]
#            print self.mean_f_path[index1]
            self.mean_f_path[index1] /= (np.pi * self.number_density)
            self.mean_col_rate[index1] *= 2.0 * self.number_density * np.sqrt(2.0 * 
                                        k * np.pi / self.reduced_mass)
            self.mean_col_time[index1] = 1 / self.mean_col_rate[index1]
#            print "col_rate = ", self.mean_col_rate[index1]