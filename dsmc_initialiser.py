# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:16:38 2015

@author: mohit
"""
import numpy as np
import dsmc_particles as dm_p
import dsmc_detector as dm_d
import dsmc_cells as dm_c
import dsmc_sampler as dm_s



class Initialiser:
    @staticmethod    
    def run(cells, gas, domain, n_particles_in_cell, ref_point):
        cell_detector = CellDetector(cells, domain.get_surface(), ref_point)
        cells_in = cell_detector.detect_all()
        num = n_particles_in_cell * len(cells_in)
        particles = dm_p.Particles(num,
                    gas.get_number_density() * domain.get_volume() / num)
        particles.setup(gas.get_mol_frac(), gas.get_species(), gas.get_mpv())
        particle_initialiser = ParticleInitialiser(particles, 
                                                   n_particles_in_cell, gas)
        particle_initialiser.run(cells, cells_in)
        distributor = dm_c.Distributor(cells, particles)
        distributor.distribute_all_particles()
        sampler = dm_s.Instant_sampler(cells, cells_in, particles, 
                                       gas.get_n_species())
        sampler.run()
#        print cells.temperature
        return [particles, cells_in]



class CellDetector:    
    def __init__(self, cells, surface, ref_point):
        self.cells = cells
        self.surface = surface
        self.cell_center = []
        self.point_detector = dm_d.PointDetector(surface, ref_point)
    
    
    def detect_all(self):
        for index in range(len(self.cells.get_temperature())):
            self.cell_center.append(self.cells.get_center(index))
        return self.point_detector.detect_all(self.cell_center, 0)



class ParticleInitialiser:
    def __init__(self, particles, n_particles_in_cell, gas):
        self.particles = particles
        self.n_particles_in_cell = int(n_particles_in_cell)
        self.gas = gas
    
    
    def run(self, cells, cells_in):
        self.init_particles(cells, cells_in)
        self._init_cells(cells, cells_in)
    
    
    def init_particles(self, cells, cells_in):
        self._init_location(cells, cells_in)
        self._init_velocity(self.particles.get_mpv())
    
    
    def _init_location(self, cells, cells_in, offset=0):
        for index1, cell_index in enumerate(cells_in):
            xc, yc = cells.get_center(cell_index)
            index = index1 * self.n_particles_in_cell + offset
            for i in range(self.n_particles_in_cell):
                constt = ((2.0 * np.random.rand() - 1) / 2.0 * 
                            cells.get_cell_length(index + i) + xc)
                self.particles.set_x(constt, index + i)
                
                constt = ((2.0 * np.random.rand() - 1) / 2.0 *
                        cells.get_cell_width(index + i) + yc)
                self.particles.set_y(constt, index + i)
        return index + self.n_particles_in_cell
    
    
    # c is the speed of sound.
    def _init_velocity(self, mpv):
        k = 1.3806488e-23
        c = np.sqrt(self.gas.get_gamma() * self.gas.get_temperature() * k / 
                    self.gas.get_mass())
        c1 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        c2 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        c3 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        
        self.particles.set_velx(c1 * mpv + self.gas.get_mach_x() * c)
        self.particles.set_vely(c2 * mpv + self.gas.get_mach_y() * c)
        self.particles.set_velz(c3 * mpv + self.gas.get_mach_z() * c)
    
    
    def _init_cells(self, cells, cells_in):
        c = np.sqrt(self.gas.get_temperature() * self.gas.get_gamma() * 8.314)
        
        n_particles = (self.n_particles_in_cell * np.array(self.gas.get_mol_frac()))
        
        for index in cells_in:
            cells.set_temperature(self.gas.get_temperature(), index)
            
            for tag in range(len(self.gas.get_species())):
                cells.set_n_particles(n_particles[tag], tag, index)
            
            cells.set_velx(self.gas.get_mach_x() * c, index)
            cells.set_vely(self.gas.get_mach_y() * c, index)