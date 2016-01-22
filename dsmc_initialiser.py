# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:16:38 2015

@author: mohit
"""
import numpy as np
import dsmc_particles as dm_p
import dsmc_cells as dm_c
import dsmc_sampler as dm_s



class Initialiser:
    @staticmethod
    def run(cells, gas, surf_group, volume, n_particles_in_cell, ref_point):
        cell_detector = CellDetector(cells, surf_group, ref_point)
        cells_in = cell_detector.detect_all()
        num = n_particles_in_cell * len(cells_in)
        particles = dm_p.Particles(num, gas.get_number_density() * volume / num)
        particles.setup(gas.get_mol_frac(), gas.get_species(), gas.get_mpv())
        particle_initialiser = ParticleInitialiser(particles, 
                                                   n_particles_in_cell, gas)
        particle_initialiser.run(cells, cells_in)
        distributor = dm_c.Distributor(cells, particles)
        distributor.distribute_all_particles()
        sampler = dm_s.Instant_sampler(cells, cells_in, particles, 
                                       gas.get_n_species())
        sampler.run()
        return [particles, cells_in]


# need to modify
class CellDetector:
    def __init__(self, cells, surf_group, ref_point):
        self.cells = cells
        self.surf_group = surf_group
    
    
    def detect_all(self):
        cells_in = []
        group_count = self.surf_group.get_group_count()
        for index in range(len(self.cells.get_temperature())):
            xc, yc = self.cells.get_center(index)
            if not (self._check_groups(xc, yc, 0, group_count - 1)):
                cells_in.append(index)
        return cells_in
    
    
    # This function checks if the cell center is inside a surface group
    # If a point is inside any groups of surfaces, then it is not in the domain
    # of the particles being considered.
    def _check_groups(self, xc, yc, start, end):
        if start == end:
            surf_count = self.surf_group.get_surf_count(start)
            return self._check_surf(xc, yc, start, 0, surf_count - 1)
        else:
            mid = int((start + end) / 2)
            return (self._check_groups(xc, yc, start, mid) or 
                    self._check_groups(xc, yc, mid + 1, end))
    
    
    # This function checks if the cell center is inside a specific group by 
    # checking over all the surfaces of the group
    # If a point is inside all the surfaces of the group then certainly it is
    # inside that object otherwise it is atleast not inside this object.
    def _check_surf(self, xc, yc, group_index, start, end):
        if start == end:
            return self._check_cell(xc, yc, group_index, start)
        else:
            mid = int((start + end) / 2)
            return (self._check_surf(xc, yc, group_index, start, mid) and
                    self._check_surf(xc, yc, group_index, mid + 1, end))
    
    
    # This function checks whether the cell center is in a direction opposite 
    # to the normal direction of the surface.
    def _check_cell(self, xc, yc, group_index, surf_index):
        normal = self.surf_group.get_surf_normal(group_index, surf_index)
        vertex1 = self.surf_group.get_surf_vertex1(group_index, surf_index)
        dot_prod = normal[0] * (xc - vertex1[0]) + normal[1] * (yc - vertex1[1])
        if dot_prod > 0:
            return False
        else:
            return True



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
        self._init_velocity()
    
    
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
    def _init_velocity(self):
        k = 1.3806488e-23
        c = np.sqrt(self.gas.get_gamma() * self.gas.get_temperature() * k / 
                    self.gas.get_mass())
        c1 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        c2 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        c3 = np.random.normal(0.0, 0.5, self.particles.get_particles_count())
        
        mpv = self.particles.get_mpv()        
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