# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:16:38 2015

@author: mohit
"""
import numpy as np
import particles as dm_p
import sampler as dm_s


class Initialiser:
    @staticmethod
    def run(cells, gas, surf_group, domain, n_par_in_cell,
            closed=False, multiphase=False):
        if closed:
            cells_in = range(len(cells.get_temperature()))
        else:
            cell_detector = CellDetector(cells, surf_group)
            cells_in = cell_detector.detect_all()
        f = open('cells_in.txt', 'w')
        np.savetxt(f, cells_in)
        f.close()
        exit()
        num = n_par_in_cell * len(cells_in)
        n_eff = gas.get_number_density() * domain.get_volume() / num
        
        particles = dm_p.Particles(num, n_eff)
        particles.setup(gas.get_mol_frac(), gas.get_species(), gas.get_mpv())
        
        particle_init = ParticleInitialiser(particles, n_par_in_cell, gas, domain)
        particle_init.run(cells, cells_in, multiphase)
        
        cells.distribute_all_particles(particles)
        sampler = dm_s.Instant_sampler(cells, cells_in, particles, 
                                       gas.get_n_species())
        sampler.run()
        print "cell temperature = ", cells.temperature.mean()
        return [particles, cells_in]


# need to modify
class CellDetector:
    def __init__(self, cells, surf_group):
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
    def __init__(self, particles, n_par_in_cell, gas, domain):
        self.particles = particles
        self.n_par_in_cell = int(n_par_in_cell)
        self.gas = gas
        self.domain = domain
    
    
    # It is assumed that only below mentioned constraints is applied during 
    # initialization. The following are the constrainsts :-
    # 1. There are different patches where different molecules may be 
    #    initialized. for e.g. shocktube
    # 2. There is a constraint of space for every molecule in the domain.
    #    for e.g. flow past a wedge.
    #
    # The above constraints can be satisfied simultaneously. The following can 
    # is the solution in this implementation.
    # Find the cell indices that satisfy the constraints.
    # For cell index, check all the patches it lies in and consequently make a
    # sequence of these cell indices for each patch.
    # Initialise the particles representing the molecule inside that patch in 
    # their respective cells.
    def run(self, cells, cells_in, multiphase):
        if self.domain.is_patched():
            if multiphase:
                
                for tag in self.domain.get_tags():
                    tag_cells_in = self._find_tag_cells(cells, cells_in, tag)
                    tag_num = self.particles.get_tag_num(tag)
                    tag_start = self.particles.get_tag_start(tag)
                    self._init_multi_loc(cells, tag_cells_in[tag],
                                         tag_num, tag_start)
                
            else:
                self._init_patched_loc()
            
        else:
            self._init_location(cells, cells_in)
        
        self._init_velocity()
    
    
    def _find_tag_cells(self, cells, cells_in, tag):
        tag_cells_in = []
        x_min, x_max = self.domain.get_patch_x(tag)
        if x_min > x_max:
            c = x_min
            x_min = x_max
            x_max = c    
        
        y_min, y_max = self.domain.get_patch_x(tag)
        if y_min > y_max:
            c = y_min
            y_min = y_max
            y_min = c
        
        for cell_index in cells_in:
            xc, yc = cells.get_center(cell_index)
            length = cells.get_cell_length(cell_index)
            width = cells.get_cell_width(cell_index)
            
            cell_x_min = xc - length * 0.5
            cell_x_max = xc + length * 0.5
            cell_y_min = yc - width * 0.5
            cell_y_max = yc + width * 0.5
            
            if (cell_x_min > x_min and cell_x_max < x_max and 
                cell_y_min > y_min and cell_y_max < y_max):
                    tag_cells_in.append(cell_index)
        
        return tag_cells_in
    
    
    def _init_multi_loc(self, cells, cells_in, tag_num, tag_start):
        sample = np.random.choice(cells_in, tag_num)
        for index, cell_index in enumerate(sample):
            xc, yc = cells.get_center(cell_index)
            particle_index = index + tag_start
            length = cells.get_cell_length(cell_index)
            width = cells.get_cell_width(cell_index)
            
            x = ((2.0 * np.random.rand() - 1) / 2.0 * length + xc)
            y = ((2.0 * np.random.rand() - 1) / 2.0 * width + yc)
            
            self.particles.set_x(x, particle_index)
            self.particles.set_y(y, particle_index)
    
    
    def _init_patched_loc(self):
        x = np.zeros(len(self.particles.get_x()))
        y = np.zeros(len(self.particles.get_y()))
        
        start = 0
        for tag in range(self.gas.get_n_species()):
            end = self.particles.get_tag_num(tag) + start
            
            x_min, x_max = self.domain.get_patch_x(tag)
            if x_min > x_max:
                c = x_min
                x_min = x_max
                x_max = c    
                
            y_min, y_max = self.domain.get_patch_x(tag)
            if y_min > y_max:
                c = y_min
                y_min = y_max
                y_min = c
                
            x[start:end] = np.random.random(end - start) * (x_max - x_min) + x_min
            y[start:end] = np.random.random(end - start) * (y_max - y_min) + y_min
            start = end
        
        self.particles.set_x(x)
        self.particles.set_y(y)
    
    
    def _init_location(self, cells, cells_in):
        sample = np.random.permutation(len(self.particles.get_x()))
        for index1, cell_index in enumerate(cells_in):
            xc, yc = cells.get_center(cell_index)
            index = index1 * self.n_par_in_cell
            length = cells.get_cell_length(cell_index)
            width = cells.get_cell_width(cell_index)
            for i in range(self.n_par_in_cell):
                x = ((2.0 * np.random.rand() - 1) / 2.0 * length + xc)
                self.particles.set_x(x, sample[index + i])
            
                y = ((2.0 * np.random.rand() - 1) / 2.0 * width + yc)
                self.particles.set_y(y, sample[index + i])
    
    
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
