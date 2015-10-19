# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:08:12 2015

@author: mohit
"""

import numpy as np

# Cells would only the instantaneous values. For sampling, these instantaneous
# values should be used. In other words, sampler should store its own set of 
# values for sampling.
# assuming cells are of equal length and width.
# only horizontal or vertical orientation may be formed.
# assuming these cells would only be used in a rectangular domain.
# vertices contain the vertices of domain in anti-clockwise direction
# starting at leftmost.
# particle_inside would contain particles in a cell. each element of 
# particle inside is a list of particles inside. For, e.g., 1st element
# would contain list particles inside 1st cell and so on.
class RectCells:
    def __init__(self, n_x, n_y, length, width, centre, gas):
        """ 
        just storing the x values in xc at y = 0 and y values in yc 
        at x = 0.
        """
        n_cells = int(n_x) * int(n_y)
        self.n_x = int(n_x)
        self.n_y = int(n_y)
        self.xc = np.zeros(n_x)
        self.yc = np.zeros(n_y)
        self.length = np.ones(n_cells) * length / self.n_x
        self.width = np.ones(n_cells) * width / self.n_y
        self.volume = np.ones(n_cells) * self.length * self.width
        self.particles_inside = [[] for i in range(n_cells)]
        self.pressure = np.zeros(n_cells)
        self.u = np.zeros(n_cells)
        self.v = np.zeros(n_cells)
        self.w = np.zeros(n_cells)
        self.mass = np.zeros(n_cells)
        self.number_density = np.zeros(n_cells)
        # n_species gives the no. of particles of different species
        self.n_species = np.zeros((n_cells, len(gas.species)), dtype = int)
        self.mach = np.zeros(n_cells)
        self.density = np.zeros(n_cells)
        self.temperature = np.zeros(n_cells)
        self.gas_temperature = np.zeros((n_cells, len(gas.species)))
        self._find_cell_location(centre, length, width)
    
    
    # this function returns the array of locations of the cells.
    def get_all_center(self):
        stack = []
        self.get_subset_center(stack, 0, self.n_x * self.n_y - 1)
        return (np.asarray(stack))
        
    
    # stack should be a list. start is the starting index of cell and end is 
    # the last index of cell.
    def get_subset_center(self, stack, start, end):
        if (start == end):
            stack.append(self.get_center(start))
        else:
            mid = int((start + end) / 2)
            self.get_subset_center(stack, start, mid)
            self.get_subset_center(stack, mid + 1, end)
    
    
    # this function generates the location of cells based on xc and yc values 
    # and its cell number.
    def get_center(self, cell_index):
        xc = self.xc[cell_index % self.n_x]
        yc = self.yc[int(cell_index / self.n_x)]
        return (xc, yc)
    
    
    # assuming equal length and width.
    def find_cell_index(self, x, y):
        datum = (self.xc[0]-self.length[0]/2.0, self.yc[0]-self.width[0]/2.0)
        x = x - datum[0]
        y = y - datum[1]
        cell_index = int(y / self.width[0]) * self.n_x
        cell_index += int(x / self.length[0])
#        print "cell_index = ", cell_index, x, y
        return cell_index
        
    
    # assuming all cells have same dimensions.
    def _find_cell_location(self, centre, length, width):
        x_lower_left = centre[0] - length / 2.0 + self.length[0] / 2.0
        y_lower_left = centre[1] - width / 2.0 + self.width[0] / 2.0
        for i in range(self.n_x):
            self.xc[i] = x_lower_left +  self.length[0] * i
        
        for i in range(self.n_y):
            self.yc[i] = y_lower_left +  self.width[0] * i



class Distributor:
    def __init__(self, cells, particles):
        self.cells = cells
        self.particles = particles
    
    
    def find_tag_particles(self, tag):
        tag_particles = []
        for index in self.cells.particles_inside:
            if self.particles.tag[index] == tag:
                tag_particles.append
        
        return (np.asarray(tag_particles))
    
    
    # this function distributes particle in a rectangular cell.
    # assumption: particles are inside the domain
    def distribute_all_particles(self):
        self.reset_particles_inside()
#        print "n_cells = ", len(self.cells.xc) * len(self.cells.yc)
        for index in range(len(self.particles.x)):
            y = self.particles.x[index] - 0.2
            if (y > self.particles.y[index]):
                vel = (self.particles.u[index], self.particles.v[index])
                loc = (self.particles.x[index] - vel[0] * 1e-5,
                       self.particles.y[index] - vel[1] * 1e-5)
                if (loc[1] + 0.2 < loc[0] or loc[0] < 0.0 or loc[1] < 0.0 or 
                    loc[0] > 1.0 or loc[1] > 1.0):
                    flag = True
                else:
                    flag = False
                print "not reflected index loc vel flag  = ", index, loc, vel, flag
            cell_index = self.cells.find_cell_index(self.particles.x[index], 
                                                    self.particles.y[index])
#           print "cell_index = ", cell_index, index
#           print "particle_x = ", self.particles.x[index]
#           print "particle_y = ", self.particles.y[index]
            self.cells.particles_inside[cell_index].append(index)
            tag = self.particles.tag[index]
            self.cells.n_species[cell_index][tag] += 1
        
    
    # after every time step this function needs to be called.
    def reset_particles_inside(self):
        n_cells = self.cells.n_x * self.cells.n_y
        n_species = len(self.cells.gas_temperature[0])
        self.cells.n_species = np.zeros((n_cells, n_species), dtype = int)
        for i in range (n_cells):
            self.cells.particles_inside[i] = []
            