# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 23:08:12 2015

@author: mohit
"""

"""
Distribution algorithm would change depending upon the geometry.
particles out is a domain geometry specific job.
So, it needs to be done in the distributor itself.
Here, I am using rectangular cells having constant length and width.
Also, it is assumed that the rectangular cells only have verticla and 
horizontal sides.
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
# It is assumed that the rect cells have only horizontal and vertical sides.
class RectCells:
    def __init__(self, n_x, n_y, length, width, centre, n_species=1):
        """ 
        storing the x values of cell centre in xc of cells at y = 0 and
        y values in yc at x = 0.
        """
        n_cells = int(n_x) * int(n_y)
        self.n_x = int(n_x)
        self.n_y = int(n_y)
        self.xc = np.zeros(n_x)
        self.yc = np.zeros(n_y)
        self.cell_length = float(length) / self.n_x
        self.cell_width = float(width) / self.n_y
        self.cell_volume = self.cell_length * self.cell_width
        self.particles_inside = [[] for i in range(n_cells)]
        self.u = np.zeros(n_cells)
        self.v = np.zeros(n_cells)
        self.w = np.zeros(n_cells)
        self.mass = np.zeros(n_cells)
        self.number_density = np.zeros((n_species, n_cells), dtype=float)
        # n_particles gives the no. of particles of each species at each cell
        # for e.g. no. of particles of species of tag 1 at cell 0 can be 
        # located by self.n_particles[0][0], (cell index is the second one)
        self.n_particles = np.zeros((n_species, n_cells), dtype = int)
        self.temperature = np.zeros(n_cells)
        self.gas_temperature = np.zeros((n_cells, n_species))
        self._find_cell_location(centre, float(length), float(width))
        self.n_species = n_species
        # datum = the lower left corner point of the domain.
        self.datum = (self.xc[0] - self.cell_length / 2.0,
                      self.yc[0] - self.cell_width / 2.0)
        self.x_max = self.datum[0] + float(length)
        self.y_max = self.datum[1] + float(width)
        self.domain_length = length
        self.domain_width = width
        self.particles_out = []
    

    # this function distributes particle in a rectangular cell.
    # assumption: particles are inside the domain.
    def distribute_all_particles(self, particles):
        self.reset_particles()
        self.particles_out = []
        for index in range(particles.get_particles_count()):            
            cell_index = self.find_cell_index(particles.get_x(index), 
                                              particles.get_y(index))
            if cell_index == None:
                self.particles_out.append(index)
            else:
                self.add_particle(cell_index, index, particles.get_tag(index))
    
    
    # randomising all particles.
    def shuffle_all_particles(self):
        for cell_particles in self.particles_inside:
            np.random.shuffle(cell_particles)
    
    
    def get_particles_out(self):
        return self.particles_out
    
    
    # assuming constant cell length and cell width.
    def find_cell_index(self, x, y):
        if (self.is_inside_domain(x, y)):
            x = x - self.datum[0]
            y = y - self.datum[1]
            cell_index = int(y / self.cell_width) * self.n_x
            cell_index += int(x / self.cell_length)
    #        print "cell_index = ", cell_index, x, y
            return cell_index
        else:
            return None
        
    
    # this function would reset the particles inside attribute of cells.
    def reset_particles(self):
        n_cells = self.n_x * self.n_y
        self.n_particles = np.zeros((self.n_species, n_cells), dtype = int)
        self.particles_inside = [[] for i in range(n_cells)]
    
    
    # this function would return a boundary cell given a boundary line segment.
    # Since creation of a boundary is a geometric feature, it is done here 
    # rather than inside the boundary module.
    def generate_cell(self, line):
        vertex1 = line[0]
        vertex2 = line[1]
        dx = vertex1[0] - vertex2[0]
        if dx * dx <= 1.0e-12:
            return self._generate_vert_cell(vertex1, vertex2)
        else:
            return self._generate_hor_cell(vertex1, vertex2)
    
    
    # this function is a helper function to generate_cell.
    # this function is used to generates boundary cells if the line segment is
    # horizontal.
    def _generate_hor_cell(self, vertex1, vertex2):
        x = (vertex1[0] + vertex2[0]) / 2.0
        dist = np.abs(vertex1[0] - vertex2[0])
        
        if vertex1[1] < self.y_max:
            y = self.datum[1] - self.cell_width / 2.0
        else:
            y = self.y_max + self.cell_width / 2.0
        return RectCells(int(round(dist / self.cell_length)), 1, dist, 
                         self.cell_width, (x, y), self.n_species)
    
    # this function is another helper function to generate_cell.
    # this function is used to generates boundary cells if the line segment is
    # vertical.
    def _generate_vert_cell(self, vertex1, vertex2):
        y = (vertex1[1] + vertex2[1]) / 2.0
        dist = np.abs(vertex1[1] - vertex2[1])
        
        if vertex1[0] < self.x_max:
            x = self.datum[0] - self.cell_length / 2.0
        else:
            x = self.x_max + self.cell_length / 2.0
        
        return RectCells(1, int(round(dist / self.cell_width)), self.cell_length,
                         dist, (x, y), self.n_species)
    
    
    # assuming all cells have same dimensions.
    # centre = domain centre, length = domain length and width = domain width.
    def _find_cell_location(self, centre, length, width):
        x_lower_left = centre[0] - length / 2.0 + self.cell_length / 2.0
        y_lower_left = centre[1] - width / 2.0 + self.cell_width / 2.0
        for i in range(self.n_x):
            self.xc[i] = x_lower_left +  self.cell_length * i
        
        for i in range(self.n_y):
            self.yc[i] = y_lower_left +  self.cell_width * i
    
    
    def is_inside_domain(self, x, y):
        if (x < self.datum[0] or y < self.datum[1] or
            x > self.x_max or y > self.y_max):
            return False
        else:
            return True
    
    
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
    
    
    # this function generates the location of cell centre given its cell index.
    def get_center(self, cell_index):
        xc = self.xc[cell_index % self.n_x]
        yc = self.yc[int(cell_index / self.n_x)]
        return (xc, yc)
    
    # this function sets the value of the xc value of a certain cell index to a 
    # new value.
    def set_cell_x_center(self, cell_index, xc):
        cell_index = cell_index / self.n_x
        self.xc[cell_index] = xc
    
    
    # this function sets the value of the yc value of a certain cell index to a
    # new value.
    def set_cell_y_center(self, cell_index, yc):
        cell_index = cell_index % self.n_y
        self.yc[cell_index] = yc
    
    
    def get_cell_length(self, cell_index=None):
        return self.cell_length
    
    
    def get_cell_width(self, cell_index=None):
        return self.cell_width
    
    
    def get_cell_volume(self):
        return self.cell_volume
    
    
    def get_mass(self, cell_index=None):
        if cell_index == None:
            return self.mass
        else:
            return self.mass[cell_index]
    
    
    def set_mass(self, mass, cell_index=None):
        if cell_index == None:
            self.mass = mass
        else:
            self.mass[cell_index] = mass
    
    
    def get_temperature(self, cell_index=None):
        if cell_index == None:
            return self.temperature
        else:
            return self.temperature[cell_index]
    
    
    def set_temperature(self, temperature, cell_index=None):
        if cell_index == None:
            self.temperature = temperature
        else:
            self.temperature[cell_index] = temperature
    
    
    def get_number_density(self, tag=None, cell_index=None):
        if tag == None and cell_index == None:
            return self.number_density
        else:
            return self.number_density[tag][cell_index]
    
    
    def set_number_density(self, number_density, tag=None, cell_index=None):
        if tag == None and cell_index == None:
            self.number_density = number_density
        else:
            self.number_density[tag][cell_index] = number_density
    
    
    def get_n_particles(self, tag=None, cell_index=None):
        if tag == None and cell_index == None:
            return self.n_particles
        else:
            return self.n_particles[tag][cell_index]
    
    
    def set_n_particles(self, n_particles, tag=None, cell_index=None):
        if tag == None and cell_index == None:
            self.n_particles = n_particles
        else:
            self.n_particles[tag][cell_index] = n_particles
    
    
    def get_particles_inside(self, cell_index=None):
        if cell_index == None:
            return self.particles_inside
        else:
            return self.particles_inside[cell_index]
    
    
    def add_particle(self, cell_index, particle_index, tag):
        self.particles_inside[cell_index].append(particle_index)
        self.n_particles[tag][cell_index] += 1
    
    
    def get_velx(self, cell_index=None):
        if cell_index == None:
            return self.u
        else:
            return self.u[cell_index]
    
    
    def get_vely(self, cell_index=None):
        if cell_index == None:
            return self.v
        else:
            return self.v[cell_index]
    
    
    def get_velz(self, cell_index=None):
        if cell_index == None:
            return self.w
        else:
            return self.w[cell_index]
    
    
    def set_velx(self, velx, cell_index=None):
        if cell_index == None:
            self.u = velx
        else:
            self.u[cell_index] = velx
    
    
    def set_vely(self, vely, cell_index=None):
        if cell_index == None:
            self.v = vely
        else:
            self.v[cell_index] = vely
    
    
    def set_velz(self, velz, cell_index=None):
        if cell_index == None:
            self.w = velz
        else:
            self.w[cell_index] = velz