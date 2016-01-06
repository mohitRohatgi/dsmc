# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:05:46 2015

@author: mohit
"""
import dsmc_cells as dm_c
import dsmc_particles as dm_p
import dsmc_detector as dm_d
import numpy as np


# assuming only open boundary conditions are applied here. Otherwise, the 
# closed boundary should be defined in surface of the domain and reflector 
# should take care of it. This class only deals with the inlet and outlet 
#  parts of domain.
# only the boundary particles entering from inlet are considered to affect the
# flow. Particles going outside through the outlet are neglected.
# this boundary condition can only be applied in aerospace or high speed flows.
# for dense flows surface generators should be used.
# this technique used for the implementation of open boundary in the simulation
# is  known as volume generator.
# It should have an adder and a deletor for adding net extra particles or 
# deleting net particle that has crossed the domain from particle array.
class Boundary:
    def __init__(self, cells, gas, domain, n_particles_in_cell, particles, ref_point):
        self.cells = cells
        self.gas = gas
        self.domain = domain
        self.n_particles_in_cell = n_particles_in_cell
        self.particles = particles
        self.ref_point = ref_point
        self.cell_gen = CellGenerator(cells, domain, ref_point)
        self.b_cells = self.cell_gen.run()
        self.slope = self.cell_gen.slope
        self.cell_detector = CellDetector(cells, self.b_cells, self.slope)
        self.adj_cells = self.cell_detector.run()
        self.particle_setter = ParticleSetter(cells, self.b_cells, 
                                              self.adj_cells, gas, particles)
        self.total_particles = np.zeros(len(self.b_cells), dtype = int)
        self.particle_detector = ParticleDetector(particles, domain, ref_point,
                                                  self.total_particles)
        self.particle_modifier = ParticleModifier(particles)
    
    
    # this function would either add or delete particles from the particles 
    # depending upon the no. of particles crossing into the domain from 
    # boundary and no. of particles crossing away from the domain through the
    # outlet.
    # the particles crossing through outlet are re-initialised with values
    # corresponding to the particles crossing into domain from boundary cells.
    def run(self, dt):
        self._setup()
        b_particles = dm_p.Particles(sum(self.total_particles))
        self.particle_setter.run(b_particles)
        b_particles.move_all(dt)
        particles_out, b_particles_in = self.particle_detector.run(b_particles)
        self.particle_modifier.run(particles_out, b_particles_in, b_particles)
    
    
    # this function would initialise particle according to the number_density of
    # adjacent cells and then find particles crossing into the domain.
    def _setup(self):
        for b_index in range(len(self.b_cells)):
            b_cells = self.b_cells[b_index]
            self.total_particles[b_index] = 0
            for cell_index in range(len(b_cells.xc)):
                adj_cell = self.adj_cells[b_index][cell_index]
                self.total_particles[b_index] += int(sum(self.cells.n_species[adj_cell]))
                b_cells.number_density[cell_index] = self.cells.number_density[adj_cell]
                b_cells.temperature[cell_index] = self.cells.temperature[adj_cell]
                b_cells.length[cell_index] = self.cells.length[adj_cell]
                b_cells.width[cell_index] = self.cells.width[adj_cell]
                b_cells.volume[cell_index] = self.cells.volume[adj_cell]
                b_cells.u[cell_index] = self.cells.u[adj_cell]
                b_cells.v[cell_index] = self.cells.v[adj_cell]



class CellGenerator:
    def __init__(self, cells, domain, ref_point):
        self.domain = domain
        self.cells = cells
        self.ref_point = ref_point
        self.slope = []
        self.b_cells = []
    
    
    def run(self):
        self._generate(0, len(self.domain.inlet) - 2)
        return self.b_cells
    
    
    def _generate(self, start, end):
        start = int(start)
        end = int(end)
        if start == end:
            self.b_cells.append(self._generate_segment(start))
        else:
            mid = int((start + end) / 2)
            self._generate(start, mid)
            self._generate(mid + 1, end)
    
    
    # assuming domain is rectangular.
    # assuming cells have constant length and width.
    # slope represents the slope of the inlet boundary.
    def _generate_segment(self, inlet_index):
        vertex1 = self.domain.inlet[inlet_index]
        vertex2 = self.domain.inlet[inlet_index + 1]
        if (vertex1[0] - vertex2[0]) != 0:
            slope = (vertex1[1] - vertex2[1]) / (vertex1[0] - vertex2[0])
        else:
            slope = 1.0e8
        self.slope.append(slope)
        if slope >= 1.0e8:
            sin = 1.0
            cos = 0.0
        else:
            if slope == 0.0:
                sin = 0.0
                cos = 1.0
            else:
                sin = slope / np.sqrt(slope ** 2.0 + 1)
                cos = sin / slope
        length = (vertex1[0] - vertex2[0]) ** 2
        length += (vertex1[1] - vertex2[1]) ** 2
        length = np.sqrt(length)
        n_x = round(length / self.cells.length[0])
        centre = [(vertex1[0] + vertex2[0]) / 2.0, (vertex1[1] + vertex2[1]) / 2.0]
        datum = self._find_datum(centre, sin, cos, self.cells.width[0])
        b_cell = dm_c.RectCells(n_x, 1, length, self.cells.width[0], datum)
        b_cell = self._transform_cell(b_cell, sin, cos, datum)
        return b_cell
    
    
    def _find_datum(self, centre, sin, cos, length):
        sine = -cos
        cos = sin
        l_x = length / 2.0 * cos
        l_y = length / 2.0 * sine
        origin1 = (centre[0] + l_x, centre[1] + l_y)
        origin2 = (centre[0] - l_x, centre[1] - l_y)
        dist1 = (origin1[0] - self.ref_point[0]) ** 2.0
        dist1 += (origin1[1] - self.ref_point[1]) ** 2.0
        dist2 = (origin2[0] - self.ref_point[0]) ** 2.0
        dist2 += (origin2[1] - self.ref_point[1]) ** 2.0
        if dist1 > dist2:
            return origin1
        else:
            return origin2
    
    
    def _transform_cell(self, b_cell, sin, cos, datum):
        y = b_cell.yc[0]
        b_cell.yc = np.zeros(b_cell.n_x)
        for i in range(b_cell.n_x):
            b_cell.yc[i] = y
        
        for i in range (len(b_cell.xc)):
            a = b_cell.xc[i] - datum[0]
            b = b_cell.yc[i] - datum[1]
            b_cell.xc[i] = a * cos - b * sin + datum[0]
            b_cell.yc[i] = a * sin + b * cos + datum[1]

        return b_cell



class CellDetector:
    def __init__(self, cells, b_cells, slope):
        self.cells = cells
        self.ref_point = (self.cells.xc[int(self.cells.n_x / 2)],
                          self.cells.yc[int(self.cells.n_y / 2)])
        self.adj_cells = []
        self.b_cells = b_cells
        self.slope = slope
        
    
    def run(self):
        self._detect(0, len(self.b_cells) - 1)
        return self.adj_cells
    
    
    def _detect(self, start, end):
        start = int(start)
        end = int(end)
        if start == end:
            self.adj_cells.append(self._detect_all_adj_cells(start, 0, 
                                        len(self.b_cells[start].xc) - 1, []))
        else:
            mid = int((start + end) / 2)
            self._detect(start, mid)
            self._detect(mid + 1, end)
    
    
    def _detect_all_adj_cells(self, index, start, end, stack):
        start = int(start)
        end = int(end)
        if start == end:
            stack.append(self._detect_adj_cell(index, start))
        else:
            mid = int((start + end) / 2)
            self._detect_all_adj_cells(index, start, mid, stack)
            self._detect_all_adj_cells(index, mid + 1, end, stack)
            return stack
    
    
    def _detect_adj_cell(self, b_index, cell_index):
        x = self.b_cells[b_index].xc[cell_index]
        y = self.b_cells[b_index].yc[cell_index]
        slope = self.slope[b_index]
        if slope == 0.0:
            sin = 1.0
            cos = 0.0
        elif slope >= 1.0e8:
            sin = 0.0
            cos = 1.0
        else:
            sin = - 1.0 / np.sqrt(slope ** 2.0 + 1.0)
            cos = sin * (-slope)
        dist = self.b_cells[b_index].width[cell_index]
        x1 = x + dist * cos
        y1 = y + dist * sin
        x2 = x - dist * cos
        y2 = y - dist * sin
        dist1 = (x1 - self.ref_point[0])** 2.0 + (y1 - self.ref_point[1]) ** 2.0
        dist2 = (x2 - self.ref_point[0])** 2.0 + (y2 - self.ref_point[1]) ** 2.0
        if dist1 > dist2:
            cell_index = self.cells.find_cell_index(x2, y2)
        else:
            cell_index = self.cells.find_cell_index(x1, y1)
        return cell_index



class ParticleDetector:
    def __init__(self, particles, domain, ref_point, total_particles):
        self.particles = particles
        self.domain = domain
        self.ref_point = ref_point
        self.total_particles = total_particles
        self.particle_detector = dm_d.PointDetector(domain.outlet, ref_point)
        self.b_particle_detector = dm_d.PointDetector(domain.inlet, ref_point)
    
    
    def run(self, b_particles):
        particles_out = self._detect_particles_out()
#        print "particles_out = ", particles_out
        b_particles_in = self._detect_b_particles_in(b_particles)
        return (particles_out, b_particles_in)
    
    
    # this function returns the particles out of domain that were inside the 
    # particles array in the domain.
    def _detect_particles_out(self):
        points = []
        for index in range(len(self.particles.x)):
            point = (self.particles.x[index], self.particles.y[index])
            points.append(point)
#        print "points = ", points
        return self.particle_detector.detect_all(points, 1)
    
    
    def _detect_b_particles_in(self, b_particles):
        points = []
        for index in range(len(b_particles.x)):
            point = (b_particles.x[index], b_particles.y[index])
            points.append(point)        
        b_particles_in = self.b_particle_detector.detect_all(points)
        return b_particles_in



# this class would take the particle array and initialise it according to its
# adj cells value. 
class ParticleSetter:
    def __init__(self, cells, b_cells, adj_cells, gas, particles):
        self.particles = particles
        self.cells = cells
        self.gas = gas
        self.b_cells = b_cells
        self.adj_cells = adj_cells
        self.count = 0
    
    
    def run(self, b_particles):
        self.count = 0
        b_particles = self._set_all(0, len(self.b_cells) - 1, b_particles)
    
    
    # this function generates particles in all of the boundaries.
    def _set_all(self, start, end, b_particles):
        start = int(start)
        end = int(end)
        if start == end:
            self._set_segment(start, 0, len(self.b_cells[start].xc) - 1, b_particles)
        else:
            mid = int(start + end) / 2
            self._set_all(start, mid, b_particles)
            self._set_all(mid + 1, end, b_particles)
    
    
    # this function generates particles in all of the cells of a boundary.
    def _set_segment(self, b_index, start, end, b_particles):
        start = int(start)
        end = int(end)
        if start == end:
            self._set_in_cell(b_index, start, b_particles)
        else:
            mid = int((start +end) / 2)
            self._set_segment(b_index, start, mid, b_particles)
            self._set_segment(b_index, mid + 1, end, b_particles)
    
    # this function generates particles in a particular cell of certain boundary.
    def _set_in_cell(self, b_index, cell_index, b_particles):
        adj_cell = self.adj_cells[b_index][cell_index]
        n_species = self.cells.n_species[adj_cell]
        for tag in range(len(self.gas.species)):
            self._set_particles(b_particles, n_species[tag], tag, b_index, 
                                cell_index)
    
    
    # this function only sets the attributes of the particles corresponding to 
    # its species and the cell.
    def _set_particles(self, b_particles, n_particle, tag, b_index, cell_index):
        end = int(self.count + n_particle)
        molecule = self.gas.species[tag]
        k = 1.3806488e-23
        T = self.b_cells[b_index].temperature[cell_index]
        mpv = np.sqrt(2.0 * k * T / molecule.mass)
        b_cell = self.b_cells[b_index]
        length = b_cell.length[cell_index]
        width = b_cell.width[cell_index]
        datum_x = b_cell.xc[cell_index] - length / 2.0
        datum_y = b_cell.yc[cell_index] - width / 2.0
        b_particles.x[self.count:end] = datum_x + np.random.rand(n_particle) * length
        b_particles.y[self.count:end] = datum_y + np.random.rand(n_particle) * width
        b_particles.tag[self.count:end] = tag
        b_particles.mass[self.count:end] = molecule.mass
        b_particles.dia[self.count:end] = molecule.dia
        b_particles.viscosity_index[self.count:end] = molecule.viscosity_index
        b_particles.ref_temperature[self.count:end] = molecule.ref_temperature
        b_particles.viscosity_coeff[self.count:end] = molecule.viscosity_coeff
        b_particles.dof[self.count:end] = molecule.dof
        b_particles.gamma[self.count:end] = molecule.gamma
        b_particles.mpv[self.count:end] = mpv
        
        mpv = mpv / np.sqrt(3.0)
        
        n = int(end - self.count)
        
        c1 = np.random.normal(0.0, 2.0, n)
        c2 = np.random.normal(0.0, 2.0, n)
        c3 = np.random.normal(0.0, 2.0, n)
        
        
        b_particles.u[self.count:end] = (2.0 * c1 - np.ones(n)) * mpv
        b_particles.v[self.count:end] = (2.0 * c2 - np.ones(n)) * mpv
        b_particles.w[self.count:end] = (2.0 * c3 - np.ones(n)) * mpv
        
        b_particles.u[self.count:end] += np.ones(n) * b_cell.u[cell_index]
        b_particles.v[self.count:end] += np.ones(n) * b_cell.u[cell_index]
        self.count += n_particle



class ParticleModifier:
    def __init__(self, particles):
        self.particles = particles
    
    
    def run(self, particles_out, b_particles_in, b_particles):
        while(True):
            try:
                particle_index = particles_out.pop()
            except:
                self._add_particles(b_particles_in, b_particles)
                break
            try:
                b_particle_index = b_particles_in.pop()
            except:
                self._delete_particles(particles_out)
                break
            
            self._modify(particle_index, b_particle_index, b_particles)
    
    
    def _add_particles(self, b_particles_index, b_particles):
        length = len(b_particles_index)
        x = np.zeros(length)
        y = np.zeros(length)
        u = np.zeros(length)
        v = np.zeros(length)
        w = np.zeros(length)
        eu = np.zeros(length)
        ev = np.zeros(length)
        ew = np.zeros(length)
        mass = np.ones(length)
        dia = np.ones(length)
        viscosity_index = np.ones(length)
        ref_temperature = np.ones(length)
        viscosity_coeff = np.ones(length)
        dof = np.ones(length)
        gamma = np.ones(length)
        mpv = np.ones(length)
        tag = np.ones(length, dtype = int)
        
#        for index in range(len())
        
        for index in range(length):
            x[index] = b_particles.x[b_particles_index[index]]
#            if x[index] < 0.0:
#                print "x = ", x[index], index
#            y[index] = b_particles.y[b_particles_index[index]]
#            if y[index] < 0.0:
#                print "y = ", y[index], index
            u[index] = b_particles.u[b_particles_index[index]]
            v[index] = b_particles.v[b_particles_index[index]]
            w[index] = b_particles.w[b_particles_index[index]]
            eu[index] = b_particles.eu[b_particles_index[index]]
            ev[index] = b_particles.ev[b_particles_index[index]]
            ew[index] = b_particles.ew[b_particles_index[index]]
            mass[index] = b_particles.mass[b_particles_index[index]]
            dia[index] = b_particles.dia[b_particles_index[index]]
            viscosity_index[index] = b_particles.viscosity_coeff[b_particles_index[index]]
            ref_temperature[index] = b_particles.ref_temperature[b_particles_index[index]]
            viscosity_coeff[index] = b_particles.viscosity_coeff[b_particles_index[index]]
            dof[index] = b_particles.dof[b_particles_index[index]]
            gamma[index] = b_particles.gamma[b_particles_index[index]]
            mpv[index] = b_particles.mpv[b_particles_index[index]]
            tag[index] = b_particles.tag[b_particles_index[index]]
        self.particles.x = np.concatenate((self.particles.x, x))
        self.particles.y = np.concatenate((self.particles.y, y))
        self.particles.u = np.concatenate((self.particles.u, u))
        self.particles.v = np.concatenate((self.particles.v, v))
        self.particles.w = np.concatenate((self.particles.w, w))
        self.particles.eu = np.concatenate((self.particles.eu, eu))
        self.particles.ev = np.concatenate((self.particles.ev, ev))
        self.particles.ew = np.concatenate((self.particles.ew, ew))
        self.particles.mass = np.concatenate((self.particles.mass, mass))
        self.particles.dia = np.concatenate((self.particles.dia, dia))
        self.particles.viscosity_index = np.concatenate((self.particles.viscosity_index,
                                                         viscosity_index))
        self.particles.ref_temperature = np.concatenate((self.particles.ref_temperature,
                                                         ref_temperature))
        self.particles.viscosity_coeff = np.concatenate((self.particles.viscosity_coeff,
                                                         viscosity_coeff))
        self.particles.dof = np.concatenate((self.particles.dof, dof))
        self.particles.gamma = np.concatenate((self.particles.gamma, gamma))
        self.particles.mpv = np.concatenate((self.particles.mpv, mpv))
        self.particles.tag = np.concatenate((self.particles.tag, tag))
    
    
    def _delete_particles(self, particle_out):
        self.particles.x = np.delete(self.particles.x, particle_out)
        self.particles.y = np.delete(self.particles.y, particle_out)
        self.particles.u = np.delete(self.particles.u, particle_out)
        self.particles.v = np.delete(self.particles.v, particle_out)
        self.particles.w = np.delete(self.particles.w, particle_out)
        self.particles.mass = np.delete(self.particles.mass, particle_out)
        self.particles.dia = np.delete(self.particles.dia, particle_out)
        self.particles.viscosity_coeff = np.delete(self.particles.viscosity_coeff
                                                    , particle_out)
        self.particles.dof = np.delete(self.particles.dof, particle_out)
        self.particles.gamma = np.delete(self.particles.gamma, particle_out)
        self.particles.mpv = np.delete(self.particles.mpv, particle_out)
        self.particles.tag = np.delete(self.particles.tag, particle_out)
    
    
    def _modify(self, particle_index, b_particle_index, b_particles):
        self.particles.x[particle_index] = b_particles.x[b_particle_index]
        self.particles.y[particle_index] = b_particles.y[b_particle_index]
        self.particles.u[particle_index] = b_particles.u[b_particle_index]
        self.particles.v[particle_index] = b_particles.v[b_particle_index]
        self.particles.w[particle_index] = b_particles.w[b_particle_index]
        self.particles.mass[particle_index] = b_particles.mass[b_particle_index]
        self.particles.dia[particle_index] = b_particles.dia[b_particle_index]
        self.particles.viscosity_index[particle_index] = b_particles.viscosity_index[b_particle_index]
        self.particles.ref_temperature[particle_index] = b_particles.ref_temperature[b_particle_index]
        self.particles.viscosity_coeff[particle_index] = b_particles.viscosity_coeff[b_particle_index]
        self.particles.dof[particle_index] = b_particles.dof[b_particle_index]
        self.particles.gamma[particle_index] = b_particles.gamma[b_particle_index]
        self.particles.mpv[particle_index] = b_particles.mpv[b_particle_index]
        self.particles.tag[particle_index] = b_particles.tag[b_particle_index]