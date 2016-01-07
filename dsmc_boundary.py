import numpy as np
import dsmc_cells as dm_c
import dsmc_particles as dm_p
import dsmc_detector as dm_d


class BoundaryManager:
    def __init__(self, cells, detector, domain, n_particles_in_cell,
                 mole_fraction, temperature, particles):
        self.boundary = BoundaryGenerator(cells, domain).get_boundary()
        self.particle_generator = ParticleGenerator(self.boundary, mole_fraction,
                                            temperature, n_particles_in_cell)
        self.in_detector = InParticleDetector(detector)
        self.modifier = ParticleModifier(particles)
        self.run([])
    
    
    def run(self, particles_out):
        b_particles = self.particle_generator.run()
        particles_in = self.in_detector.detect(b_particles)
        self.modifier.run(particles_in, particles_out)



# this class takes care of boundary cell generation of the domain and it would 
# also maps the cells to their adacent cells.
# adj_map is a map between a line segment and the map between boundary cell on 
# that line segment and their respective adjacent cell in the domain. 
class BoundaryGenerator:
    def __init__(self, cells, domain):
        self.ref_point = cells.get_center(0)
        self.domain = domain
        self.cells = cells
    
    
    def get_boundary(self):
        adj_map = {}
        inlet_cells = []
        zero_grad_cells = []
        if self.domain.is_inlet_on():
            for line in self.domain.get_inlet():
                inlet_cell = self.cells.generate_cell(line)
                inlet_cells.append(inlet_cell)
                adj_map[inlet_cell] = self._create_map(inlet_cell, line)
        
        if self.domain.is_zero_grad_on():
            for line in self.domain.get_zero_grad():
                zero_grad_cell = self.cells.generate_cell(line)
                zero_grad_cells.append(zero_grad_cell)
                adj_map[zero_grad_cell] = self._create_map(zero_grad_cell, line)
        
        return Boundary(inlet_cells, zero_grad_cells, adj_map)
    
    
    # this function would map each cell index of inlet/ zero gradient cell to 
    # its adjacent cell in the domain.
    # bound = boundary
    def _create_map(self, bound_cell, line):
        cell_map = {}
        for i in range(len(bound_cell.u)):
            xc, yc = bound_cell.get_center(i)
            normal_vector = self._find_normal_vector(line, (xc, yc))
            xc += normal_vector[0] * 1.01
            yc += normal_vector[1] * 1.01
            index = self.cells.find_cell_index(xc, yc)
            cell_map[i] = index
        return cell_map
    
    
    # this function returns a normal vector to the line from boundary cell
    # to the ajacent cell inside the domain.
    def _find_normal_vector(self, line, point):
        vertex1 = line[0]
        vertex2 = line[1]
        dx = vertex2[0] - vertex1[0]
        dy = vertex2[1] - vertex1[1]
        t = (point[0] - vertex1[0]) * dx + (point[1] - vertex1[1]) * dy
        t /= (dx * dx + dy * dy)
        return (vertex1[0] + t * dx - point[0], vertex1[1] + t * dy - point[1])
        


# this class would generate particles based on whether its inlet boundary or
# the zero gradient boundary. For inlet boundary, particles are generated at
# initial boundary coondition while, for zero gradient boundary, particles are
# generated at the properties of the adjacent cells.
class ParticleGenerator:
    def __init__(self, boundary, mole_fraction, temperature, n_particles_in_cell):
        pass



# this class would detect the particles that went inside the domain.
class InParticleDetector:
    def __init__(self, detector):
        pass



# this class would modify the particle array being used for computation 
# according to the boundary conditions.
class ParticleModifier:
    def __init__(self, particles):
        pass



# this class contains the information about different boundary cells, e.g., 
# inlet boundary, zero gradient boundary and oulet boundary.
# reference to all the cells being generated over the inlet line segments (can 
# be multiple line segments) are being stored as a list of references.
# same goes for the zero_grad_cells. outlet cells are not being stored as
# particles do not need to be generated in that place. Moreover, if it does 
# need to be it could always be thought of as a zero_grad_cells.
class Boundary:
    def __init__(self, inlet_cells, zero_grad_cells, adj_map):
        self.inlet_cells = inlet_cells
        self.zero_grad_cells = zero_grad_cells
        self.adj_map = adj_map
    
    
    def get_inlet_cells(self):
        return self.inlet_cells
    
    
    def get_zero_grad_cells(self):
        return self.zero_grad_cells
    
    
    # this method returns the map between the boundary cell on a particular 
    # line segment and the cell inside the domain.
    def get_adj_map(self, cell_ref):
        return self.adj_map[cell_ref]
    
    
    # this method returns the cell index of the cell inside the domain adjacent
    # to the cell being refered by the cell_ref(reference variable of cell 
    # object) and  cell index.
    def get_adj_cell_index(self, cell_ref, cell_index):
        Map = self.adj_map[cell_ref]
        return Map[cell_index]