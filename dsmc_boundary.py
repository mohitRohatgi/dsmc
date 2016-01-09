import numpy as np
import dsmc_particles as dm_p


# this is the top level manager which woud be responsible for the boundary 
# condition implementation.
class BoundaryManager:
    def __init__(self, cells, detector, domain, n_particles_in_cell,
                 mole_fraction, temperature, particles, gas):
        self.gas = gas
        self.boundary = BoundaryGenerator(cells, domain).get_boundary()
        
        self.particle_generator = ParticleGenerator(self.boundary, gas, cells, 
                                    n_particles_in_cell, particles.get_n_eff())
        
        self.in_detector = InParticleDetector(cells)
        self.modifier = ParticleModifier(particles)
        self.run([])
    
    
    def run(self, particles_out):
        self.particle_generator.run()
        inlet_particles = self.particle_generator.get_inlet_particles()
        zero_grad_particles = self.particle_generator.get_zero_grad_particles()
        self.in_detector.detect(inlet_particles)
        self.in_detector.detect(zero_grad_particles)
        self.modifier.run(self.in_detector, particles_out)



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
    def __init__(self, boundary, gas, cells, n_particles_in_cell, n_eff):
        k = 1.3806488e-23
        c = np.sqrt(gas.get_gamma() * gas.get_temperature() * k / gas.get_mass())
        self.n_particles_in_cell = n_particles_in_cell
        self.boundary = boundary
        self.n_eff = n_eff
        self.gas = gas
        self.cells = cells
        n_particles = self._find_inlet_particles()
        self.inlet_particles = dm_p.Particles(n_particles, self.n_eff)
        self.inlet_particles.setup(gas.get_mol_frac(), gas.get_species(),
                                   gas.get_mpv())
        self.zero_grad_particles = dm_p.Particles(1)
        self.mean_inlet_velx = self.gas.get_mach_x() * c 
        self.mean_inlet_vely = self.gas.get_mach_x() * c
        self.mean_inlet_velz = self.gas.get_mach_x() * c
    
    
    def run(self):
        start = 0
        for inlet_cell in self.boundary.get_inlet_cells():
            start = self._init_inlet_particles(inlet_cell, start)
        
        n_particles = self._find_zero_grad_particles()
        self.zero_grad_particles = dm_p.Particles(n_particles, self.n_eff)
        start = 0
        for zero_grad_cell in self.boundary.get_zero_grad_cells():
            start = self._init_zero_grad_particles(zero_grad_cell, start)
    
    
    def get_inlet_particles(self):
        return self.inlet_particles
    
    
    def get_zero_grad_particles(self):
        return self.zero_grad_particles
    
    
    # this function initialises the inlet particles.
    def _init_inlet_particles(self, inlet_cell, start):
        offset = start
        for cell_index in range(len(inlet_cell.get_temperature())):
            i = 0
            while i < self.n_particles_in_cell:
                x, y = self._find_rand_loc(inlet_cell, cell_index)
                self.inlet_particles.set_x(x, offset + i)
                self.inlet_particles.set_y(y, offset + i)
                
                v1, v2, v3 = self._find_rand_vel(self.inlet_particles, 
                                                 offset + i)
                self.inlet_particles.set_velx(v1 + self.mean_inlet_velx)
                self.inlet_particles.set_vely(v2 + self.mean_inlet_vely)
                self.inlet_particles.set_velz(v3 + self.mean_inlet_velz)
                i += 1
            offset += i
        
        return offset
    
    
    def _init_zero_grad_particles(self, zero_grad_cell, start):
        offset = start
        for cell_index in range(len(zero_grad_cell.get_temperature())):
            i = 0
            index = self.boundary.get_adj_cell_index(zero_grad_cell, cell_index)
            while i < len(self.cells.get_particles_inside(index)):
                x, y = self._find_rand_loc(zero_grad_cell, cell_index)
                self.zero_grad_particles.set_x(x, offset + i)
                self.zero_grad_particles.set_y(y, offset + i)
                
                v1, v2, v3 = self._find_rand_vel(self.zero_grad_particles, 
                                                 offset + i)
                self.zero_grad_particles.set_velx(v1 + self.cells.get_velx(cell_index))
                self.zero_grad_particles.set_vely(v2 + self.cells.get_vely(cell_index))
                self.zero_grad_particles.set_velz(v3 + self.cells.get_velz(cell_index))
                i += 1
            offset += i
        
        return offset
    
    
    # this function generates a random location inside the cell provided.
    def _find_rand_loc(self, cells, cell_index):
        xc, yc = cells.get_center(cell_index)
        
        x = ((2.0 * np.random.rand() - 1) / 2.0 * 
            cells.get_cell_length(cell_index) + xc)
        
        y = ((2.0 * np.random.rand() - 1) / 2.0 *
            cells.get_cell_width(cell_index) + yc)
        
        return (x, y)
    
    
    # this function finds and sets the particle velocity.
    # mean vel is the velocity vector of the cell in which particle is.
    def _find_rand_vel(self, particles, index):
        v1 = np.random.normal(0.0, 0.5) * particles.get_mpv(index)
        v2 = np.random.normal(0.0, 0.5) * particles.get_mpv(index)
        v3 = np.random.normal(0.0, 0.5) * particles.get_mpv(index)
        return (v1, v2, v3)
    
    
    # this is a helper function to find the total number of particles needed 
    # to be introduced inside all the inlet cells.
    def _find_inlet_particles(self):
        count = 0
        for inlet_cell in self.boundary.get_inlet_cells():
            count += len(inlet_cell.get_temperature())
        
        return count * self.n_particles_in_cell
    
    
    # this is also a helper function to find the total no. of particles needed
    # for zero gradient boundaries.
    def _find_zero_grad_particles(self):
        count = 0
        for zero_grad_cells in self.boundary.get_zero_grad_cells():
            for index in range(len(zero_grad_cells.get_temperature())):
                index = self.boundary.get_adj_cell_index(zero_grad_cells, index)
                count += len(self.cells.get_particles_inside(index))
        
        return count



# this class would detect the particles that went inside the domain.
class InParticleDetector:
    def __init__(self, detector):
        pass



# this class would modify the particle array being used for computation 
# according to the boundary conditions.
class ParticleModifier:
    def __init__(self, particles):
        pass
    
    
    def run(in_detector, particles_out):
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