import numpy as np
import particles as dm_p


# this is the top level manager which woud be responsible for the boundary 
# condition implementation.
class BoundaryManager:
    def __init__(self, cells, domain, gas, particles, n_particles_in_cell):
        self.gas = gas
        self.boundary = BoundaryGenerator(cells, domain).get_boundary()
        
        self.particle_generator = ParticleGenerator(self.boundary, gas, cells, 
                                    n_particles_in_cell, particles.get_n_eff())
        
        self.in_detector = InParticleDetector(cells)
        self.modifier = ParticleModifier(particles)
    
    
    def run(self, particles_out, dt):
        self.particle_generator.run()
        inlet_particles = self.particle_generator.get_inlet_particles()
        zero_grad_particles = self.particle_generator.get_zero_grad_particles()
        self.in_detector.detect(inlet_particles, dt)
        self.in_detector.detect(zero_grad_particles, dt)
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
    def __init__(self, cells):
        self.cells = cells
        self.particle_map = {}
    
    
    def detect(self, particles, dt):
        particles_out = []
        for index in range(len(particles.get_x())):
            x = particles.get_x(index) + particles.get_velx(index) * dt
            y = particles.get_x(index) + particles.get_vely(index) * dt
            if self.cells.is_inside_domain(x, y):
                particles_out.append(index)
        if len(particles_out) > 0:
            self.particle_map[particles] = particles_out
    
    
    def get_particles_in(self, particles):
        return self.particle_map[particles]
    
    
    def get_particles(self):
        return self.particle_map.keys()



# this class would modify the particle array being used for computation 
# according to the boundary conditions.
class ParticleModifier:
    def __init__(self, particles):
        self.particles = particles
    
    
    def run(self, in_detector, particles_out):
        count = 0.0
        for b_particles in in_detector.get_particles():
            particles_in = in_detector.get_particles_in(b_particles)
            count += len(particles_in)
            while (len(particles_in) > 0):
                particle_in = particles_in.pop()
                try:
                    particle_out = particles_out.pop()
                    self._modify_particles(b_particles, particle_in, particle_out)
                except:
                    particles_in.append(particle_in)
                    self._add_particles(b_particles, particles_in)
                    break
    
    
    def _modify_particles(self, b_particles, particle_in, particle_out):
        self.particles.x[particle_out] = b_particles.x[particle_in]
        self.particles.y[particle_out] = b_particles.y[particle_in]
        self.particles.u[particle_out] = b_particles.u[particle_in]
        self.particles.v[particle_out] = b_particles.v[particle_in]
        self.particles.w[particle_out] = b_particles.w[particle_in]
        self.particles.mpv[particle_out] = b_particles.mpv[particle_in]
        self.particles.tag[particle_out] = b_particles.tag[particle_in]
    
    
    def _add_particles(self, b_particles, particles_in):
        length = len(particles_in)
        x = np.zeros(length)
        y = np.zeros(length)
        u = np.zeros(length)
        v = np.zeros(length)
        w = np.zeros(length)
        eu = np.zeros(length)
        ev = np.zeros(length)
        ew = np.zeros(length)
        mpv = np.ones(length)
        tag = np.ones(length, dtype = int)
        
        # energy would be computed when required. So, there is no point in 
        # saving it here.
        for index in range(length):
            x[index] = b_particles.x[index]
            y[index] = b_particles.y[index]
            u[index] = b_particles.u[index]
            v[index] = b_particles.v[index]
            w[index] = b_particles.w[index]
            mpv[index] = b_particles.mpv[index]
            tag[index] = b_particles.tag[index]
        self.particles.x = np.concatenate((self.particles.x, x))
        self.particles.y = np.concatenate((self.particles.y, y))
        self.particles.u = np.concatenate((self.particles.u, u))
        self.particles.v = np.concatenate((self.particles.v, v))
        self.particles.w = np.concatenate((self.particles.w, w))
        self.particles.eu = np.concatenate((self.particles.eu, eu))
        self.particles.ev = np.concatenate((self.particles.ev, ev))
        self.particles.ew = np.concatenate((self.particles.ew, ew))
        self.particles.mpv = np.concatenate((self.particles.mpv, mpv))
        self.particles.tag = np.concatenate((self.particles.tag, tag))
                



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
    # to the cell being referred by the cell_ref(reference variable of cell 
    # object) and  cell index.
    def get_adj_cell_index(self, cell_ref, cell_index):
        Map = self.adj_map[cell_ref]
        return Map[cell_index]