# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 02:02:44 2016

@author: mohit
"""

import unittest
import dsmc_boundary as dm_b
from dsmc_cells import RectCells
from dsmc_geometry import Domain
import dsmc_particles as dm_p

class test_ParticleGenerator(unittest.TestCase):
    def test_box(self):
        inlet = [((-2.0, -1.0), (-2.0, 1.0))]
        zero_grad = [((-2.0, -1.0), (0.0, -1.0)), ((-2.0, 1.0), (2.0, 1.0))]
        outlet = []
        volume = 8.0
        domain = Domain(volume, inlet, zero_grad, outlet)
        cells = RectCells(10, 10, 4.0, 2.0, (0.0, 0.0), 1)
        for i in range(100):
            p = []
            for j in range(i):
                p.append(j)
            cells.particles_inside[i] = p
        boundary = dm_b.BoundaryGenerator(cells, domain).get_boundary()
        
        dof = 3.0
        mass = 66.3e-27
        viscosity_coeff = 2.117
        viscosity_index = 0.81
        mole_fraction = [0.5, 0.5]
        dia = 4.17e-10
        mach = [0.0, 0.0, 0.0]
        temperature = 500.0
        ref_temperature = 273.0
        number_density = 1.699e19
        gamma = 5.0 / 3.0
        volume = 1.0
        argon = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 0,
                    ref_temperature, gamma, volume, number_density)
        argon1 = dm_p.Molecules(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                    ref_temperature, gamma, volume, number_density)
        gas = dm_p.Gas([argon, argon1], mole_fraction, mach, temperature)
        gas.setup()
        particle_generator = dm_b.ParticleGenerator(boundary, gas, cells, 10, 1.0)
        particle_generator.run()
        inlet_particles = particle_generator.get_inlet_particles()
        zero_grad_particles = particle_generator.get_zero_grad_particles()
        
        offset = -10
        for inlet_cells in boundary.get_inlet_cells():
            for index in range(len(inlet_cells.get_temperature())):
                i = 0
                offset += 10
                length = inlet_cells.get_cell_length()
                width = inlet_cells.get_cell_width()
                while i < 10:
                    x = inlet_particles.get_x(offset + i)
                    y = inlet_particles.get_y(offset + i)
                    xc, yc = inlet_cells.get_center(index)
                    self.assertTrue(x <= xc + length)
                    self.assertTrue(y <= yc + width)
                    self.assertTrue(x >= xc - length)
                    self.assertTrue(y >= yc - width)
                    i += 1
        self.assertTrue(offset + 10 == len(inlet_particles.get_x()))
        
        offset = 0
        for zero_grad_cells in boundary.get_zero_grad_cells():
            for index in range(len(zero_grad_cells.get_temperature())):
                cell_index = boundary.get_adj_cell_index(zero_grad_cells, index)
                i = 0
                length = zero_grad_cells.get_cell_length(index)
                width = zero_grad_cells.get_cell_width(index)
                while i < len(cells.get_particles_inside(cell_index)):
                    x = zero_grad_particles.get_x(offset + i)
                    y = zero_grad_particles.get_y(offset + i)
                    xc, yc = zero_grad_cells.get_center(index)
                    self.assertTrue(x <= xc + length)
                    self.assertTrue(y <= yc + width)
                    self.assertTrue(x >= xc - length)
                    self.assertTrue(y >= yc - width)
                    i += 1
                offset += i



if __name__ == '__main__':
    unittest.main()