# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 17:24:05 2015

@author: mohit
"""

import unittest
import numpy as np
import particles as dm_p
import reflector as dm_r
from geometry import SurfaceGroup

class TestReflector(unittest.TestCase):
    def test_horizontal_line(self):
        surf = [((0.0, 0.0), (5.0, 0.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (0.0, -0.5))
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.y = np.ones(100) * 0.5
        particles.u = np.ones(100) * 0.0
        particles.v = np.ones(100) * -2.0
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.x[index], 0.0)
            reflector.reflect(index, 1, 0, 0, por)
            self.assertTrue(particles.y[index] == 0.0)
            self.assertTrue(particles.u[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)
    
    
    def test_horizontal_line_2(self):
        surf = [((0.0, 0.0), (5.0, 0.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (0.0, -0.5))
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.x -= 0.5
        particles.y = np.ones(100) * 0.5
        particles.u = np.ones(100) * 2.0
        particles.v = np.ones(100) * -2.0
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.x[index], 0.0)
            reflector.reflect(index, 1, 0, 0, por)
            self.assertTrue(particles.y[index] == 0.0)
            self.assertTrue(particles.u[index] == 2.0)
            self.assertTrue(particles.v[index] == 2.0)
    
    
    def test_vertical_line(self):
        surf = [((0.0, 0.0), (0.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (-0.5, 0.5))
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.ones(100) * 0.5
        particles.v = np.zeros(100, dtype=float)
        particles.u = np.ones(100) * -2.0
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (0.0, particles.y[index])
            reflector.reflect(index, 1, 0, 0, por)
            self.assertTrue(particles.x[index] == 0.0)
            self.assertTrue(particles.v[index] == 0.0)
            self.assertTrue(particles.u[index] == 2.0)
    
    
    def test_vertical_line_2(self):
        surf = [((0.0, 0.0), (0.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (-0.5, 0.5))
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.y -= 0.5
        particles.x = np.ones(100) * 0.5
        particles.v = np.ones(100) * 2.0
        particles.u = np.ones(100) * -2.0
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (0.0, particles.y[index])
            reflector.reflect(index, 1, 0, 0, por)
            self.assertTrue(particles.x[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)
            self.assertTrue(particles.u[index] == 2.0)
    
    
    def test_45_degree_line(self):
        surf = [((1.0, 0.0), (6.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (2.0, 0.0))
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.zeros(100, dtype=float)
        particles.v = np.zeros(100, dtype=float)
        particles.u = np.ones(100, dtype=float) * 2.0
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.y[index] + 1.0, particles.y[index])
            reflector.reflect(index, 1, 0, 0, por)
            self.assertAlmostEqual(particles.u[index], 0.0, 6)
            self.assertAlmostEqual(particles.v[index], 2.0, 6)
    
    
    def test_45_degree_line_2(self):
        surf = [((0.0, 0.0), (5.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (1.0, 0.0))
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.y = np.ones(100) * 5.0
        particles.v = np.ones(100) * -2.0
        particles.u = np.zeros(100, dtype=float)
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            reflector.reflect(index, 1, 0, 0,por)
            self.assertAlmostEqual(particles.v[index], 0.0, 6)
            self.assertAlmostEqual(particles.u[index], -2.0, 6)
    
    
    def test_neg_45_degree_line(self):
        surf = [((0.0, 0.0), (-5.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (-1.0, 0.5))
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * -5.0
        particles.y = np.ones(100) * 5.0
        particles.v = np.ones(100) * -2.0
        particles.u = np.zeros(100, dtype=float)
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            reflector.reflect(index, 1, 0, 0, por)
            self.assertAlmostEqual(particles.v[index], 0.0, 6)
            self.assertAlmostEqual(particles.u[index], 2.0, 6)
    
    
    def test_box(self):
        surf1 = [((0.0, 0.0), (0.0, 1.0))] 
        surf2 = [((0.0, 1.0), (1.0, 1.0))] 
        surf3 = [((1.0, 1.0), (1.0, 0.0))] 
        surf4 = [((1.0, 0.0), (0.0, 0.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf1, (-1.0, 0.5))
        surface.add_new_group(surf2, (0.0, 1.5))
        surface.add_new_group(surf3, (1.2, 0.5))
        surface.add_new_group(surf4, (0.0, -0.5))
        particles = dm_p.Particles(100)
        reflector = dm_r.Reflector(surface, particles)
        particles.x[0] = 0.0
        particles.y[0] = 0.26082442984289667
        particles.u[0] = -356.44993069494706
        particles.v[0] = -110.17505903656894
        reflector.__init__(surface, particles)
        por = (0.0, 0.26082442984289667)
        reflector.reflect(0, 1, 0, 0, por)
        self.assertAlmostEqual(particles.u[0], 356.44993069494706, 5)
        self.assertAlmostEqual(particles.v[0], -110.17505903656894, 5)
    
    
    def test_neg_45_degree_line_2(self):
        surf = [((0.0, 0.0), (-5.0, 5.0))]
        surface = SurfaceGroup()
        surface.add_new_group(surf, (1.0, 0.0))
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.zeros(100, dtype=float)
        particles.u = np.ones(100) * -2.0
        particles.v = np.zeros(100, dtype=float)
        reflector = dm_r.Reflector(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            reflector.reflect(index, 1, 0, 0, por)
            self.assertAlmostEqual(particles.u[index], 0.0, 6)
            self.assertAlmostEqual(particles.v[index], 2.0, 6)



if __name__ == '__main__':
    unittest.main()