# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 17:24:05 2015

@author: mohit
"""

import unittest
import numpy as np
import dsmc_particles as dm_p
import dsmc_reflector as dm_r

class TestSpecular(unittest.TestCase):
    def test_horizontal_line(self):
        surface = [(0.0, 0.0), (5.0, 0.0)]
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.y = np.ones(100) * 0.5
        particles.u = np.ones(100) * 0.0
        particles.v = np.ones(100) * -2.0
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.x[index], 0.0)
            specular.run(index, 0, por)
            self.assertTrue(particles.y[index] == 0.0)
            self.assertTrue(particles.u[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)
    
    
    def test_horizontal_line_2(self):
        surface = [(0.0, 0.0), (5.0, 0.0)]
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.x -= 0.5
        particles.y = np.ones(100) * 0.5
        particles.u = np.ones(100) * 2.0
        particles.v = np.ones(100) * -2.0
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.x[index], 0.0)
            specular.run(index, 0, por)
            self.assertTrue(particles.y[index] == 0.0)
            self.assertTrue(particles.u[index] == 2.0)
            self.assertTrue(particles.v[index] == 2.0)
    
    
    def test_vertical_line(self):
        surface = [(0.0, 0.0), (0.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.ones(100) * 0.5
        particles.v = np.zeros(100, dtype=float)
        particles.u = np.ones(100) * -2.0
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (0.0, particles.y[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.x[index] == 0.0)
            self.assertTrue(particles.v[index] == 0.0)
            self.assertTrue(particles.u[index] == 2.0)
    
    
    def test_vertical_line_2(self):
        surface = [(0.0, 0.0), (0.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.y -= 0.5
        particles.x = np.ones(100) * 0.5
        particles.v = np.ones(100) * 2.0
        particles.u = np.ones(100) * -2.0
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (0.0, particles.y[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.x[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)
            self.assertTrue(particles.u[index] == 2.0)
    
    
    def test_45_degree_line(self):
        surface = [(1.0, 0.0), (6.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.zeros(100, dtype=float)
        particles.v = np.zeros(100, dtype=float)
        particles.u = np.ones(100, dtype=float) * 2.0
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.y[index] + 1.0, particles.y[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.u[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)
    
    
    def test_45_degree_line_2(self):
        surface = [(0.0, 0.0), (5.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * 5.0
        particles.y = np.ones(100) * 5.0
        particles.v = np.ones(100) * -2.0
        particles.u = np.zeros(100, dtype=float)
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.v[index] == 0.0)
            self.assertTrue(particles.u[index] == -2.0)
    
    
    def test_neg_45_degree_line(self):
        surface = [(0.0, 0.0), (-5.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.x = np.random.random(100) * -5.0
        particles.y = np.ones(100) * 5.0
        particles.v = np.ones(100) * -2.0
        particles.u = np.zeros(100, dtype=float)
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.v[index] == 0.0)
            self.assertTrue(particles.u[index] == 2.0)
        
        surface = ((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0))
        particles.x[0] = 0.53518289225500093
        particles.y[0] = 0.99967784164144602
        particles.u[0] = 65.12334130538774
        particles.v[0] = 141.06527432325359
        specular.__init__(surface, particles)
        por = (particles.x[index], particles.x[index])
        specular.run(0, 1, por)
        self.assertTrue(particles.v[0] == -141.06527432325359)
        self.assertTrue(particles.u[0] == 65.12334130538774)
    
    
    def test_neg_45_degree_line_2(self):
        surface = [(0.0, 0.0), (-5.0, 5.0)]
        particles = dm_p.Particles(100)
        particles.y = np.random.random(100) * 5.0
        particles.x = np.zeros(100, dtype=float)
        particles.u = np.ones(100) * -2.0
        particles.v = np.zeros(100, dtype=float)
        specular = dm_r.Specular(surface, particles)
        for index in range(100):
            por = (particles.x[index], particles.x[index])
            specular.run(index, 0, por)
            self.assertTrue(particles.u[index] == 0.0)
            self.assertTrue(particles.v[index] == 2.0)



if __name__ == '__main__':
    unittest.main()