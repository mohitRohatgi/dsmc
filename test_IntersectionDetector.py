# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 06:37:21 2015

@author: mohit
"""

import unittest
import dsmc_detector as dm_d
import dsmc_reflector as dm_r
import dsmc_particles as dm_p
import dsmc_initialiser as dm_i
import dsmc_cells as dm_c

class TestIntersectionDetector(unittest.TestCase):
    def test_points_inside(self):
        surface = [(0.0, 0.2), (1.0, 0.8)]
        detector = dm_d.IntersectionDetector(surface)
        points = [(0.5,0.51), (0.5,0.54), (0.5,0.55), (0.5, 0.6), (0.5, 0.71)]
        u = [0.2, 0.3, .4, -5.0, -2.0]
        dt = 1.0
        detector.detect_all(points, u, u, dt)
        detected_points = detector.get_points()
        for index in detected_points:
            self.assertTrue(points[index][1] >= 0.6 * points[index][0] + 0.2)
    
    
    def test_points_outside(circle):
        pass
    

if __name__ == '__main__':
    unittest.main()