# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:24:06 2016

@author: mohit
"""

import unittest
from geometry import Surface
from numpy import sqrt

class test_Surface(unittest.TestCase):
    def test_triangle(self):
        ref_point = (0.0, 0.5)
        surf = Surface(ref_point)
        triangle = [((-1.0, 0.0), (0.0, 1.0)), ((0.0, 1.0), (1.0, 0.0)), 
                    ((1.0, 0.0), (-1.0, 0.0))]
        for surface in triangle:
            surf.add_surface(surface)
        
        c = 1.0 / sqrt(2.0)
        self.assertTrue(surf.get_surf_count() == 3)
        
        tangent = surf.get_surf_tangent(0)
        normal = surf.get_surf_normal(0)
        self.assertAlmostEqual(tangent[0], c, 6)
        self.assertAlmostEqual(tangent[1], c, 6)
        self.assertAlmostEqual(normal[0], -c, 6)
        self.assertAlmostEqual(normal[1], c, 6)
        
        tangent = surf.get_surf_tangent(1)
        normal = surf.get_surf_normal(1)
        self.assertAlmostEqual(tangent[0], c, 6)
        self.assertAlmostEqual(tangent[1], -c, 6)
        self.assertAlmostEqual(normal[0], c, 6)
        self.assertAlmostEqual(normal[1], c, 6)
        
        tangent = surf.get_surf_tangent(2)
        normal = surf.get_surf_normal(2)
        self.assertAlmostEqual(tangent[0], -1.0, 6)
        self.assertAlmostEqual(tangent[1], 0.0, 6)
        self.assertAlmostEqual(normal[0], 0.0, 6)
        self.assertAlmostEqual(normal[1], -1.0, 6)
    
    
    def test_box(self):
        ref_point = (0.0, 0.0)
        surf_temp = 300.0
        surf = Surface(ref_point)
        box = [((-1.0, -1.0), (-1.0, 1.0)), ((-1.0, 1.0), (1.0, 1.0)), 
                    ((1.0, 1.0), (1.0, -1.0)), ((1.0, -1.0), (-1.0, -1.0))]
        for surface in box:
            surf.add_surface(surface,surf_temp)
        
        self.assertTrue(surf.get_surf_count() == 4)
        
        tangent = surf.get_surf_tangent(0)
        normal = surf.get_surf_normal(0)
        self.assertAlmostEqual(tangent[0], 0.0, 6)
        self.assertAlmostEqual(tangent[1], 1.0, 6)
        self.assertAlmostEqual(normal[0], -1.0, 6)
        self.assertAlmostEqual(normal[1], 0.0, 6)
        
        tangent = surf.get_surf_tangent(1)
        normal = surf.get_surf_normal(1)
        self.assertAlmostEqual(tangent[0], 1.0, 6)
        self.assertAlmostEqual(tangent[1], 0.0, 6)
        self.assertAlmostEqual(normal[0], 0.0, 6)
        self.assertAlmostEqual(normal[1], 1.0, 6)
        
        tangent = surf.get_surf_tangent(2)
        normal = surf.get_surf_normal(2)
        self.assertAlmostEqual(tangent[0], 0.0, 6)
        self.assertAlmostEqual(tangent[1], -1.0, 6)
        self.assertAlmostEqual(normal[0], 1.0, 6)
        self.assertAlmostEqual(normal[1], 0.0, 6)
        
        tangent = surf.get_surf_tangent(3)
        normal = surf.get_surf_normal(3)
        self.assertAlmostEqual(tangent[0], -1.0, 6)
        self.assertAlmostEqual(tangent[1], 0.0, 6)
        self.assertAlmostEqual(normal[0], 0.0, 6)
        self.assertAlmostEqual(normal[1], -1.0, 6)



if __name__ == '__main__':
    unittest.main()