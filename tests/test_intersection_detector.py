# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 06:37:21 2015

@author: mohit
"""

import numpy as np
import unittest
import detector as dm_d
import geometry as dm_g


class test_IntersectionDetector(unittest.TestCase):
    def test_horizontal_line_quad_1(self):
        surf = [((1.0, 1.0), (6.0, 1.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (0.0, 0.5))
        x = np.linspace(1.0, 6.0, 100)
        y = np.ones(100, dtype=float) * 3.0
        detector = dm_d.IntersectionDetector(surface)
        for i in range(100):
            point = (x[i], y[i])
            truthValue = detector.detect_point(point, 0, -2.5, 1.0)
            self.assertTrue(truthValue)
            por = detector.get_por()
            intersect_time = detector.get_intersect_time()
            surface_index = detector.get_surface_index()
            if (truthValue):
                self.assertAlmostEqual(por[1], 1.0, 5)
                self.assertAlmostEqual(por[0], x[i], 5)
                self.assertAlmostEqual(intersect_time, 0.8, 5)
                self.assertEqual(surface_index, 0)
            
    
    def test_horizontal_line_quad_2(self):        
        surf = [((-1.0, 1.0), (-6.0, 1.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (-3.0, 0.5))
        x = np.linspace(-1.0, -6.0, 100)
        y = np.ones(100, dtype=float) * 3.0
        detector = dm_d.IntersectionDetector(surface)
        for i in range(100):
            point = (x[i], y[i])
            truthValue = detector.detect_point(point, 0, -2.5, 1.0)
            self.assertTrue(truthValue)
            por = detector.get_por()
            intersect_time = detector.get_intersect_time()
            surface_index = detector.get_surface_index()
            if (truthValue):
                self.assertAlmostEqual(por[1], 1.0, 5)
                self.assertAlmostEqual(por[0], x[i], 5)
                self.assertAlmostEqual(intersect_time, 0.8, 5)
                self.assertEqual(surface_index, 0)
    
    def test_horizontal_line_quad_3(self):
        surf = [((-1.0, -1.0), (-6.0, -1.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (-3.0, -0.5))
        x = np.linspace(-1.0, -6.0, 100)
        y = np.ones(100, dtype=float) * -3.0
        detector = dm_d.IntersectionDetector(surface)
        for i in range(100):
            point = (x[i], y[i])
            truthValue = detector.detect_point(point, 0, 2.5, 1.0)
            self.assertTrue(truthValue)
            por = detector.get_por()
            intersect_time = detector.get_intersect_time()
            surface_index = detector.get_surface_index()
            if (truthValue):
                self.assertAlmostEqual(por[1], -1.0, 5)
                self.assertAlmostEqual(por[0], x[i], 5)
                self.assertAlmostEqual(intersect_time, 0.8, 5)
                self.assertEqual(surface_index, 0)
    
    def test_horizontal_line_quad_4(self):
        surf = [((1.0, -1.0), (6.0, -1.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (3.0, -0.5))
        x = np.linspace(1.0, 6.0, 100)
        y = np.ones(100, dtype=float) * -3.0
        detector = dm_d.IntersectionDetector(surface)
        for i in range(100):
            point = (x[i], y[i])
            truthValue = detector.detect_point(point, 0, 2.5, 1.0)
            self.assertTrue(truthValue)
            por = detector.get_por()
            intersect_time = detector.get_intersect_time()
            surface_index = detector.get_surface_index()
            if (truthValue):
                self.assertAlmostEqual(por[1], -1.0, 5)
                self.assertAlmostEqual(por[0], x[i], 5)
                self.assertAlmostEqual(intersect_time, 0.8, 5)
                self.assertEqual(surface_index, 0)
    
    
    def test_vertical_line_quad1(self):
        surf = [((1.0, 1.0), (1.0, 6.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (-1.0, 3.0))
        y = np.linspace(1.0, 6.0, 100)
        x = np.ones(100, dtype=float) * 3.0
        detector = dm_d.IntersectionDetector(surface)
        for i in range(100):
            point = (x[i], y[i])
            truthValue = detector.detect_point(point, -2.5, 0, 1.0)
            self.assertTrue(truthValue)
            por = detector.get_por()
            intersect_time = detector.get_intersect_time()
            surface_index = detector.get_surface_index()
            if (truthValue):
                self.assertAlmostEqual(por[1], y[i], 5)
                self.assertAlmostEqual(por[0], 1.0, 5)
                self.assertAlmostEqual(intersect_time, 0.8, 5)
                self.assertEqual(surface_index, 0)
    
    
    def test_wedge_surface(self):
        point = (0.88018521325675159, 0.68087403676204228)
        vel = (92.480485409376826, 23.598134880295735)
        surf = [((0.2, 0.0), (1.0, 0.8))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (1.0, 0.0))
        detector = dm_d.IntersectionDetector(surface)
        self.assertTrue(detector.detect_point(point, vel[0], vel[1], 1.0e-5))
    
    
    def test_box(self):
        point = (0.0, 0.26082442984289667)
        vel = (-356.44993069494706, -110.17505903656894)
        surf = [((0.0, 1.0), (0.0, 0.0)), ((0.0, 0.0), (1.0, 0.0)),
            ((1.0, 0.0), (1.0, 1.0)), ((1.0, 1.0), (0.0, 1.0))]
        surface = dm_g.SurfaceGroup()
        surface.add_new_group(surf, (0.5,0.5))
        detector = dm_d.IntersectionDetector(surface)
        self.assertTrue(detector.detect_point(point, vel[0], vel[1], 1.0e-5))

if __name__ == '__main__':
    unittest.main()

#checked
# (0.82319447822331415, 0.99963676235754495) (271.97790864799919, 414.41783401180476)
# (0.88018521325675159, 0.68087403676204228) (92.480485409376826, 23.598134880295735)
# (0.93106813984527659, 0.88285720661444866) (6751.7010953089721, -8427.2055816082375)
#(0.15639635794227405, 0.22414708673128031) (11636.452923094244, -15138.619955806384)
# (0.52196234397064634, 0.91245049533196632) (20618.160664687257, -38430.654471444737)
#checked
# (0.36780581630936537, 0.16780581630936536) (459.52391938712805, 121.73626569573594) 
#checked
# (0.89304393029964202, 0.69304393029964206) (10331.019848720225, -29222.057326943985)
# (0.89167198172446993, 0.80658473212221715) (1990.8380659039306, -9500.4369738707974)
#(0.90071882770595824, 0.70381719014869371)
# (0.8332386498831329, 0.63790829842922325) (284.42023122140074, -182.54462338763454)
# (0.64171205496424, 0.480575195874131) (-780.22555221407049, -4666.5396432031794)
# (0.85327555245854436, 0.78125248781921353) (11964.554660722222, -833.13887534469814)
# (0.34017267193020262, 0.14225050771418671), (-157.69662478879243, -365.48020318720859)
# (0.93251772954806555, 0.73771578581266317) (529.93817182751604, 10.132545367751229)
# (0.34017267193020262, 0.14225050771418671) (-157.69662478879243, -365.48020318720859)