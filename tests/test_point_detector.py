# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:21:37 2015

@author: mohit
"""

import unittest
from dsmc_unused import PointDetector

class TestPointDetector(unittest.TestCase):
    
    def test_points_inside(self):
        points = [(0.079, -0.074), (-0.090, 0.018), (0.095, -0.107), (0.058, -0.089),
                  (1.01709521975, 1.06880980177)]
        surface = ((0.2, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.8))
        ref_point = (0.1, 0.5)
        detector = PointDetector(surface, ref_point)
        detected = detector.detect_all(points)
        for index in detected:
            self.assertTrue(points[index][0] <= 1.0 and points[index][0] >= 0.0)
            self.assertTrue(points[index][1] <= 1.0 and points[index][1] >= 0.0)
    
    
    def test_points_outside(self):
        points = [(0.079, -0.074), (-0.090, 0.018), (0.095, -0.107), (0.058, -0.089),
                  (1.01709521975, 1.06880980177), (1.1770727752927648, 0.79876967042159386), 
                (1.65720671645, 0.968803880149)]
        surface = ((0.2, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.8))
        ref_point = (0.1, 0.5)
        detector = PointDetector(surface, ref_point)
        detected = detector.detect_all(points, 1)
        for index in detected:
            self.assertTrue((points[index][0] > 1.0) or (points[index][0] < 0.0) 
                    or (points[index][1] > 1.0) or (points[index][1] < 0.0))


if __name__ == '__main__':
    unittest.main()