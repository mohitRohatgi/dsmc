# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:06:30 2016

@author: mohit
"""

import numpy as np

class SurfaceLocator:
    def __init__(self, surface):
        self.surface = surface
        last_index = len(self.surface) - 1
        self.surface_slope = np.zeros(last_index)
        self.surface_constant = np.zeros(last_index)
        self._find_all_eqn(0, last_index - 1)
    
    
    def get_slope(self):
        return self.surface_slope
        
    
    def get_constant(self):
        return self.surface_constant
    
    
    # this function finds the slope for all the surfaces
    def _find_all_eqn(self, start, end):
        start = int(start)
        end = int(end)
        if (start == end):
            self._find_surface_slope(start)
            self._find_surface_constant(start)
        else:
            mid = (start + end) / 2
            self._find_all_eqn(start, mid)
            self._find_all_eqn(mid + 1, end)
    
    
    # this function find the slope of segment denoted by segment index.
    # segment index is the index of its first vertex in surface.
    def _find_surface_slope(self, index):
        constt1 = self.surface[index + 1][0] - self.surface[index][0]
        constt2 = self.surface[index + 1][1] - self.surface[index][1]
        if np.abs(constt1) <= 1.0e-8:
            slope = 1.0e8
        elif np.abs(constt2) <= 1.0e-8:
            slope = 0.0
        else:
            slope = (constt2 / constt1)
        self.surface_slope[index] = slope
    
    
    # this function finds the constant needed to define the eqn. of line passing
    # through the surface segment.
    def _find_surface_constant(self, index):
        constt = self.surface_slope[index] * self.surface[index][0]
        constt = self.surface[index][1] - constt
        self.surface_constant[index] = constt



class PointDetector:
    def __init__(self, surface, ref_point):
        self.points = []
        self.surface = surface
        self.ref_point = ref_point
        last_index = len(self.surface) - 1
        self.locator = SurfaceLocator(surface)
        self.ref_sign = np.zeros(last_index)
        self.surface_slope = self.locator.get_slope()
        self.surface_constant = self.locator.get_constant()
        self._find_all_sign(0, last_index - 1)
        self.flag = 0
    
    
    # this function would find cell center outside of the surface.
    # flag represent if the out of domain is considered or inside the domain.
    # flag = 0 for inside and flag = True for outside , i.e.,
    # for flag = 0, inside domain points would be returned.
    def detect_all(self, points, flag=0):
        self.flag = flag
        if len(points) > 0:
            self.points = points
            points_in = []
            self._detect(points_in, 0, int(len(self.points) - 1))
            return points_in
        else:
            return []
    
    
    # this function finds the slope for all the surfaces
    def _find_all_sign(self, start, end):
        start = int(start)
        end = int(end)
        if (start == end):
            self.ref_sign[start] = self._find_sign(start, self.ref_point)
        else:
            mid = (start + end) / 2
            self._find_all_sign(start, mid)
            self._find_all_sign(mid + 1, end)
    
    
    # this function provides recursion to above function
    def _detect(self, points_in, start, end):
        start = int(start)
        end = int(end)
        if (start == end):
            if(self._check_all_segments(start, 0, len(self.surface) - 2)):
                points_in.append(start)
        else:
            mid =  int(start + end) / 2
            self._detect(points_in, start, mid)
            self._detect(points_in, mid + 1, end)
    
    
    # this function checks the particle with all the segment
    # returns true if particle is out of domain.
    def _check_all_segments(self, point_index, start, end):
        start, end = int(start), int(end)
        if (start ==  end):
            point = self.points[point_index]
            return self._check_segment(start, point)
        else:
            mid = int(start + end) / 2
            if self.flag == 1:
                if (self._check_all_segments(point_index, start, mid) or 
                self._check_all_segments(point_index, mid + 1, end)):
                    return True
                else:
                    return False
            
            elif self.flag == 0:
                if (self._check_all_segments(point_index, start, mid) and 
                self._check_all_segments(point_index, mid + 1, end)):
                    return True
                else:
                    return False
    
    
    # segment index is the index of first vertex of segment in surface.
    # returns True if point is found on the other side of the ref_point.
    def _check_segment(self, index, point):
        sign = (self._find_sign(index, point) * self.ref_sign[index])
        if (self.flag == 0):
            if (sign < 0):
                return False
            else:
                return True
        else:
            if (sign < 0):
                return True
            else:
                return False
    
    
    # this function finds the sign of point w.r.t. the segments of 
    # surface.
    def _find_sign(self, segment_index, point):
        slope = self.surface_slope[segment_index]
        if np.abs(slope) <= 1.0e-8:
            if (point[1] < self.surface[segment_index][1]):
                return -1
            else:
                return 1
        elif np.abs(slope) >= 1.0e8:
            if (point[0] < self.surface[segment_index][0]):
                return -1
            else:
                return 1
        else:
            y = point[0] * slope
            y += self.surface_constant[segment_index]
            if (y > point[1]):
                return -1
            else:
                return 1