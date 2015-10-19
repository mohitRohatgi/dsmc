# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 14:48:54 2015

@author: mohit
"""

import numpy as np


class NormalFinder:
    def __init__(self, surface):
        self.surface = surface
        self.surface_normal = []
        self.surface_tangent = []
        self._find_all()
    
    
    def get_normal(self):
        return self.surface_normal
    
    
    def get__tangent(self):
        return self.surface_tangent
    
    
    def _find_all(self):
        self._find(0, len(self.surface) - 2)
    
    
    def _find(self, start, end):
        start, end = int(start), int(end)
        if start == end:
            dx = self.surface[start + 1][0] - self.surface[start][0]
            dy = self.surface[start + 1][1] - self.surface[start][1]
            ds = np.sqrt(dx * dx + dy * dy)
            dx /= ds
            dy /= ds
            normal = (-dy, dx)
            tangent = (dx, dy)
            self.surface_normal.append(normal)
            self.surface_tangent.append(tangent)
        else:
            mid = int((start + end) / 2)
            self._find(start, mid)
            self._find(mid + 1, end)



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



# this checks whether the particular particle would intersect with the surface.
# If it intersects, it stores the por, intersecting time and the index of the
# surface it gets reflected from. 
class IntersectionDetector:
    def __init__(self, surface):
        self.surface = surface
        self.dt = 0.0
        self.normal_finder = NormalFinder(surface)
        self.surface_tangent = self.normal_finder.get__tangent()
        self.surface_index = None
        self.intersect_time = None
        self.por = None
    
    
    def get_intersect_time(self):
        return self.intersect_time
    
    
    def get_por(self):
        return self.por
    
    
    def get_surface_index(self):
        return self.surface_index
    
    
    def detect_point(self, point, u, v, dt):
        self.dt = dt
        self.intersect_time = None
        self.por = None
        self.surface_index = None
        if (self._detect_point(point, u, v, 0, len(self.surface) - 2)):
            return True
        else:
            return False
    
    
    # por represents the por of the particle and surface being considered while
    # self.por is the surface it would reflect after checking with all the surfaces.
    def _detect_point(self, point, u, v, start, end):
        start, end = int(start), int(end)
        if start == end:
            por = self._find_por(point, u, v, start)
            print "por = ", por 
            if por == None and self.por == None:
                return False
#            elif (np.sqrt((por[0] - point[0]) ** 2.0 + 
#                        (por[1] - point[1]) ** 2.0) < 1e-6):
#                self.intersect_time = 0.0
#                self.por = por
#                self.surface_index = start
#                print "yes"
#                return True
            else:
                intersect_time = self._find_intersect_time(point, u, v, por)
                print "intersect time = ", intersect_time
                if intersect_time >= -1.0e-3 * self.dt:
                    if intersect_time <= self.dt:
                        self._set_data(por, intersect_time, start)
                        return True
                    elif np.abs(intersect_time - self.dt) < self.dt * 1.0e-2:
                        self.por = por
                        self.intersect_time = self.dt * (0.99)
                        self.surface_index = start
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            mid = int((start + end) / 2)
            truthValue1 = self._detect_point(point, u, v, start, mid)
            truthValue2 = self._detect_point(point, u, v, mid + 1, end)
            if truthValue1 or truthValue2:
                return True
            else:
                return False
    
    
    # por = None means that the particle is moving parallel to surface.
    def _find_por(self, point, u, v, index):
        constt = self.surface_tangent[index]
        constt = constt[1] * u - v * constt[0]
        if constt == 0.0:
            return None
        else:
            t = u * point[1] + v * (self.surface[index][0] - point[0])
            t -= u * self.surface[index][1]
            t /= constt
            x = self.surface[index][0] + t * self.surface_tangent[index][0]
            y = self.surface[index][1] + t * self.surface_tangent[index][1]
            return (x, y)
    
    
    def _find_intersect_time(self, point, u, v, por):
        vel = np.sqrt(u ** 2.0 + v ** 2.0)
        dx, dy = por[0] - point[0], por[1] - point[1]
        dt =  (u * dx + v * dy) / vel / vel
        return dt
    
    
    def _set_data(self, por, intersect_time, surface_index):
        if intersect_time < self.intersect_time or self.intersect_time == None:
            self.intersect_time = intersect_time
            self.por = por
            self.surface_index = surface_index