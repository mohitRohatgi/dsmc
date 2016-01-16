# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 14:48:54 2015

@author: mohit
"""

import numpy as np


# this checks whether the particular particle would intersect with the surface.
# If it intersects, it stores the por, intersecting time and the index of the
# surface it gets reflected from.
class IntersectionDetector:
    def __init__(self, surf_group):
        self.surf_group = surf_group
        self.dt = 0.0
        self.group_index = None
        self.surface_index = None
        self.intersect_time = None
        self.por = None
        self.FUZZ = 1.0e-6
        self.FUZZ_SQ = self.FUZZ * self.FUZZ
    
    
    def get_intersect_time(self):
        return self.intersect_time
    
    
    def get_por(self):
        return self.por
    
    
    def get_surf_group_index(self):
        return self.group_index
    
    
    def get_surface_index(self):
        return self.surface_index
    
    
    def detect_point(self, point, u, v, dt):
        self.dt = dt
        self.intersect_time = None
        self.por = None
        self.group_index = None
        self.surface_index = None
        return self._detect_group(point, u, v, 0,
                               self.surf_group.get_group_count() - 1)
    
#    # por represents the por of the particle and surface being considered while
#    # self.por is the surface it would reflect after checking with all the surfaces.
#    def _detect_point(self, point, u, v, start, end):
#        start, end = int(start), int(end)
#        if start == end:
#            por = self._find_por(point, u, v, start)
#            # checking for parallelism and finite line segment case.
#            if por == None:
#                return False
#            else:
#                intersect_time = self._find_intersect_time(point, u, v, por)
#                # checking if intersect time gets negative signifying no
#                # intersection of particle trajectory and surface.
#                if intersect_time < 0.0:
#                    return False
#                # accounting for the fuzz
#                elif self._nearby(point, start):
#                    return True
#                elif intersect_time <= self.dt:
#                    return self._set_data(por, intersect_time, start)
#                else:
#                    return False
#        else:
#            mid = int((start + end) / 2)
#            truthValue1 = self._detect_point(point, u, v, start, mid)
#            truthValue2 = self._detect_point(point, u, v, mid + 1, end)
#            if truthValue1 or truthValue2:
#                return True
#            else:
#                return False
    
    
    # this function detects the group of surfaces from which particle would
    # reflect.
    def _detect_group(self, point, u, v, start, end):
        if (start == end):
            surf_count = self.surf_group.get_surf_count(start)
            return self._detect_surf(point, u, v, start, 0, surf_count - 1)
        else:
            mid = int((start + end) / 2)
            return (self._detect_group(point, u, v, start, mid) or 
                    self._detect_group(point, u, v, mid + 1, end))
    
    
    # this function detects the surface from which particle would reflect
    def _detect_surf(self, point, u, v, group_index, start, end):
        if start == end:
            return self._check_surf(point, u, v, group_index, start)
        else:
            mid = int((start + end) / 2)
            return (self._detect_surf(point, u, v, group_index, start, mid) or
                    self._detect_surf(point, u, v, group_index, mid + 1, end))
    
    
    # this function checks the particle for reflection from the surface 
    # specified by the group index and the surface index
    # por represents the por of the particle and surface being considered while
    # self.por is the por of the particle and surface it would actually reflect
    # which is determined after checking over all groups and their surfaces.
    def _check_surf(self, point, u, v, group_index, surf_index):
        vertex1 = self.surf_group.get_surf_vertex1(group_index, surf_index)
        tangent = self.surf_group.get_surf_tangent(group_index, surf_index)
        por = self._find_por(point, u, v, vertex1, tangent)
        # checking for parallelism and finite line segment case.
        if por == None:
            return False
        else:
            intersect_time = self._find_intersect_time(point, u, v, por)
            # checking if intersect time gets negative signifying no
            # intersection of particle trajectory and surface.
            if intersect_time < 0.0:
                return False
            # accounting for the fuzz
            elif self._nearby(point, group_index, surf_index):
                return True
            elif intersect_time <= self.dt:
                return self._set_data(por, intersect_time,
                                      group_index, surf_index)
            else:
                return False
    
    
    # this function checks if particle is inside the fuzz area considered to
    # account for floating point no.
    def _nearby(self, point, group_index, surf_index):
        tangent = self.surf_group.get_surf_tangent(group_index, surf_index)
        vertex1 = self.surf_group.get_surf_vertex1(group_index, surf_index)
        t = self._find_t(point, vertex1, tangent)
        dx = vertex1[0] - point[0] + t * tangent[0]
        dy = vertex1[1] - point[1] + t * tangent[1]
        ds = dx * dx + dy * dy
        if (ds <= self.FUZZ_SQ):
            por = (vertex1[0] + t * tangent[0], vertex1[1] + t * tangent[1])
            return self._set_data(por, self.FUZZ * 1.0e-3,
                                  group_index, surf_index)
        else:
            return False
    
    
    # this function is a helper function for the _nearby function above.
    # it finds the parameter to define the point which is perpendicular to particle
    # and lies on the surface. Using this parameter length of perpendicular is
    # found.
    def _find_t(self, point, vertex1, tangent):
        t = (point[0] - vertex1[0]) * tangent[0]
        t += (point[1] - vertex1[1]) * tangent[1]
        t /= tangent[0] ** 2.0 + tangent[1] ** 2.0
        
        return t
    
    
    # por = None means that the particle would not intersect.
    # finite line segment case has also been considered here.
    def _find_por(self, point, u, v, vertex1, tangent):
        constt = tangent[1] * u - v * tangent[0]
        # parallelism case
        if constt == 0.0:
            return None
        else:
            t = u * point[1] + v * (vertex1[0] - point[0])
            t -= u * vertex1[1]
            t /= constt
            # finite line segment case.
            if (t < -self.FUZZ or t > (1 + self.FUZZ)):
                return None
            x = vertex1[0] + t * tangent[0]
            y = vertex1[1] + t * tangent[1]
            return (x, y)
    
    
    # this function would give the intersection time.
    def _find_intersect_time(self, point, u, v, por):
        vel = np.sqrt(u ** 2.0 + v ** 2.0)
        dx, dy = por[0] - point[0], por[1] - point[1]
        dt =  (u * dx + v * dy) / vel / vel
        return dt
    
    
    # this function would check if the intersection time is less than the minimum
    # intersection time detected uptill now. Depending upon which it would 
    # modify the minimum intersection time. If modified would return True
    # signifying a surface having lesser intersection time is detected.
    def _set_data(self, por, intersect_time, surface_index, group_index):
        if intersect_time < self.intersect_time or self.intersect_time == None:
            self.intersect_time = intersect_time
            self.por = por
            self.surface_index = surface_index
            self.group_index = group_index
            return True
        else:
            return False



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