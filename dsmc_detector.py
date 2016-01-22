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
    def _set_data(self, por, intersect_time, group_index, surface_index):
        if intersect_time < self.intersect_time or self.intersect_time == None:
            self.intersect_time = intersect_time
            self.por = por
            self.surface_index = surface_index
            self.group_index = group_index
            return True
        else:
            return False