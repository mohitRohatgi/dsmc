# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 23:23:36 2016

@author: mohit
"""


# assuming domain is rectangular and surface is a polygon.
# inlet, outlet and surface need to be wrapped in a list even if they are just 
# represented by a single point.

class Domain:
    def __init__(self, volume, inlet=None, zero_grad=None, outlet=None):
        self.inlet = inlet
        self.zero_grad = zero_grad
        self.outlet = outlet
        self.volume = volume
        if (inlet == None or len(inlet) == 0) and (
            outlet == None or len(outlet == 0)):
            self.bool = False
        else:
            self.bool = True
    
    def get_inlet(self):
        return self.inlet
    
    def get_outlet(self):
        return self.outlet
    
    def get_zero_grad(self):
        return self.zero_grad
    
    def get_volume(self):
        return self.volume
    
    def get_count_zero_grad(self):
        return len(self.zero_grad)
    
    def set_outlet(self, outlet):
        self.outlet = outlet
    
    def set_inlet(self, inlet):
        self.inlet = inlet
    
    def set_zero_grad(self, zero_grad):
        self.zero_grad = zero_grad
    
    def is_open(self):
        return self.bool
    
    def is_inlet_on(self):
        return self.inlet != None
    
    def is_outlet_on(self):
        return self.outlet != None
    
    def is_zero_grad_on(self):
        return self.zero_grad != None



# this class of objects would do a book-keeping of different surface groups.
# Each group refers to a a group of surfaces which would constitute a specific
# object in the simulation. for e.g. suppose you have two cylinders. Each 
# cylinder would constitute a group of surfaces. Hence, here, two groups would 
# be created.
class SurfaceGroup:
    def __init__(self):
        self.group = []
    
    
    # surfaces are a list of surfaces having same temperature constituting an 
    # object or a part of an object.
    def add_new_group(self, surfaces, ref_point, surf_temp=None):
        surf = Surface()
        for surface in surfaces:
            surf.add_surface(surface, ref_point, surf_temp)
        self.group.append(surf)
    
    
    def add_surf(self, group_index, surface, surf_temp=None):
        surf_group = self.get_surf_group(group_index)
        surf_group.add_surface(surface, surf_temp)
    
    
    def get_group_count(self):
        return len(self.group)
    
    
    def get_surface(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surface(surf_index)
    
    
    def get_surf_count(self, group_index):
        return self.get_surf_group(group_index).get_surf_count()
    
    
    def get_surf_group(self, group_index):
        return self.group[group_index]
    
    
    def get_surf_vertex1(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surf_vertex1(surf_index)
    
    
    def get_surf_vertex2(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surf_vertex2(surf_index)
    
    
    def get_surf_temp(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surf_temp(surf_index)
    
    
    def get_surf_normal(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surf_normal(surf_index)
    
    
    def get_surf_tangent(self, group_index, surf_index):
        return self.get_surf_group(group_index).get_surf_tangent(surf_index)



# This class of objects would have information about the surfaces
class Surface:
    def __init__(self):
        self.surfaces = []
        self.surf_temp = None
        self.surf_tangent = {}
        self.surf_normal = {}
    
    
    def add_surface(self, surface, ref_point, surf_temp=None):
        self.surfaces.append(surface)
        self._set_surf_tangent()
        self._set_surf_normal(ref_point)
        self.surf_temp = surf_temp
    
    
    def get_surf_count(self):
        return len(self.surfaces)
    
    
    def get_surf_vertex1(self, index):
        return self.surfaces[index][0]
    
    
    def get_surf_vertex2(self, index):
        return self.surfaces[index][1]
    
    
    def get_surface(self, index):
        return self.surfaces[index]
    
    
    def get_surf_temp(self, surf_index):
        return self.surf_temp
    
    
    def get_surf_normal(self, index):
        return self.surf_normal[index]
    
    
    def get_surf_tangent(self, index):
        return self.surf_tangent[index]
    
    
    def _set_surf_tangent(self):
        index = self.get_surf_count() - 1
        vertex1 = self.get_surf_vertex1(index)
        vertex2 = self.get_surf_vertex2(index)
        tangent = (vertex2[0] - vertex1[0], vertex2[1] - vertex1[1])
        self.surf_tangent[index] = tangent
    
    
    def _set_surf_normal(self, ref_point):
        index = self.get_surf_count() - 1
        tangent = self.get_surf_tangent(index)
        vertex = self.get_surf_vertex1(-1)
        dx = ref_point[0] - vertex[0]
        dy = ref_point[1] - vertex[1]
        if dx * tangent[1] - dy * tangent[0] > 0:
            self.surf_normal[index] =  (-tangent[1], tangent[0])
        else:
            self.surf_normal[index] =  (tangent[1], -tangent[0])