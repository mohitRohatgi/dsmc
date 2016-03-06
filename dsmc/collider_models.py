# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 12:16:14 2016

@author: mohit
"""

import numpy as np

class VhsCollider():
    def __init__(self, particles):
        self.particles = particles
        self.collision_pairs = []
        self.rel_speed = []

    
    def run(self, collision_pairs, rel_speed):
        self.collision_pairs = collision_pairs
        self.rel_speed = rel_speed
        if len(self.collision_pairs) > 0:
            self._vhs(0, len(self.collision_pairs) - 1)
            print ("collided in the time step", " # of particles collided = ", 
                   len(self.collision_pairs), " # of particles = ", 
                    len(self.particles.get_x()))
        else: 
            print "no collision in the time step"
    
    
    def _vhs(self, start, end):
        if (start == end):
            self._vhs_pair(start)
        else:
            mid = int((start + end) / 2)
            self._vhs(start, mid)
            self._vhs(mid + 1, end)
    
    
    def _vhs_pair(self, pair_index):
        pair_index = int(pair_index)
        avg_vel = self._find_pair_avg_vel(pair_index)
        rel_vel = self._find_vhs_post(self.rel_speed[pair_index])
        vel11, vel21 = self._find_velocity(avg_vel[0], rel_vel[0])
        vel12, vel22 = self._find_velocity(avg_vel[1], rel_vel[1])
        vel13, vel23 = self._find_velocity(avg_vel[2], rel_vel[2])
        index1 = self.collision_pairs[pair_index][0]
        index2 = self.collision_pairs[pair_index][1]
        particle = self.particles
        particle.set_velx(vel11, index1)
        particle.set_vely(vel12, index1)
        particle.set_velz(vel13, index1)
        particle.set_velx(vel11, index2)
        particle.set_vely(vel12, index2)
        particle.set_velz(vel13, index2)
    
    
    def _find_vhs_post(self, rel_speed):
        random_multiplication_factor = 2.0 * np.random.random() - 1.0
        post_u_rel = random_multiplication_factor * rel_speed
        Random_yz_angle = 2.0 * np.pi * np.random.random()
        yz_speed = (1.0 - random_multiplication_factor ** 2)
        yz_speed = np.sqrt(yz_speed) * rel_speed
        post_v_rel = yz_speed * np.sin(Random_yz_angle)
        post_w_rel = yz_speed * np.cos(Random_yz_angle)
        return (post_u_rel, post_v_rel, post_w_rel)
    
    
    def _find_pair_avg_vel(self, index):
        particle = self.particles
        index1 = self.collision_pairs[index][0]
        index2 = self.collision_pairs[index][1]
        u_avg = (particle.get_velx(index1) + particle.get_velx(index2)) * 0.5
        v_avg = (particle.get_vely(index1) + particle.get_vely(index2)) * 0.5
        w_avg = (particle.get_velz(index1) + particle.get_velz(index2)) * 0.5
        return (u_avg, v_avg, w_avg)
    
    
    def _find_velocity(self, avg_vel, rel_vel):
        vel1 = (2.0 * avg_vel + rel_vel) / 2.0
        vel2 = (2.0 * avg_vel - rel_vel) / 2.0
        return (vel1, vel2)