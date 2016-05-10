# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 18:58:38 2016

@author: mohit
"""
import numpy as np
import matplotlib.pyplot as plt

temperature = np.loadtxt('wedge_super_temperature.txt').view(float)
plt.contourf(temperature)
plt.colorbar()
plt.savefig('temp.png', bbox_inches='tight')
plt.clf()

number_density_0 = np.loadtxt('wedge_super_number_density_0.txt').view(float)
plt.contourf(number_density_0)
plt.colorbar()
plt.savefig('num_den_0.png', bbox_inches='tight')
plt.clf()

number_density_1 = np.loadtxt('wedge_super_number_density_1.txt').view(float)
plt.contourf(number_density_1)
plt.colorbar()
plt.savefig('num_den_1.png', bbox_inches='tight')
plt.clf()