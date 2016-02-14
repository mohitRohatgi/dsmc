# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 05:00:25 2015

@author: mohit
"""

import numpy as np
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self, x, y, dimension):
        self.x, self.y = np.meshgrid(x, y)
        self.dimension = dimension
    
    
    def plot(self, cell_x, cell_y, array, filename, title):
        dl = np.abs(self.dimension[0] - self.dimension[1]) / cell_x
        dw = np.abs(self.dimension[2] - self.dimension[3]) / cell_y
        
        x = np.arange(self.dimension[0], self.dimension[1], dl)
        y = np.arange(0.0, 1.0, dw)
        X, Y = np.meshgrid(x, y)
        plt.figure()
        CS = plt.contour(X, Y, array)
        plt.clabel(CS, inline=1, fontsize=10)
        plt.title(title)
        plt.colorbar(shrink=0.75)
        plt.savefig(filename)
        plt.clf()