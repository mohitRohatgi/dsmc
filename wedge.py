# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 22:49:35 2015

@author: mohit
"""

import numpy as np
import dsmc_solver as dm_s
import dsmc_particles as dm_p
import dsmc_initialiser as dm_i
import dsmc_cells as dm_c
import dsmc_output as dm_out


"""
'_col_' denotes collision while '_f_' denotes free
"""

def main():
    inlet = ((0.0, 0.0), (0.0, 1.0))
    outlet1 = ((0.0, 0.0), (0.2, 0.0))
    outlet2 = ((0.0, 1.0), (1.0, 1.0))
    outlet3 = ((1.0, 1.0), (1.0, 0.8))
    outlet = (outlet1, outlet2, outlet3)
    surface = ((0.2, 0.0), (1.0, 0.8))
    volume = 0.36
    center = (0.5,0.5)
    length = 1.0
    width = 1.0
    ensemble_sample = 100
    time_av_sample = 1000
    time = 2.0
    dof = 3.0
    mass = 66.3e-27
    viscosity_coeff = 2.117
    viscosity_index = 0.81
    mole_fraction = 1.0
    dia = 4.17e-10
    mach = [2.5, 0.0, 0.0]
    temperature = 500.0
    surface_temperature = 300.0
    ref_temperature = 273.0
    number_density = 1.699e19
    gamma = 5.0 / 3.0
    n_particles_in_cell = 20
    cell_x, cell_y = 10, 10
    ref_point = (0.1, 0.5)
    domain = dm_p.Domain(inlet, outlet, surface, volume)
    argon = dm_p.Molecule()
    argon.setup(dia, viscosity_index, mass, viscosity_coeff, dof, 1,
                ref_temperature, gamma, volume, [mole_fraction])
    gas = dm_p.Gas(array(argon), number_density, mach, temperature, surface_temperature)
    gas.setup()
    cells = dm_c.Rect_cells(cell_x, cell_y)
    cells.setup(length, width, center)
    cell_detector = dm_i.Cell_detector(cells, surface, ref_point)
    cells.cell_in = cell_detector.detect_all()
    solver = dm_s.Dsmc_solver(gas, time, ensemble_sample, time_av_sample, 
                              domain, cells)
    solver.setup()
    solver.solve()
    output = solver.extract()
    output_manager = dm_out.Output_manager()
    output_manager.show(output)



#if __name__ == '__main__':
#    main()

a = np.array([2, 3, 4, 5, 6, 7])
a = np.delete(a, list([2, 3, 4]))
print "a = ", a

a = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
index = [2, 3, 6]

new_a = np.delete(a, index)

print(new_a)