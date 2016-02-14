# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 23:23:56 2016

@author: mohit
"""

import unittest
import numpy as np
from dsmc_boundary import BoundaryGenerator
from dsmc_cells import RectCells
from dsmc_geometry import Domain


class TestCreateCell(unittest.TestCase):
    def test_box(self):
        cells = RectCells(10, 10, 4.0, 2.0, (0.0, 0.0), 1)
        line = ((-2.0, -1.0), (-2.0, 0.0))
        b_cells = cells.generate_cell(line)
        self.assertEqual(b_cells.n_y, 5)
        self.assertEqual(b_cells.cell_length, cells.cell_length)
        self.assertEqual(b_cells.cell_width, cells.cell_width)
        
        line = ((-2.0, -1.0), (0.0, -1.0))
        b_cells = cells.generate_cell(line)
        self.assertEqual(b_cells.n_y, 1)
        self.assertEqual(b_cells.n_x, 5)
        self.assertEqual(b_cells.cell_length, cells.cell_length)
        self.assertEqual(b_cells.cell_width, cells.cell_width)
        
        box = [((-2.0, -1.0), (-2.0, 1.0)), ((-2.0, 1.0), (2.0, 1.0)),
               ((2.0, 1.0), (2.0, -1.0)), ((2.0, -1.0), (-2.0, -1.0))]
        
        for line in box:
            b_cells = cells.generate_cell(line)
            if np.abs(line[0][0] - line[1][0]) < 1.0e-6:
                self.assertEqual(b_cells.n_y, cells.n_y)
                self.assertEqual(b_cells.n_x, 1)
            else:
                self.assertEqual(b_cells.n_x, cells.n_x)
                self.assertEqual(b_cells.n_y, 1)
            self.assertEqual(b_cells.cell_length, cells.cell_length)
            self.assertEqual(b_cells.cell_width, cells.cell_width)


class TestBoundaryGenerator(unittest.TestCase):
    def test_box(self):
        inlet = [((-2.0, -1.0), (-2.0, 1.0))]
        zero_grad = [((-2.0, -1.0), (0.0, -1.0)), ((-2.0, 1.0), (2.0, 1.0))]
        outlet = []
        volume = 1.5
        domain = Domain(volume, inlet, zero_grad, outlet)
        cells = RectCells(10, 10, 4.0, 2.0, (0.0, 0.0), 1)
        boundary = BoundaryGenerator(cells, domain).get_boundary()
        
        inlet_cells = boundary.get_inlet_cells()
        zero_grad_cells = boundary.get_zero_grad_cells()
        
        self.assertEqual(inlet_cells[0].cell_length, cells.cell_length)
        self.assertEqual(inlet_cells[0].n_y, cells.n_y)
        self.assertEqual(inlet_cells[0].n_x, 1)
        self.assertAlmostEqual(inlet_cells[0].get_center(4)[0], -2.2, 6)
        self.assertAlmostEqual(inlet_cells[0].get_center(4)[1], -0.1, 6)
        
        self.assertEqual(zero_grad_cells[0].cell_length, cells.cell_length)
        self.assertEqual(zero_grad_cells[0].n_y, 1)
        self.assertEqual(zero_grad_cells[0].n_x, cells.n_x / 2)
        self.assertAlmostEqual(zero_grad_cells[0].get_center(4)[0], -0.2, 6)
        self.assertAlmostEqual(zero_grad_cells[0].get_center(4)[1], -1.1, 6)
        
        self.assertEqual(zero_grad_cells[1].cell_length, cells.cell_length)
        self.assertEqual(zero_grad_cells[1].n_y, 1)
        self.assertEqual(zero_grad_cells[1].n_x, cells.n_x)
        self.assertAlmostEqual(zero_grad_cells[1].get_center(4)[0], -0.2, 6)
        self.assertAlmostEqual(zero_grad_cells[1].get_center(4)[1], 1.1, 6)


if __name__ == '__main__':
    unittest.main()