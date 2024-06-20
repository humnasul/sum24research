
import numpy as np
from math import cos, pi

def compute_meshgrid_dct(grid_dct, sx, sy):
    w_x = np.zeros(sx)
    w_y = np.zeros(sy)

    for ix in range(sx):
        w_x[ix] = 2.0 * cos((pi * ix) / sx) - 2.0
    for iy in range(sy):
        w_y[iy] = 2.0 * cos((pi * iy) / sy) - 2.0
    # laplacian split into parts

    # Initialize the grid_dct
    # ! we need to avoid (ix=0,iy=0) term, since then we may have division by zero
    grid_dct[0][0] = 0.0

    # Compute grid_dct values
    for iy in range(1):
        for ix in range(1, sx):
            grid_dct[iy][ix] = 1.0 / ((w_x[ix] + w_y[iy]) * (4.0 * sx * sy))

    for iy in range(1, sy):
        for ix in range(sx):
            grid_dct[iy][ix] = 1.0 / ((w_x[ix] + w_y[iy]) * (4.0 * sx * sy))
    '''
    4*sx*sy is normalization parameter that is incorporated in grid_dct for optimization purposes
    '''

'''
sx = 10  # Replace with your desired value
sy = 10  # Replace with your desired value
grid_dct = [[0.0] * sx for _ in range(sy)]
compute_meshgrid_dct(grid_dct, sx, sy)
print(grid_dct)  # Display the computed grid_dct
'''