import numpy as np
import pyfftw

def compute_poisson_rec(input_image, output_image, binary_M, grid_dct, sx, sy, Lap, Lap_dct, plan_forward, plan_backward):
    Lap[0, 0] = -input_image[0, 0] + input_image[0, 1]

    if binary_M[0, 0] != 0:
                Lap[0, 0] += input_image[1, 0] - input_image[0, 0]
    for ix in range(1, sx - 1):
        Lap[0, ix] = input_image[0, ix-1] - 2.0*input_image[0, ix] + input_image[0, ix+1]
        if binary_M[0, ix]!=0: 
              Lap[0, ix] += input_image[1, ix] - input_image[0, ix]
    
    Lap[0, sx-1] = input_image[0, sx-2] - input_image[0, sx-1]
    if binary_M[0, sx-1]!=0: 
              Lap[0, sx-1] += input_image[1, sx-1] - input_image[0, sx-1]

    return

