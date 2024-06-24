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

    for iy in range (1, sy-1):
        Lap[iy, 0] = -input_image[iy, 0] + input_image[iy, 1]
        if binary_M[iy, 0] != 0:
               Lap[iy, 0] += input_image[iy+1, 0] - input_image[iy, 0]
        if binary_M[iy-1, 0] != 0:
               Lap[iy, 0] += input_image[iy-1, 0] - input_image[iy, 0]

        for ix in range (1, sx-1):
                Lap[iy, ix] = input_image[iy, ix-1] - 2.0*input_image[iy, ix] + input_image[iy, ix+1]; 
                if binary_M[iy, ix] != 0:
                    Lap[iy, ix] += input_image[iy+1, ix] - input_image[iy, ix]
                if binary_M[iy-1, ix] != 0:
                    Lap[iy, ix] += input_image[iy-1, ix] - input_image[iy, ix]

        Lap[iy, sx-1] = input_image[iy][sx-2] - input_image[iy][sx-1];
        if(binary_M[iy  ][sx-1]!=0) Lap[iy][sx-1] += input_image[iy+1][sx-1] - input_image[iy][sx-1];
        if(binary_M[iy-1][sx-1]!=0) Lap[iy][sx-1] += input_image[iy-1][sx-1] - input_image[iy][sx-1];
    
    return

