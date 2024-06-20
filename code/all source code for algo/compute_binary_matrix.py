def compute_binary_matrix(input_image, binary_M, thresh_x, thresh_y, Tmin, Tmax, sx, sy):
    for iy in range(sy):
        for ix in range(sx):
            if ix < (sx - 1):
                gx = input_image[iy][ix + 1] - input_image[iy][ix]
            else:
                gx = 0.0
            if iy < (sy - 1):
                gy = input_image[iy + 1][ix] - input_image[iy][ix]
            else:
                gy = 0.0
            if (abs(gx) > thresh_x) or (abs(gy) > thresh_y) or (input_image[iy][ix] < Tmin) or (input_image[iy][ix] > Tmax):
                binary_M[iy][ix] = 1
            else:
                binary_M[iy][ix] = 0