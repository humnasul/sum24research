import numpy as np
import pyfftw

def compute_poisson_rec(input_image, output_image, binary_M, grid_dct, sx, sy, Lap, Lap_dct, plan_forward, plan_backward):
    for iy in range(sy):
        for ix in range(sx):
            Lap[iy, ix] = 0.0

    # Calculate Laplacian
    for iy in range(1, sy - 1):
        for ix in range(1, sx - 1):
            Lap[iy, ix] = input_image[iy, ix - 1] - 2.0 * input_image[iy, ix] + input_image[iy, ix + 1]
            if binary_M[iy, ix] != 0:
                Lap[iy, ix] += input_image[iy + 1, ix] - input_image[iy, ix]
            if binary_M[iy - 1, ix] != 0:
                Lap[iy, ix] += input_image[iy - 1, ix] - input_image[iy, ix]

    # Solve Discrete Poisson Equation
    Lap_dct[0, :] = pyfftw.interfaces.numpy_fft.rfft(Lap[0, :])
    Lap_dct[0, :] *= grid_dct[0, :]
    output_image[0, :] = pyfftw.interfaces.numpy_fft.irfft(Lap_dct[0, :])

    return

# Example input data (replace with your actual data)
sx = 100  # Image width
sy = 80   # Image height
input_image = np.random.rand(sy, sx)
binary_M = np.random.randint(0, 2, size=(sy, sx))
grid_dct = np.random.rand(sy, sx)

# Initialize output arrays
output_image = np.zeros_like(input_image)
Lap = np.zeros_like(input_image)
Lap_dct = np.zeros_like(input_image)

# Lap[0] = pyfftw.empty_aligned(128, dtype='complex128')
# Lap_dtc[0] = pyfftw.empty_aligned(128, dtype='complex128')

# Create FFTW plans
plan_forward = pyfftw.FFTW(Lap[0], Lap_dct[0], direction='FFTW_FORWARD', flags=('FFTW_MEASURE',))
plan_backward = pyfftw.FFTW(Lap_dct[0], Lap[0], direction='FFTW_BACKWARD', flags=('FFTW_MEASURE',))

# Call the compute_poisson_rec function
compute_poisson_rec(input_image, output_image, binary_M, grid_dct, sx, sy, Lap, Lap_dct, plan_forward, plan_backward)

# Print or use the resulting output_image
print(output_image)
