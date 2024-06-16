import numpy as np
# D = min [ alpha * D sub 0.99, Dmax]
"""
according to diagram shown on document, i will be using the following estimated values:
 # 99th percentile for x gradient: 0.38
 # 99th percentile for y gradient: 0.35
 # max threshold for x gradient: 0.45
 # max threshold for y gradient: 0.42
 # alpha = 1.2 (can be changed)
"""
alpha = 1.2
ninX = 0.38
ninY = 0.35
maxGrX = 0.45
maxGrY = 0.42
calcDx = alpha * ninX
Dx = min (calcDx, maxGrX)

calcDy = alpha * ninY
Dy = min (calcDy, maxGrY)

print(Dx)
print(Dy)

"""
M Matrix (mask)

"""
# Laplacian based on restricted gradients
# Algebraic Eqs for Different Fourier Components
# Fourier Transform Algo from C --> Python
# Remove artifacts from striped image component - using nonlinear Gaussian Type Filter
"""
Filter parameter info
Calculate...
- delta r
- sigma0, sigma0^2
- sigma
- value of filter parameter
"""