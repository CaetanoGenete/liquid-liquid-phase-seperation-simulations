from pde_tools.finite_difference import laplacian, pde_solve2D
from pde_tools.video import save_as_mp4

import numpy as np
random_seed = 69
np.random.seed(random_seed)


def noiseless_model_b(phi, div_kernel, a, b, k):
    return laplacian(phi * (a + b * phi * phi) - k * laplacian(phi, div_kernel), div_kernel)


#Model constants:
a=-1
b=1

#PDE solver constants
dx = 1
dy = 1

#Time to solve up to:
t_max = 100

t, phi_t = pde_solve2D(
    noiseless_model_b, 
    phi0  = np.random.randn(256, 256), 
    args  = [a, b, k], 
    dx    = dx, 
    dy    = dy,
    t_max = t_max,
    t_samples = 2500, 
    error_order = 6)

save_as_mp4(
    t, phi_t, 
    file_name=f"../Solutions/Noiseless model B/noiseless_model_b(a={a:.2f} b={b:.2f} k={k:.2f} tmax={t_max}).mp4",
    xlabel   =f"$i\Delta x$ ($\Delta x=${dx})",
    ylabel   =f"$j\Delta y$ ($\Delta y=${dy})",
    title    =f"Noiseless Model B with params: a = {a:.2f}, b = {b:.2f}, $\kappa$ = {k:.2f}",
    interval=20)