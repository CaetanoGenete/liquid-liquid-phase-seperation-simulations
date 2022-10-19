from pde_tools.finite_difference import laplacian, gradient2D, divergence2D, pde_solve2D
from pde_tools.video import save_as_mp4
from pde_tools.other import min_max_normalise

import scipy.constants as constants

import numpy as np
random_seed = 69
np.random.seed(random_seed)

def model_b(phi, div_kernel, a, b, k, M, D, dx, dy, dt):
    height, width = phi.shape

    grad_mu_x, grad_mu_y = gradient2D(phi * (a + b * phi * phi) - k * laplacian(phi, div_kernel), div_kernel)

    noise_coeff = np.sqrt(2 * D/(dx*dy*dt))
    return divergence2D(
        M * grad_mu_x - noise_coeff * np.random.randn(height, width), 
        M * grad_mu_y - noise_coeff * np.random.randn(height, width),
        div_kernel)


#Model constants:
a = -1
b = 1
k = 1
M = 1 
# equal to K_b * T * M
D = 3e-2

#PDE solver constants
dx = 1
dy = 1
dt = 1e-3

#Time to solve up to:
t_max = 100

t, phi_t = pde_solve2D(
    model_b, 
    phi0  = np.random.randn(256, 256), 
    args  = [a, b, k, M, D, dx, dy, dt], 
    dx    = dx, 
    dy    = dy,
    dt    = dt,
    t_max = t_max,
    t_samples = 2500, 
    error_order = 8)

save_as_mp4(
    t, phi_t, 
    file_name = f"../Solutions/Model B/model_b(a={a:.2f} b={b:.2f} k={k:.2f} tmax={t_max}).mp4",
    xlabel    = f"$i\Delta x$ ($\Delta x=${dx})",
    ylabel    = f"$j\Delta y$ ($\Delta y=${dy})",
    title     = f"Model B with params: a = {a:.2f}, b = {b:.2f}, $\kappa$ = {k:.2f},\nM={M:.2f}",
    interval  = 20)