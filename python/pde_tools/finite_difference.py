import numpy as np
from scipy.signal import convolve2d

def _vandermonde_row2_denominator(coeffs, j):
    """
    Computes the denominator of the entry at the second row and jth column of the inverse
    vandermonde matrix, which is given by: product of all pairs (1<=l<k<=N) (a_k - a_l) with l=j, which may be written
    as the product from  k=0->j-1: (a_j - a_k) and k=j+1->N (a_k - a_j).

    Args:
        coeffs (np.array): The coefficients of the polynomial related to the vandermonde
        matrix.
        j (int): The column of the inverse vandermonde matrix

    Returns:
        Integer: The denominator of the entry at the second row and jth column of the inverse
        vandermonde matrix
    """

    multiple_j = coeffs[j]

    denominator = 1
    for k, sample_k in enumerate(coeffs):
        if(j != k):
            denominator *= (multiple_j - sample_k) * np.sign(j-k)

    return denominator


def _vandermonde_row2_numerator(coeffs, j, product):
    """
    Computes the numerator of the entry at the second row and jth column of the inverse
    vandermonde matrix, which is given by the sum of the ordered products of the coeffients,
    which does not contain the jth coeffcient and one other. 

    Args:
        coeffs (np.array): The coefficients of the polynomial related to the vandermonde
        matrix.
        j (int): The column of the inverse vandermonde matrix
        product (_type_): Optimisation. Precomputed product of the coefficients.

    Returns:
        Integer: The numerator of the entry at the second row and jth column of the inverse
        vandermonde matrix
    """

    partial_numerator = 0

    for k, sample_k in enumerate(coeffs):
        if k != j:
            partial_numerator += product // sample_k

    return partial_numerator


def inv_vandermonde_row2(coeffs):
    """
    O(n^2) algorithm for computing the second row of the inverse vandermonde matrix.

    Note: currently very slow due to python nonsense, rewrite in C++.

    Args:
        coeffs (np.array): The coefficients of the polynomial related to the vandermonde
        matrix.
    Returns:
        The second row of the inverse vandermonde matrix
    """
    
    #Vandermonde matrix has no inverse if coefficients aren't unique
    unique_offsets = np.unique(coeffs)
    zeros = unique_offsets == 0

    N, = unique_offsets.shape#
    assert N > 1, "Atleast two sample points are required!"

    result = np.zeros(N)

    #Entry sign coefficient
    parity = -1
    #Special care must be taken if one of the coefficients is zero
    if(zeros.any()):
        non_zero_prod = np.prod(unique_offsets, where=np.logical_not(zeros))

        for j, multiple_j in enumerate(unique_offsets):
            numerator = 0
            #All other elements except at the jth index are non-zero
            if multiple_j == 0:
                numerator = _vandermonde_row2_numerator(unique_offsets, j, non_zero_prod)
            #Product will always contain a zero, except when the kth term is zero, which occurs only once
            else:
                numerator = non_zero_prod // multiple_j 

            result[j] = parity * (numerator/_vandermonde_row2_denominator(unique_offsets, j))
            #Flip the parity
            parity *= -1
        
    else:
        product = np.prod(unique_offsets)

        for j, multiple_j in enumerate(unique_offsets):
            numerator = _vandermonde_row2_numerator(unique_offsets, j, product) // multiple_j

            result[j] = parity * (numerator/_vandermonde_row2_denominator(unique_offsets, j))
            #Flip the parity
            parity *= -1

    return result


def finite_difference_stencil(grid_offsets, order, dx=1):
    """
    Calculates the stencil for the finite difference approximation of the
    'order'th derivative at the point x, in the base given by:
    [f(x + grid_offsets * dx)].

    Args:
        grid_offsets (np.array): If a_i is the i^th entry of this variable, then
        the finite difference approximation will be calculated in the base given
        by: [f(x + a_1 * dx), f(x + a_2 * dx), ...].
        order (int): The order of the derivative.
        dx (float, optional): The spacing between grid points. Defaults to 1.

    Returns:
        np.array: The stencil for the finite difference approximation of the
        'order'^th derivative at the point x
    """

    #Vandermonde matrix is only invertible if coefficients are unique,
    #here the grid_offsets will be used as the coefficients of the 
    #polynomial. 
    unique_offsets = np.unique(grid_offsets)

    N, = unique_offsets.shape
    assert N > order, "Atleast 'order + 1' number of sample points are required!"

    mat = np.zeros((N, N))

    #Construct vandermonde matrix
    column = np.ones(N)
    for power in range(N):
        mat[:, power] = column
        column *= unique_offsets

    return (np.linalg.inv(mat)[order] / np.power(dx, order)) * np.math.factorial(order)


def gradient2D(phi, div_kernel):
    """
    Computes the gradient of the vector 2D vector field.

    Args:
        phi (np.array): A 2D vector field
        div_kernel (np.array): A 2D array, where the real part and imaginary parts
        represents differentiation with respect to x and y, respectively.

    Returns:
        np.array, np.array: The gradient of the vector field
    """
    result = convolve2d(phi, div_kernel, mode="same", boundary="wrap")

    return result.real, result.imag


def divergence2D(phi_x, phi_y, div_kernel):
    m, n = div_kernel.shape
    
    return convolve2d(phi_x, div_kernel[m//2, :].real.reshape(1, -1), mode="same", boundary="wrap") +\
           convolve2d(phi_y, div_kernel[:, n//2].imag.reshape(-1, 1), mode="same", boundary="wrap")


def laplacian(phi, div_kernel):
    return divergence2D(*gradient2D(phi, div_kernel), div_kernel)


def calculate_central_div_kernel(order, dx, dy):
    assert order % 2 == 0, "Central difference should have even order"

    half_order = order//2
    #Generate order 8 stencil for partial derivative
    stencil = finite_difference_stencil(np.arange(-half_order, half_order+1), 1)

    div_kernel = np.zeros((order+1, order+1), dtype=np.complex128)
    #Real and imaginary part of kernel represent differentiation with respect to x and y, respectively
    div_kernel[half_order, :] += stencil / dx
    div_kernel[:, half_order] += stencil * 1j / dy

    return div_kernel


def pde_solve2D(time_derivative, phi0, dx=1, dy=1, dt=1e-3, args = [], t_max=100, t_samples=1500, error_order=8):
    div_kernel = calculate_central_div_kernel(error_order, dx, dy)

    #The t-interval between stored values of phi
    sample_rate = int(t_max / (dt * t_samples))
    #The number of steps in the finite difference method
    step_count = t_samples * sample_rate

    #Start with random noise
    phi   = phi0
    phi_t = np.zeros((t_samples, *phi0.shape))
    t     = np.arange(t_samples) * dt * sample_rate

    print("[INFO] Integration started:")
    print("progress: 0%", end="\r")

    for step in range(step_count):
        if not step % sample_rate:
            phi_t[step // sample_rate] = phi

        if not step % 100:
            print(f"progress: {step * 100 / step_count:.3f}%", end="\r")

        phi += dt * time_derivative(phi, div_kernel, *args)

    print("progress: 100.000% DONE!")
    return t, phi_t


