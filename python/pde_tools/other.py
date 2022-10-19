import numpy as np


def min_max_normalise(phi, new_min=0, new_max=1):
    phi_min = np.min(phi)
    phi_max = np.max(phi)

    return (phi - phi_min) * (new_max - new_min) / (phi_max - phi_min) + new_min