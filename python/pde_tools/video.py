import matplotlib.pyplot as plt
import matplotlib.animation as anim

import numpy as np

import os

def save_as_mp4(t, phi_t, file_name, xlabel = "", ylabel = "", title = "", cmap="gnuplot2", **kwargs):
    fig, ax = plt.subplots()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    vmin = np.min(phi_t)
    vmax = np.max(phi_t)

    frames = []
    for j, (time, frame) in enumerate(zip(t, phi_t)):
        image = plt.imshow(
            frame, 
            animated = j!=0, 
            vmin = vmin, vmax = vmax, 
            cmap = cmap)

        time_text = plt.text(
            0.02, 0.98, 
            f"time = {time:.3f}", 
            color = "white", 
            va = "top", ha = "left",
            transform = ax.transAxes)

        frames.append([image, time_text])

    animation = anim.ArtistAnimation(fig, frames, **kwargs)

    cbar = plt.colorbar()
    cbar.set_label('$\phi(t, x, y)$', rotation=270)

    curr_directory = os.path.dirname(__file__)
    full_file_path = os.path.normpath(os.path.join(curr_directory , file_name))

    print("[INFO] Saving solution to:", str(full_file_path))
    animation.save(str(full_file_path), bitrate=-1)
    print("[INFO] saving completed!")