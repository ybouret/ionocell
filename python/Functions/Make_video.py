#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 16:11:14 2025

@author: jleclezio
"""

import numpy as np
import matplotlib.pyplot as ppl
from moviepy import ImageSequenceClip



def generate_images_new(Final_Grid_L, times, spacesL, Respace, Mb_dict, SPECIES, thickness, c_log, nb_frames=50):
    
    if c_log == None: 
        log_indices = np.logspace(0, np.log10(len(times)-1), nb_frames, dtype=int) 
        times_courbelog = np.array(times)[log_indices.astype(int)]
        
    elif c_log == False:
        
        indices = np.linspace(0, len(times) - 1, nb_frames, dtype=int)
        times_courbelog = np.array(times)[indices]
        
        
        log_indices = indices

    else : 
        mid_idx = np.searchsorted(times, c_log)
        
        log_indices_1 = np.logspace(0, np.log10(mid_idx-1), int(nb_frames/2), dtype=int)
        log_indices_2 = np.logspace(0, np.log10(len(times) - mid_idx - 1), int(nb_frames/2), dtype=int) + mid_idx
        
        # log_indices_1 = np.logspace(0, np.log10((len(times)/2)-1), int(nb_frames/2), dtype=int)
        # log_indices_2 = np.logspace(0, np.log10((len(times)/2)-1), int(nb_frames/2), dtype=int) + round(len(times)/2)
        
        log_indices = np.concatenate([log_indices_1, log_indices_2])
        times_courbelog = np.array(times)[log_indices]
        
    


    # Préparer la position réelle des compartiments dans l'espace
    real_spacesL = []
    real_Rspace = []
    for i in range(len(spacesL)):
        offset = sum(Respace[:i])
        real_spacesL.append(spacesL[i] + offset)
        real_Rspace.append(Respace[i] + offset)

    # Couleurs
    cmap = ppl.colormaps['Dark2']  
    discrete_cmap = cmap.resampled(nb_frames) 
    colors = [discrete_cmap(i) for i in range(nb_frames)]

    images = []

    for idx, t_index in enumerate(log_indices):
        fig, ax = ppl.subplots()
        
        for space in range(len(spacesL)):
            
            to_plot = Final_Grid_L[t_index][space]  # [species, points]
            x_values = real_spacesL[space]
            
            for species in range(len(SPECIES)):
                label = SPECIES[species].name if space == 0 else None
                ax.plot(x_values, to_plot[species], color=colors[idx], marker=SPECIES[species].marker,
                        markersize=4, label=label)

        # Membranes
        for Mb in Mb_dict:
            i = Mb_dict[Mb][0]
            x1 = real_Rspace[i]
            x2 = x1 + thickness
            ax.axvspan(x1, x2, color='grey', alpha=0.5)

        ax.set_title(f"t = {times_courbelog[idx]:.4f}")
        ax.set_xlabel("Position dans l'espace (x)")
        ax.set_ylabel("Concentration")
        ax.set_ylim(0, 1)
        ax.grid()
        ax.legend(loc='upper right')

        fig.canvas.draw()
        image = np.array(fig.canvas.renderer.buffer_rgba())
        images.append(image)
        ppl.close(fig)

    return images


def create_video_from_images(images, output_file, fps=8):
    clip = ImageSequenceClip(images, fps=fps)
    clip.write_videofile(output_file, codec='libx264')
