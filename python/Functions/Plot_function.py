#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:17:53 2025

@author: jleclezio
"""
import numpy as np
import matplotlib.pyplot as ppl


def plot_inTime_andSpace_linear(Final_Grid_L, params_plot, spacesL, SPECIES, Mb_pos, thickness):
    
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    indices_lin = params_plot["indices"]
    times_courbelin = params_plot["times_courbelin"]
    
    NS = len(SPECIES)
    
    Res = []

    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
        
        
    ppl.figure()
    sp_name=[]
    
    for space in range (0, len(spacesL)):
        
        to_plot = Res[space]
        # print(len(to_plot))
        x_values = real_spacesL[space]
        
        for species in range(0, NS):
            # pour t = 0 
            label = SPECIES[species].name if space == 0 else None
            
            ppl.plot(x_values, to_plot[0][species], color = 'black', marker=SPECIES[species].marker, markersize=4, label=label)
            
            for i, t in enumerate(times_courbelin):
                to_plot_log = [to_plot[i] for i in indices_lin]   
                
                ppl.plot(x_values, to_plot_log[i][species], color = colors[i], marker=SPECIES[species].marker, markersize=4)
    
    # griser les membranes : 
    for Mb in Mb_pos:
        pos_Mb = Mb[1]
        i = pos_Mb[0]
        
        space1 = real_Rspace[i]
        space2 = space1 + thickness
         
        ppl.axvspan(space1, space2, color='grey', alpha=0.5) 
        
    ppl.plot([], [], color='black', label='t=0')  
    for i,t in enumerate(times_courbelin):
        ppl.plot([], [], color=colors[i], label=f't={t:.4f}')
    
    ppl.title(f"Concentration des especes {sp_name} dans l'espace au cours du temps")
    ppl.xlabel("Position dans l'espace (x)")
    ppl.ylabel("Concentration")
    ppl.legend(loc='upper right')
    ppl.grid()
    
    ppl.show()

def plot_inTime_andSpace_log(Final_Grid_L, params_plot, spacesL, SPECIES, Mb_pos, thickness):
    
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    log_indices = params_plot["log_indices"]
    times_courbelog = params_plot["times_courbelog"]
    
    NS = len(SPECIES)
    
    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
        

    ppl.figure()
    
    for space in range (0, len(spacesL)):
        
        to_plot = Res[space]
        # print(len(to_plot))
        x_values = real_spacesL[space]
        
        for species in range(0, NS):
            # pour t = 0 
            label = SPECIES[species].name if space == 0 else None
            
            ppl.plot(x_values, to_plot[0][species], color = 'black', marker=SPECIES[species].marker, markersize=3, label=label)
    
            for i, t in enumerate(times_courbelog):
                to_plot_log = [to_plot[i] for i in log_indices]   
                # print(np.shape(to_plot_log))
                ppl.plot(x_values, to_plot_log[i][species], color = colors[i], marker=SPECIES[species].marker, markersize=3) # colors[i], marker=SPECIES[species].marker, markersize=4
    
    # griser les membranes : 
    for Mb in Mb_pos:
        pos_Mb = Mb[1]
        i = pos_Mb[0]
        
        space1 = real_Rspace[i]
        space2 = space1 + thickness
         
        ppl.axvspan(space1, space2, color='grey', alpha=0.5) 
        
    ppl.plot([], [], color='black', label='t=0')  
    for i,t in enumerate(times_courbelog):
        ppl.plot([], [], color=colors[i], label=f't={t:.5f}')
    
    ppl.title(f"Evolution of {SPECIES[species].name[0]} concentration throught an osmotic membrane, in space and time", fontsize=18)
    ppl.xlabel("Space position (x)", fontsize=16)
    ppl.ylabel("Concentration", fontsize=16)
    ppl.legend(loc='upper right')
    ppl.grid()
    
    # ppl.savefig("CSI/osmotic_membrane.pdf")
    
    ppl.show()
    
    
def plot_TimeandSpace_oneSpecies(Final_Grid_L, species_to_plot, params_plot, spacesL, SPECIES, Mb_pos, real_Rspace, thickness):
    
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    log_indices = params_plot["log_indices"]
    times_courbelog = params_plot["times_courbelog"]
    
    
    for i, species_l in enumerate(SPECIES):
        if species_l.name[0] == species_to_plot:
            pos_sp = i
    
    ppl.figure()
    
    sp_name=[]

    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
        

    for space in range (0, len(spacesL)):
        
        to_plot = Res[space]
        x_values = real_spacesL[space]
        
            # pour t = 0 
        label = SPECIES[pos_sp].name if space == 0 else None
        sp_name.append(SPECIES[pos_sp].name)if space == 0 else None
            
        ppl.plot(x_values, to_plot[0][pos_sp], color = 'black', marker=SPECIES[pos_sp].marker, markersize=4, label=label)
            
        for i, t in enumerate(times_courbelog):
            to_plot_log = [to_plot[i] for i in log_indices]   
                
            ppl.plot(x_values, to_plot_log[i][pos_sp], color = colors[i], marker=SPECIES[pos_sp].marker, markersize=4)
    # griser les membranes : 

    for Mb in Mb_pos:
        pos_Mb = Mb[1]
        i = pos_Mb[0]
            
        space1 = real_Rspace[i]
        space2 = space1 + thickness
         
        ppl.axvspan(space1, space2, color='grey', alpha=0.5) 
        
    ppl.plot([], [], color='black', label='t=0')  
    for i,t in enumerate(times_courbelog):
        ppl.plot([], [], color=colors[i], label=f't={t:.4f}')

    ppl.title(f"Concentration des especes {sp_name} dans l'espace au cours du temps")
    ppl.xlabel("Position dans l'espace (x)")
    ppl.ylabel("Concentration")
    ppl.legend(loc='upper right')
    ppl.grid()

    # savefig("electro_osmo_impl_Na_t1_Vm_p70.pdf")

    ppl.show()
    
def plot_t0(Final_Grid_L, spacesL, params_plot, SPECIES, Mb_pos, thickness):
    
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    log_indices = params_plot["log_indices"]
    times_courbelog = params_plot["times_courbelog"]
    
    NS = len(SPECIES)
    
    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
        
    
    ppl.figure()
    sp_name=[]

    for space in range (0, len(spacesL)):
        # type(space)
        # print(np.shape(Final_Grid_L))
        npFinal_Grid_L = np.array(Final_Grid_L)
        to_plot = npFinal_Grid_L[:, space, :, :]
        x_values = real_spacesL[space]
        
        for species in range(0, NS):
            # pour t = 0 
            
            label = SPECIES[species].name if space == 0 else None
            sp_name.append(SPECIES[species].name)if space == 0 else None
            
            ppl.plot(x_values, to_plot[0][species], color = colors[species], marker=SPECIES[species].marker, markersize=4, label=label)
     # griser les membranes : 

    for Mb in Mb_pos:
         pos_Mb = Mb[1]
         i = pos_Mb[0]
             
         space1 = real_Rspace[i]
         space2 = space1 + thickness
          
         ppl.axvspan(space1, space2, color='grey', alpha=0.5) 

    ppl.title(f"Concentration initiale (times = 0 ) de A en fonction de l'espace")
    ppl.xlabel("Position dans l'espace (x)", fontsize=14)
    ppl.ylabel("Concentration initale", fontsize=14)
    ppl.legend(loc='upper right')
    ppl.grid()

    # ppl.savefig("CI_10points.pdf")

    ppl.show()
    
def plot_tLast(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness):
    
     
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    log_indices = params_plot["log_indices"]
    times_courbelog = params_plot["times_courbelog"]
    
    NS = len(SPECIES)
    
    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
        
    ppl.figure()
    sp_name=[]

    for space in range (0, len(spacesL)):
        
        to_plot = Res[space]
        x_values = real_spacesL[space]
        
        for species in range(0, NS):
            # pour t = 0 
            label = SPECIES[species].name if space == 0 else None
            
            ppl.plot(x_values, to_plot[len(times)-1][species], color = 'black', marker=SPECIES[species].marker, markersize=4, label=label)
            print(to_plot[len(times)-1][species])
    # griser les membranes : 

    for Mb in Mb_pos:
        pos_Mb = Mb[1]
        i = pos_Mb[0]
            
        space1 = real_Rspace[i]
        space2 = space1 + thickness
         
        ppl.axvspan(space1, space2, color='grey', alpha=0.5) 

    ppl.title(f"Concentration finale (times = {times[len(times)-1]} ) de {sp_name} dans l'espace", fontsize=1)
    ppl.xlabel("Position dans l'espace (x)")
    ppl.ylabel("Concentration finale")
    ppl.legend(loc='upper right')
    ppl.grid()

    # ppl.savefig("elosm_NaK_t1_dtC09_D_np.pdf")

    ppl.show()
    
    
def plot_conservation(result_int, variations_relative, times, SPECIES):
    
    NS = len(SPECIES)
    
    # integral

    ppl.figure()
    for species in range(0, NS):
        ppl.plot(times, [array[species] for array in result_int], label=SPECIES[species].name)
    ppl.xlabel("Temps")
    ppl.ylabel("Integrale sur l'espace")
    ppl.title("Vérification de la conservation")
    ppl.grid()
    ppl.legend()


    # variations relatives
    ppl.figure()
    for species in range(0,NS): 
        ppl.plot(times, [array[species] for array in variations_relative], label=SPECIES[species].name)
    ppl.xlabel("Temps")
    ppl.ylabel("Variations relative")
    ppl.title("Vérification de la conservation")
    ppl.legend()
    ppl.grid()

    # ppl.savefig("conservation_IS_small_middle_t1_diff_bmailles_Ptrue.pdf")

    ppl.show()
    

def plot_3D(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness ):
     
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    times_courbelog = params_plot["times_courbelog"]
    
    NS = len(SPECIES)
    
    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
    

    #def nombre de courbe
    nb_courbe_3D = 10
    log_indices_3D = np.logspace(0, np.log10(len(times)-1), nb_courbe_3D, dtype=int) 
    times_courbelog_3D = np.array(times)[log_indices_3D.astype(int)]

    #def couleur des courbes
    cmap = ppl.get_cmap('viridis')  
    colors = [cmap(i / (nb_courbe_3D - 1)) for i in range(nb_courbe_3D)]


    fig = ppl.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    # ax.set_xscale('log')

    for space in range(len(spacesL)):
        to_plot = Res[space]         
        x_values = real_spacesL[space] 

        for species in range(NS):
            for j, t_idx in enumerate(log_indices_3D):
                t_real = times_courbelog_3D[j]            
                y_t = to_plot[t_idx][species]
                x_t = x_values
                z_t = np.full_like(x_t, np.log(t_real))

                ax.plot(z_t, x_values, y_t,
                        color=colors[j],
                        marker=SPECIES[species].marker,
                        markersize=2)  
                

    for Mb in Mb_pos:
        i = Mb[1][0]
        x1 = real_Rspace[i]
        x2 = x1 + thickness

        for t_val in times_courbelog:
            ax.plot([t_val, t_val], [x1, x2], [0, 0], color='grey', alpha=0.3)

    ax.set_xlabel('Time (t)', fontsize=14)
    ax.set_ylabel('Space (x)', fontsize=14)
    ax.set_zlabel('Concentration', fontsize=14)

    ax.set_xticks([])
    ax.set_xticklabels([])


    ax.set_title(f"Evolution of {SPECIES[0].name[0]} concentration in space and time", fontsize=16)

    for j, t in enumerate(times_courbelog_3D):
        ax.plot([], [], [], color=colors[j], label=f't = {t:.1e}')

    ax.legend(loc='upper right', bbox_to_anchor=(0.85, 1.0), fontsize=12)

    # ppl.savefig("CSI/diff_simple_timlog.pdf")

    ppl.show()
    
    
def plot_heatmap(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness ):
     
    colors = params_plot["colors"]
    real_spacesL = params_plot["real_spacesL"]
    real_Rspace = params_plot["real_Rspace"]
    times_courbelog = params_plot["times_courbelog"]
    
    NS = len(SPECIES)
    
    Res = []
    for i in range(0, len(Final_Grid_L[0])):
        modif_Final_Grid_L = []
        for i1 in range(0, len(Final_Grid_L)):
            modif_Final_Grid_L.append(Final_Grid_L[i1][i])
        Res.append(modif_Final_Grid_L)
    
    
    # from matplotlib.colors import LogNorm  # pour l’échelle logarithmique
    
    nb_courbe_hm = 100
    
    raw_log_indices = np.logspace(0, np.log10(len(times) - 1), nb_courbe_hm * 10)
    log_indices = np.unique(np.round(raw_log_indices).astype(int))
    log_indices = log_indices[log_indices < len(times)]
    times_courbelog = np.array(times)[log_indices.astype(int)]
    
    
    # Choix du sous-ensemble : un seul "space"
    space = 0
    species = 0
    
    # Extraire les données
    to_plot = Res[space]  # shape: (Nt, NS, Nx)
    x_values = real_spacesL[space]
    Nt = len(log_indices)
    Nx = len(x_values)
    
    # Construction de la matrice de concentration [temps, espace]
    data = np.zeros((len(log_indices), len(x_values)))
    time_labels = [times_courbelog[j] for j in range(len(log_indices))]
    
    
    for j, t_idx in enumerate(log_indices):
        data[j, :] = to_plot[t_idx][species]
        
    fig, ax = ppl.subplots(figsize=(8, 5))
    im = ax.imshow(data, aspect='auto', origin='lower',
                   extent=[x_values[0], x_values[-1],
                           time_labels[0], time_labels[-1]],
                   cmap='viridis')
    
    
    ax.set_yscale('log')
    ax.set_xlabel("Espace (x)")
    ax.set_ylabel("Temps (t)")
    fig.colorbar(im, ax=ax, label="Concentration")
    ppl.tight_layout()
    
    # ppl.savefig("CSI/heatmap_diff_simple_timelog.pdf")
    
    ppl.show()
        
