# -*- coding: utf-8 -*-
"""
Created on Tue Dec 2 11:38:17 2024

@author: jeann

Resolution de l'equation de diffusion pour une concentration 
et la methode des elemenst finis'
"""

import numpy as np
import matplotlib.pyplot as ppl

import ionocell

import sys
sys.path.insert(0, '/Users/jleclezio/Documents/ionocell/python/Functions/')

from Initiation_functions import calcul
from Make_video import generate_images_new, create_video_from_images


ppl.rcParams['figure.figsize'] = [13, 6]



""" Definition des données d'entrées """

'''Def membrane'''

# nombre d'espace pour la simulation 
space_nb = 3

Membrane_list = ["osmotique implicite", "osmotique analytique", "electro-osmotique implicite", "non permeable"]

# consequences: 
allMb_dict = {
    "Mb1" : [0,1],
    "Mb2" : [1,2], 
    "Type" : Membrane_list[1]}
    # {"type" : Membrane_list[1]}
    # }

# type de membrane
# Membrane_list = ["osmotique implicite", "osmotique analytique", "electro-osmotique implicite", "non permeable"]
# Type_mb = Membrane_list[1]

# Params_calcul = [Type_mb, space_nb]


'''parametres modifiables'''

# 2**5, round(100*(4/9))
nbmailles_x1 = round(100*(4/9))
nbmailles_x2 = round(100*(1/9))
nbmailles_x3 = round(100*(4/9))

allmailles_x = [nbmailles_x1, nbmailles_x2, nbmailles_x3]
# nbmailles_x = [nbmailles_x1, nbmailles_x2]
# nbmailles_x = [nbmailles_x1]

# ratio x1 vs x2: taille dans la cuve
Respace1 = 4/9
Respace2 = 1/9
Respace3 = 4/9
all_Respace = [Respace1, Respace2, Respace3]
# Respace = [Respace1, Respace2]
# Respace = [Respace1]

# TEMPS: de 0 à 1
ti = 0
tf = 1

# CONSERVATION : seuil de tolérance
tolerance = 0.05

# CI GAUSSIENNE: especes decritent à t=0 par gaussienne
var_g = 0.005 # largeur de la gaussienne


'''Calcul pas de temps et d'espaces'''

# PAS D'ESPACE
delta_x1 = Respace1/(nbmailles_x1-1)
delta_x2 = Respace2/(nbmailles_x2-1)
delta_x3 = Respace3/(nbmailles_x3-1)
alldelta_x = [delta_x1, delta_x2, delta_x3]
# delta_x = [delta_x1, delta_x2]
# delta_x = [delta_x1]

# espace
x1 = np.arange(nbmailles_x1)* delta_x1
x2 = np.arange(nbmailles_x2)* delta_x2 
x3 = np.arange(nbmailles_x3)* delta_x3
allspacesL = [x1, x2, x3]
# spacesL = [x1, x2]
# spacesL = [x1]

# PAS DE TEMPS A MODIFIER : AJOUTER D POUR CHAQUE ESPECE
delta_t1 = 0.9 * (delta_x1**2) / 2 # 0.9 * delta_t max 
delta_t2 = 0.9 * (delta_x2**2) / 2 # 0.9 * delta_t max 
delta_t3 = 0.9 * (delta_x3**2) / 2 # 0.9 * delta_t max 
alldelta_t = [delta_t1,delta_t2, delta_t3] # recuperer le delta_t min = celui avec le plus de precisions
# delta_t = min(delta_t1,delta_t2)
# delta_t = delta_t1

thickness = min(delta_x1, delta_x2, delta_x3)/2 #largeur de la membrane pour le plot
# thickness = delta_x1/2 #largeur de la membrane pour le plot

x1_info = ionocell.Space(_number = 0, _type = "Cell", _two_wall = False, _one_inje = True ) 
x2_info = ionocell.Space(_number = 1, _type = "Extra", _two_wall = True, _one_inje = False ) 
x3_info =  ionocell.Space(_number = 2, _type = "Cell", _two_wall = True, _one_inje = False ) 


allInfo_space = [x1_info, x2_info, x3_info]

'''nombre d'espace'''

spacesL = allspacesL[:space_nb]
nbmailles_x = allmailles_x[:space_nb]
delta_t = min(alldelta_t[:space_nb])
delta_x = alldelta_x[:space_nb]
Info_space = allInfo_space[:space_nb]
Respace = all_Respace[:space_nb]

Mb_pos = list(allMb_dict.items())[:space_nb-1]
type_mb = list(allMb_dict.items())[-1]

Mb_dict = dict(Mb_pos + [type_mb])

'''description d'especes'''

# test avec Class

Na = ionocell.Specie(_specie_name = "Na", _diff_coeff = 1, _perm_coeff = 1, _charge = +1, _marker= 'o')
Cl = ionocell.Specie(_specie_name = "Cl", _diff_coeff = 1, _perm_coeff = 1, _charge = -1, _marker= 's')
NaCl = ionocell.Specie(_specie_name = "NaCl", _diff_coeff = 1, _perm_coeff = 10, _charge = 0, _marker= 'd')
OH = ionocell.Specie(_specie_name = "OH", _diff_coeff = 1, _perm_coeff = 1, _charge = -1, _marker= 'd')
H = ionocell.Specie(_specie_name = "H", _diff_coeff = 1, _perm_coeff = 1, _charge = +1, _marker= 's')
H20 = ionocell.Specie(_specie_name = "H2O", _diff_coeff = 1, _perm_coeff = 1, _charge = 0, _marker= 'v')

# liste des especes utilisées poru la simulation (decrites: "Na", "Cl", "OH", "NaCl")
# si pas decrite: rajouter les caracteristiques dans les dict + les CI
# SEULE LISTE A MODIFIER POUR AJOUTER/ENLEVER DES ESPECES POUR UNE SIMULATION
SPECIES = [Cl] 

Marker = ["o", "d", "s", "h", "v"]

'''reaction'''

# TFreaction = True

# d_react = {
#     "reaction0": {"reactif" : ["Na", "Cl"] , \
#                   "produit" : ["NaCl"],  \
#                   "cst_eq" : 1} }

TFreaction = False
d_react = {}




'''Condition initiale des elements'''

# longueur de la cuve: arbitraire
L1 = Respace1
L2 = Respace2

Gauss_1_4_x1 = np.exp(-((x1-(L1*1/4))**2) / var_g) 
Gauss_1_4_x2 = np.exp(-((x2-(L2*1/4))**2) / var_g)

Gauss_3_4_x1 = np.exp(-((x1-(L1*3/4))**2) / var_g) 
Gauss_3_4_x2 = np.exp(-((x2-(L2*3/4))**2) / var_g)

cst_1_x1 = np.ones(nbmailles_x1)
cst_1_x2 = np.ones(nbmailles_x2)

Null_x1 = np.zeros(nbmailles_x1) 
Null_x2 = np.zeros(nbmailles_x2)
Null_x3 = np.zeros(nbmailles_x3)

IC_dict = {
    "Na_x1": Null_x1, 
    "Na_x2": Null_x2, 
    "Na_x3": Null_x3, 
    
    "Cl_x1": Null_x1, 
    "Cl_x2": Null_x2,
    "Cl_x3": Null_x3,
    
    "NaCl_x1": Null_x1, 
    "NaCl_x2": Null_x2,
    "OH_x1": Null_x1, 
    "OH_x2": Gauss_3_4_x2, 
    "H_x1": Null_x1, 
    "H_x2": Gauss_3_4_x2, 
    "H2O_x1": Null_x1, 
    "H2O_x2": Null_x2
    }

'''Params model pour calcul'''

# # Calcul_gradient_avant = ["reaction", "diffusion", "transport"]
# # Gradient = Calcul_gradient_avant[2]

# Membrane_list = ["osmotique implicite", "osmotique analytique", "electro-osmotique implicite", "non permeable"]
# Type_mb = Membrane_list[1]

# Params_calcul = [Type_mb]


# # Mb_dict = {"Mb1" : [0,1]}

# Mb_dict = {
#     "Mb1" : [0,1],
#     "Mb2" : [1,2]}


""" Calculs """

(Final_Grid_L, integrals_grid_L, times) = calcul(SPECIES, spacesL, IC_dict, d_react, nbmailles_x, ti, tf, delta_t, delta_x, TFreaction, Mb_dict, Info_space)


integrals_grid_L = np.array(integrals_grid_L)

NS = len(SPECIES)


""" Analyse de la conservation """  

result_int = np.zeros(integrals_grid_L[:, 0, :].shape) 

for i in range(integrals_grid_L.shape[1]): 
    result_int += integrals_grid_L[:, i, :]

init_int = result_int[0]

variations_relative = [abs((s - init_int) / init_int) for s in result_int]

for species in range(0, NS):    
    variations_relative_specie = [array[species] for array in variations_relative]
    non_conservative = any(var > tolerance for var in variations_relative_specie)
    
    if non_conservative:
        print("Le modèle avec", SPECIES[species].name, "n'est pas conservatif.")
    else:
        print("Le modèle", SPECIES[species].name, "est conservatif.")



""" Plot resultats """

#pas de temps:
nb_courbe_t = 8    # PEUT ETRE MODIFIÉ POUR PLOT PLUS DE COURBE = PLUS DE TEMPS


# description des couleurs en fct du nb de courbe
cmap = ppl.colormaps['Dark2']  
discrete_cmap = cmap.resampled(nb_courbe_t) 
colors = [discrete_cmap(i) for i in range(nb_courbe_t)]

'''parametres plot'''

#recuperation des elements espacés en log: sur tout le temps
log_indices = np.logspace(0, np.log10(len(times)-1), nb_courbe_t, dtype=int) 
times_courbelog = np.array(times)[log_indices.astype(int)]


#recuperation des elements espacés en log: sur la moitié du temps
log_indices_1 = np.logspace(0, np.log10((len(times)/2)-1), int(nb_courbe_t/2), dtype=int)
log_indices_2 = np.logspace(0, np.log10((len(times)/2)-1), int(nb_courbe_t/2), dtype=int) + round(len(times)/2)

log_indices_2log = np.concatenate([log_indices_1, log_indices_2])
times_courbelog_2log = np.array(times)[log_indices_2log]


real_spacesL = []
real_Rspace = []
for i in range(len(spacesL)):
    offset = sum(Respace[:i])
    real_spacesL.append(spacesL[i] + offset)
    real_Rspace.append(Respace[i] + offset)

'''Plot final'''


Res = []
for i in range(0, len(Final_Grid_L[0])):
    modif_Final_Grid_L = []
    for i1 in range(0, len(Final_Grid_L)):
        modif_Final_Grid_L.append(Final_Grid_L[i1][i])
    Res.append(modif_Final_Grid_L)
    
# 1 log
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
        
        for i, t in enumerate(times_courbelog):
            to_plot_log = [to_plot[i] for i in log_indices]   
            
            ppl.plot(x_values, to_plot_log[i][species], color = colors[i], marker=SPECIES[species].marker, markersize=4)

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


# # 2 log
# ppl.figure()
# sp_name=[]

# for space in range (0, len(spacesL)):
    
#     to_plot = Final_Grid_L[:, space, :, :]
#     x_values = real_spacesL[space]
    
#     for species in range(0, NS):
#         # pour t = 0 
#         label = SPECIES[species].name if space == 0 else None
        
#         ppl.plot(x_values, to_plot[0][species], color = 'black', marker=SPECIES[species].marker, markersize=4, label=label)
        
#         for i, t in enumerate(times_courbelog_2log):
#             to_plot_log = [to_plot[i] for i in log_indices_2log]   
            
#             ppl.plot(x_values, to_plot_log[i][species], color = colors[i], marker=SPECIES[species].marker, markersize=4)
# # griser les membranes : 

# for Mb in Mb_pos:
#     pos_Mb = Mb[1]
#     i = pos_Mb[0]
        
#     space1 = real_Rspace[i]
#     space2 = space1 + thickness
     
#     ppl.axvspan(space1, space2, color='grey', alpha=0.5) 
    
# ppl.plot([], [], color='black', label='t=0')  
# for i,t in enumerate(times_courbelog):
#     ppl.plot([], [], color=colors[i], label=f't={t:.4f}')

# ppl.title(f"Concentration des especes {sp_name} dans l'espace au cours du temps")
# ppl.xlabel("Position dans l'espace (x)")
# ppl.ylabel("Concentration")
# ppl.legend(loc='upper right')
# ppl.grid()


# # savefig("electro_osmo_impl_Na_t1_Vm_p70.pdf")

# ppl.show()




'''plot seulement une espece (definir species=)'''

# ppl.figure()
# sp_name=[]

# for space in range (0, len(spacesL)):
    
#     to_plot = Final_Grid_L[:, space, :, :]
#     x_values = real_spacesL[space]
    
#     species = 1
#         # pour t = 0 
#     label = SPECIES[species].name if space == 0 else None
#     sp_name.append(SPECIES[species].name)if space == 0 else None
        
#     ppl.plot(x_values, to_plot[0][species], color = 'black', marker=SPECIES[species].marker, markersize=4, label=label)
        
#     for i, t in enumerate(times_courbelog):
#         to_plot_log = [to_plot[i] for i in log_indices]   
            
#         ppl.plot(x_values, to_plot_log[i][species], color = colors[i], marker=SPECIES[species].marker, markersize=4)
# # griser les membranes : 

# for Mb in Mb_pos:
#     pos_Mb = Mb[1]]
#     i = pos_Mb[0]
        
#     space1 = real_Rspace[i]
#     space2 = space1 + thickness
     
#     ppl.axvspan(space1, space2, color='grey', alpha=0.5) 
    
# ppl.plot([], [], color='black', label='t=0')  
# for i,t in enumerate(times_courbelog):
#     ppl.plot([], [], color=colors[i], label=f't={t:.4f}')

# ppl.title(f"Concentration des especes {sp_name} dans l'espace au cours du temps")
# ppl.xlabel("Position dans l'espace (x)")
# ppl.ylabel("Concentration")
# ppl.legend(loc='upper right')
# ppl.grid()


# # savefig("electro_osmo_impl_Na_t1_Vm_p70.pdf")

# ppl.show()

'''plot au premier temps (CI)'''


# ppl.figure()
# sp_name=[]

# for space in range (0, len(spacesL)):
#     to_plot = Final_Grid_L[:, space, :, :]
#     x_values = real_spacesL[space]
    
#     for species in range(0, NS):
#         # pour t = 0 
        
#         label = SPECIES[species].name if space == 0 else None
#         sp_name.append(SPECIES[species].name)if space == 0 else None
        
#         ppl.plot(x_values, to_plot[0][species], color = colors[species], marker=SPECIES[species].marker, markersize=4, label=label)
#  # griser les membranes : 

# for Mb in Mb_pos:
#      pos_Mb = Mb[1]
#      i = pos_Mb[0]
         
#      space1 = real_Rspace[i]
#      space2 = space1 + thickness
      
#      ppl.axvspan(space1, space2, color='grey', alpha=0.5) 

# ppl.title(f"Concentration initiale (times = 0 ) de A en fonction de l'espace")
# ppl.xlabel("Position dans l'espace (x)")
# ppl.ylabel("Concentration initale")
# ppl.legend(loc='upper right')
# ppl.grid()


# ppl.show()


'''Plot au dernier temps'''


# Res = []
# for i in range(0, len(Final_Grid_L[0])):
#     modif_Final_Grid_L = []
#     for i1 in range(0, len(Final_Grid_L)):
#         modif_Final_Grid_L.append(Final_Grid_L[i1][i])
#     Res.append(modif_Final_Grid_L)


ppl.figure()
sp_name=[]

for space in range (0, len(spacesL)):
    
    to_plot = Res[space]
    x_values = real_spacesL[space]
    
    for species in range(0, NS):
        # pour t = 0 
        label = SPECIES[species].name if space == 0 else None
        
        ppl.plot(x_values, to_plot[len(times)-1][species], color = 'black', marker=SPECIES[species].marker, markersize=4, label=label)
# griser les membranes : 

for Mb in Mb_pos:
    pos_Mb = Mb[1]
    i = pos_Mb[0]
        
    space1 = real_Rspace[i]
    space2 = space1 + thickness
     
    ppl.axvspan(space1, space2, color='grey', alpha=0.5) 

ppl.title(f"Concentration finale (times = {times[len(times)-1]} ) de {sp_name} dans l'espace")
ppl.xlabel("Position dans l'espace (x)")
ppl.ylabel("Concentration finale")
ppl.legend(loc='upper right')
ppl.grid()

# ppl.savefig("IS_euler_t01.pdf")

ppl.show()



'''Plot de la conservation'''
# integral

# ppl.figure()
# for species in range(0,NS):
#     ppl.plot(times, [array[species] for array in result_int], label=SPECIES[species].name)
# ppl.xlabel("Temps")
# ppl.ylabel("Integrale sur l'espace")
# ppl.title("Vérification de la conservation")
# ppl.grid()
# ppl.legend()


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


"""Video variation du temps params classique"""

# Plot : 2 posibilités: 
    # Log : choix des coubres en log dans le temps
    # Equal : choix des courbe de facon equidistante dans le temps
    
change_log_time = "Equal"

images = generate_images_new(Final_Grid_L, times, spacesL, Respace, Mb_pos, SPECIES, thickness, change_log_time, nb_frames=200)

# create_video_from_images(images, "Results/video/IS_osmo_t5_P10.mp4", fps=8)



