# -*- coding: utf-8 -*-
"""
Created on Tue Dec 2 11:38:17 2024

@author: jeann

Resolution de l'equation de diffusion pour une concentration 
par la methode des differences finis'
"""

import numpy as np
import matplotlib.pyplot as ppl

import ionocell

import sys
# sys.path.insert(0, '/Users/jleclezio/Documents/ionocell/python/Functions/')

from Initiation_functions import calcul
from Make_video import generate_images_new, create_video_from_images

from Plot_function import plot_inTime_andSpace_linear, plot_inTime_andSpace_log, plot_TimeandSpace_oneSpecies, plot_t0, plot_tLast, plot_conservation, plot_3D, plot_heatmap

ppl.rcParams['figure.figsize'] = [13, 6]



""" Definition des données d'entrées """

'''Temps de simulation'''

ti = 0
tf = 1

'''Espace'''

# 2**5, round(100*(4/9))
nbmailles_x1 = 20
nbmailles_x2 = 50
nbmailles_x3 = round(100*(4/9))

allmailles_x = [nbmailles_x1, nbmailles_x2, nbmailles_x3]

# ratio x1 vs x2: taille dans la cuve
Respace1 = 1
Respace2 = 1/2
Respace3 = 4/9
all_Respace = [Respace1, Respace2, Respace3]


# PAS D'ESPACE
delta_x1 = Respace1/(nbmailles_x1-1)
delta_x2 = Respace2/(nbmailles_x2-1)
delta_x3 = Respace3/(nbmailles_x3-1)
alldelta_x = [delta_x1, delta_x2, delta_x3]

# espace
x1 = np.arange(nbmailles_x1)* delta_x1
x2 = np.arange(nbmailles_x2)* delta_x2 
x3 = np.arange(nbmailles_x3)* delta_x3
allspacesL = [x1, x2, x3]


thickness = min(delta_x1, delta_x2)/2 #largeur de la membrane pour le plot

x1_info = ionocell.Space(_number = 0, _type = "Cell", _two_wall = True, _one_inje = False ) 
x2_info = ionocell.Space(_number = 1, _type = "Extra", _two_wall = True, _one_inje = False ) 
x3_info =  ionocell.Space(_number = 2, _type = "Cell", _two_wall = True, _one_inje = False ) 
allInfo_space = [x1_info, x2_info, x3_info]

'''Membrane'''

# nombre d'espace pour la simulation 
space_nb = 1

Membrane_list = ["osmotique implicite", "osmotique analytique", "electro-osmotique implicite", "non permeable"]

# consequences: 
allMb_dict = {
    "Mb1" : [0,1],
    # "Mb2" : [1,2], 
    "Type" : Membrane_list[1],
    "thickness" : thickness}

'''Especes'''

# descriptions
Na = ionocell.Specie(_specie_name = "Na", _diff_coeff = 1.33, _perm_coeff = 10, _charge = +1, _marker= 'o')
# Na = ionocell.Specie(_specie_name = "Na", _diff_coeff = 1.33, _perm_coeff = 1, _charge = +1, _marker= 'o')
K = ionocell.Specie(_specie_name = "K", _diff_coeff = 1, _perm_coeff = 15, _charge = +1, _marker= 's') # d= 1.96
Cl = ionocell.Specie(_specie_name = "Cl", _diff_coeff = 1, _perm_coeff = 1, _charge = -1, _marker= 's')
NaCl = ionocell.Specie(_specie_name = "NaCl", _diff_coeff = 10, _perm_coeff = 10, _charge = 0, _marker= 'd')
OH = ionocell.Specie(_specie_name = "OH", _diff_coeff = 1, _perm_coeff = 1, _charge = -1, _marker= 'd')
H = ionocell.Specie(_specie_name = "H", _diff_coeff = 1, _perm_coeff = 1, _charge = +1, _marker= 's')
H20 = ionocell.Specie(_specie_name = "H2O", _diff_coeff = 1, _perm_coeff = 1, _charge = 0, _marker= 'v')

Marker = ["o", "d", "s", "h", "v"]


# Type conditions initiales
L1 = Respace1 # longueur de la cuve: arbitraire
L2 = Respace2

var_g = 0.005 # largeur de la gaussienne
Gauss_1_4_x1 = np.exp(-((x1-(L1*1/4))**2) / var_g) 
Gauss_1_4_x2 = np.exp(-((x2-(L2*1/4))**2) / var_g)

Gauss_1_2_x1 = np.exp(-((x1-(L1*1/2))**2) / var_g) 
Gauss_1_2_x2 = np.exp(-((x2-(L2*1/2))**2) / var_g) 


Gauss_3_4_x1 = np.exp(-((x1-(L1*3/4))**2) / var_g) 
Gauss_3_4_x2 = np.exp(-((x2-(L2*3/4))**2) / var_g)

cst_1_x1 = np.ones(nbmailles_x1)
cst_1_x2 = np.ones(nbmailles_x2)

Na_i = np.ones(nbmailles_x2)*0.01
Na_e = np.ones(nbmailles_x1)*0.14

K_i = np.ones(nbmailles_x2)*0.14
K_e = np.ones(nbmailles_x1)*0.004

Null_x1 = np.zeros(nbmailles_x1) 
Null_x2 = np.zeros(nbmailles_x2)
Null_x3 = np.zeros(nbmailles_x3)

cst_05_x1 = np.ones(nbmailles_x1)*0.5
cst_05_x2 = np.ones(nbmailles_x2)*0.5


# CI pour les especes
IC_dict = {
    "Na_x1": Gauss_1_4_x1, 
    "Na_x2": Null_x2, 
    "Na_x3": Null_x3, 
    
    
    "K_x1": K_i, 
    "K_x2": K_e, 
    
    "Cl_x1": Gauss_3_4_x1, 
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



#pour notre simu : 
    # liste des especes utilisées poru la simulation (decrites: "Na", "Cl", "OH", "NaCl")
    # si pas decrite: rajouter les caracteristiques dans les dict + les CI
    # SEULE LISTE A MODIFIER POUR AJOUTER/ENLEVER DES ESPECES POUR UNE SIMULATION
SPECIES = [Na, Cl, NaCl] 


'''Temps'''

# choix du dt 

prec = 0.1

dtmax_allsp = []
for sp in SPECIES:
    D = sp.diff
    dtmax_sp = []
    for es in range(0,space_nb):
    
        dt = (alldelta_x[es]**2) / (2*D)
        
        dtmax_sp.append(prec * dt)
    dtmax_allsp.append(dtmax_sp)
    

delta_t = min(min(sublist) for sublist in dtmax_allsp)


'''nombre d'espace'''

spacesL = allspacesL[:space_nb]
nbmailles_x = allmailles_x[:space_nb]
# delta_t = min(alldelta_t[:space_nb])
delta_x = alldelta_x[:space_nb]
Info_space = allInfo_space[:space_nb]
Respace = all_Respace[:space_nb]

Mb_pos = list(allMb_dict.items())[:space_nb-1]
type_mb = list(allMb_dict.items())[-2]
Thickness = list(allMb_dict.items())[-1]


Mb_dict = dict(Mb_pos + [type_mb] + [Thickness])



'''reaction'''

TFreaction = True

d_react = {
    "reaction0": {"reactif" : ["Na", "Cl"] , \
                  "produit" : ["NaCl"],  \
                  "cst_eq" : 1} }

# TFreaction = False
# d_react = {}


""" Calculs """

(Final_Grid_L, integrals_grid_L, times) = calcul(SPECIES, spacesL, IC_dict, d_react, nbmailles_x, ti, tf, delta_t, delta_x, TFreaction, Mb_dict, Info_space)


integrals_grid_L = np.array(integrals_grid_L)

NS = len(SPECIES)


""" Analyse de la conservation """  
# seuil de tolérance
tolerance = 0.05

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


'''parametres plot'''

#pas de temps:
nb_courbe_t = 8    # PEUT ETRE MODIFIÉ POUR PLOT PLUS DE COURBE = PLUS DE TEMPS


# description des couleurs en fct du nb de courbe
cmap = ppl.colormaps['Dark2']  
discrete_cmap = cmap.resampled(nb_courbe_t) 
colors = [discrete_cmap(i) for i in range(nb_courbe_t)]


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
    
# recup element en lispace   
indices = np.linspace(0, (len(times) - 1), nb_courbe_t, dtype=int)
times_courbelin = np.array(times)[indices]


params_plot = {
    "colors" : colors, 
    "real_spacesL" : real_spacesL, 
    "real_Rspace" : real_Rspace,
    "indices" : indices,
    "times_courbelin" : times_courbelin,
    "log_indices" : log_indices,
    "times_courbelog" : times_courbelog, 
    "log_indices_2log" : log_indices_2log,
    "times_courbelog_2log" : times_courbelog_2log
    }

'''Plot final'''

# plot all species, in time and space, 2 scale of time : linear or log
plot_inTime_andSpace_linear(Final_Grid_L, params_plot, spacesL, SPECIES, Mb_pos, thickness)

plot_inTime_andSpace_log(Final_Grid_L, params_plot, spacesL, SPECIES, Mb_pos, thickness)

# plot for one species, scale time in log
species_to_plot = "Na"
plot_TimeandSpace_oneSpecies(Final_Grid_L, species_to_plot, params_plot, spacesL, SPECIES, Mb_pos, real_Rspace, thickness)

# plot all species, initial condition, t=0 
plot_t0(Final_Grid_L, spacesL, params_plot, SPECIES, Mb_pos, thickness)

# plot all species, final condition t=tfin
plot_tLast(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness)
    
# analyse de la conservation
plot_conservation(result_int, variations_relative, times, SPECIES)


"""Video variation du temps params classique"""

# Plot : 2 posibilités: 
    # Log : choix des coubres en log dans le temps
    # Equal : choix des courbe de facon equidistante dans le temps
    
change_log_time = "Log"

images = generate_images_new(Final_Grid_L, times, spacesL, Respace, Mb_pos, SPECIES, thickness, change_log_time, nb_frames=100)

# create_video_from_images(images, "Results/video/Diff_simple_Na.mp4", fps=8)


""" 3D plot"""

plot_3D(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness)


"""HEAT MAP"""

plot_heatmap(Final_Grid_L, times, spacesL, params_plot, SPECIES, Mb_pos, thickness )
