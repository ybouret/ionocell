# -*- coding: utf-8 -*-
"""
Created on Tue Dec 2 11:38:17 2024

@author: jeann

Resolution de l'equation de diffusion pour une concentration 
et la methode des elemenst finis'
"""

import numpy as np
import matplotlib.pyplot as ppl

import sys

from moviepy import ImageSequenceClip


from Functions import Initiation_functions
from Functions import Reaction_functions
from Functions import Diffusion_functions
from Functions import Transport_functions
from Functions import Make_video

import ionocell

ppl.rcParams['figure.figsize'] = [13, 6]



# """ Definitions des fonctions """

# '''CALCUL'''

# def calcul(list_Csp, list_space, dict_CI, react_dict, mesh_len, t_init, t_fin, delta_t, delta_x, TFreaction, Mb_dict, Space_info):
#     """
#     Simule l'évolution des concentrations chimiques dans un espace 2D, 
#     en prenant en compte les reactions possibles, la diffusion d'especes et le transport à travers une membrane

#     Entrées :
#     ----------
#     - list_Csp : list
#         Liste des espèces étudiées : en class
#     - list_space : list
#         List contenant les arrays: Tableaux des valeurs des mailles spatiales pour les deux espaces
#     - dict_CI : dict
#         Dictionnaire des conditions initiales pour chaque espèce
#         Format : {espèce: {"espace1": array, "espace2": array}}
#     - react_dict : dict
#         Dictionnaire des reactions
#         Format : soit pas de reaction = {}
#                 soit reaction : {"reactionN": {"reactif" : espece(s), "produit" : espece(s)} }
#     - all_dict : dict
#         Dictionnaire des differentes caracteristiques des especes :
#         Format : {"D": dict_coeffdiff, "P": dict_perm, "Z": dict_charge}
#             {espèce: valeur} pour chaque dict
#     - mesh_len: list
#         Nombre de mailles par espace
#     - t_init, t_fin : float
#         Temps initial et final de la simulation
#     - delta_t : float
#         Pas de temps pour la simulation
#     - delta_x : float
#         Pas d'espace pour la simulation
#     - TFreaction : bool
#         Indique si il y a reaction (True/False)
#     # - params_mod : list
#     #     liste des parametres de notre simulation
#     #     Format : [Type_mb]

#     Sorties :
#     ----------
#     - CS_x1_dt, CS_x2_dt : array
#         Tableaux des concentrations pour chaque espèce dans le temps et l'espace
#         Format : array([nombre de pas de temps, nombre d'espèces, nombre de mailles d'espace])
#     - integrals_TSx1, integrals_TSx2 : list
#         Intégrales spatiales des concentrations dans espace1 et espace2 pour chaque pas de temps
#         Format : list([nombre de pas de temps, nombre d'espèces])
#     - times : list
#         Liste contenant les valeurs des pas de temps

#     Fonctionnement :
#     ----------------
#     1. Initialise les conditions de départ pour chaque espèce.
#     2. Boucle sur les pas de temps :
#         - Calcule les concentrations pour le pas suivant avec reaction, diffusion et transport à travers la membrane
#         - Stocke les résultats dans les tableaux de concentrations
#         - Calcule les intégrales spatiales des concentrations
#     3. Retourne les résultats sous forme de tableaux et listes
#     """

#     Grid_CI_L = condition_initiale(list_Csp, dict_CI, list_space)  # cherche les conditions initiales des especes
    
#     t = t_init
#     times = [] 
    
#     print("La membrane modelisée est ", Mb_dict["Type"])
    
#     integrals_Grid_L = [] # analyse conservation
    
#     CS_dt = []
#     GridFinal_L = Grid_CI_L
    
#     while t < t_fin:  
        
#         CS_dt.append(GridFinal_L)
        
#         (GridFinal_L) = iteration(list_Csp, GridFinal_L, mesh_len, delta_t, delta_x, TFreaction, react_dict, Mb_dict, Space_info, t)
                
        
#         # calcul des intégrales: methodes des trapezes
#         integral_Grid_L = []
#         for space in range(0, len(Grid_CI_L)): #nb d'espace
            
#             integral_Grid = np.trapz(GridFinal_L[space], list_space[space])
#             integral_Grid_L.append(integral_Grid)
#         integrals_Grid_L.append(integral_Grid_L)
        
        
#         times.append(t)
#         t += delta_t # passe au temps suivant (a la maille de temps d'apres)   

#         print(t)
#     CS_dt.append(Grid_CI_L)
#     return (CS_dt, integrals_Grid_L, times)



# '''CONDITION INITIALES'''

# def condition_initiale(list_Csp, dict_CI, list_space):
    
#     """
#     Initialise les conditions de départ des espèces dans deux espaces.

#     Entrées :
#     ----------
#     - list_Csp : list
#         Liste des espèces étudiées , en Class
#     - dict_CI : dict
#         Dictionnaire des conditions initiales pour chaque espèce
#         Format : {espèce_x1: array, espèce_x2: array}.

#     Sorties :
#     ----------
#     - TCI_x1, TCI_x2 : list
#         Listes contenant les conditions initiales pour chaque espèce dans espace1 et espace2 respectivement
#         Format : array([nombre d'espèces, nombre de mailles d'espace])

#     Fonctionnement :
#     ----------------
#     1. Parcourt la liste des espèces (`list_sp`)
#     2. Récupère les conditions initiales de chaque espèce depuis `dict_CI` pour espace1 (`x1`) et espace2 (`x2`)
#     3. Ajoute les conditions initiales aux tableaux `TCI_x1` et `TCI_x2`
#     4. Retourne les tableaux contenant les conditions initiales pour les deux espaces
#     """
#     CI_L = []
#     for i in range (1, len(list_space)+1):
#         CI_specie = []
#         for specie in list_Csp:
#             name_sp = specie.name[0]
            
#             Val_CI = dict_CI[f"{name_sp}_x{i}"]
            
#             CI_specie.append(Val_CI)
#         CI_L.append(CI_specie)

#     return (CI_L)


# '''ITERATION'''

# def iteration(list_Csp, Grid_to_analyse_L, mesh_len, delta_t, delta_x, TFreaction, react_dict, Mb_dict, Space_info, t):
#     """
#     Calcule les concentrations dans l'espace pour le pas de temps suivant.

#     Entrées :
#     ----------
#     - list_Csp : list
#         Liste des espèces étudiées, en Class
#     - TSx1, TSx2 : array
#         Tableaux contenant les concentrations dans espace1 et espace2 au pas de temps précédent
#         Format : array([nombre d'espèces, nombre de mailles d'espace])
#     - N, M : int
#         Nombre de mailles dans espace1 (N) et espace2 (M)
#     - delta_t : float
#         Pas de temps pour la simulation
#     - delta_x : float
#         Pas d'espace pour la simulation
#     - TFreaction : bool
#         Indique si il y a reaction (True/False)
#     - react_dict : dict
#         Dictionnaire décrivant les réactions chimiques impliquant les espèces
#     - all_dict : dict
#         Dictionnaire des differentes caracteristiques des especes :
#         Format : {"D": dict_coeffdiff, "P": dict_perm, "Z": dict_charge}
#             {espèce: valeur} pour chaque dict
#     - params_mod : list
#         liste des parametres de notre simulation
#         Format : [Gradient, Type_mb]
#     - t : int
#         temps de la simulation

#     Sorties :
#     ----------
#     - F_x1, F_x2 : array
#         Tableaux des concentrations calculées pour chaque espèce dans espace1 et espace2 respectivement, au pas de temps suivant
#         Format : array([nombre d'espèces, nombre de mailles d'espace])

#     Fonctionnement :
#     ----------------
#     1. Recupère les parametres necessaires
#     2. Si il y a des réactions (`TFreaction` == True) :
#         - Calcule les nouvelles concentrations apres réaction(s) via la méthode `reaction_euler`; donne TSRx1, TSRx2
#     3. Pour chaque espèce dans `list_sp` :
#         a. Calcule la diffusion dans espace1 et espace2 avec la fonction duffision(); donne  TF_x1, TF_x2
#         b. Calcule le transport à travers la membrane (3 scenarios : osmose implicite, osmose analytique, electro-osmotique implicite)
#             choix de la methode en fonction du params_mod[Gradient]; donne F_x1, F_x2 
#     4. Retourne les tableaux `F_x1` et `F_x2` contenant les concentrations pour le pas de temps suivant
#     """
    
#     NS = len(list_Csp)
    
#     Type_mb = Mb_dict["Type"]
    
#     # step analyse : 1. reaction 2. diffusion 3. transport
#     Greact_L = []
    
    
#     # REACTION
#     for space1 in range (0, len(Grid_to_analyse_L)):
#         Grid_to_analyse = Grid_to_analyse_L[space1]

#         if TFreaction:
#             (Greact_space) = Reaction_functions.reaction_euler(list_Csp, Grid_to_analyse , react_dict)
            
#         else : # if no reaction described
#             Greact_space = Grid_to_analyse
#         Greact_L.append(Greact_space)

#     Gdiff_L = []
#     # DIFFUSION
#     for space2 in range (0, len(Greact_L)):
        
#         Grid_to_diff = Greact_L[space2]

#         nbmaille_space = mesh_len[space2]
#         delta_x_space = delta_x[space2]
        
#         Grid_diff = []
#         for num_sp1 in range(0, NS):
            
#             coeffdiff_sp = list_Csp[num_sp1].diff
            
#             if Space_info[space2].walls :
#                 (Gdiff_sp) = Diffusion_functions.diffusion_two_wall(num_sp1, Grid_to_diff, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)
            
#             elif Space_info[space2].inje :
#                 if t < 0.5 : 
#                     (Gdiff_sp) = Diffusion_functions.diffusion_one_wall_one_inject(num_sp1, Grid_to_diff, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)
#                 else : 
#                     (Gdiff_sp) = Diffusion_functions.diffusion_one_wall_one_supp(num_sp1, Grid_to_diff, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)
            
#         Grid_diff = Gdiff_sp

#         Gdiff_L.append(np.array(Grid_diff))
        
#     # TRANSPORT à la membrane
    
#     #position mb : 
#     dict_pos_mb = dict(list(Mb_dict.items())[:-1])
    
#     Gtranp_L = []
    
#     modif_c = []
#     for num_sp2 in range(0, NS):
        
#         perm_sp = list_Csp[num_sp2].perm
#         charge_sp = list_Csp[num_sp2].charge
        
#         # Grid_to_transp_sp = np.copy(Gdiff_L[:,num_sp2,:])
#         # Grid_to_transp_sp = Gdiff_L[num_sp2].copy()
#         # Grid_to_transp_sp = Gdiff_L[:,num_sp2,:]
        
#         Grid_to_transp_sp = [arr[num_sp2, :] for arr in Gdiff_L]
                
#         for Mb in dict_pos_mb:
#             pos_Mb = dict_pos_mb[Mb]
            
#             if Type_mb == "osmotique implicite":
#                 (Grid_to_transp_sp) = Transport_functions.transport_mb_osmo(Grid_to_transp_sp, pos_Mb, mesh_len, delta_t, delta_x, perm_sp)
                
#                 # (Grid_to_transp_sp) = Transport_functions.transport_mb_osmo(list_Csp, Grid_to_transp_sp, pos_Mb, mesh_len, delta_t, delta_x, perm_sp)
                
#             elif Type_mb == "osmotique analytique":
#                 (Grid_to_transp_sp) = Transport_functions.transport_mb_osmo_analyt(Grid_to_transp_sp, pos_Mb, mesh_len, perm_sp, delta_t, delta_x)

#             elif Type_mb == "electro-osmotique implicite":
#                 # a verifier
#                 (Grid_to_transp_sp) = Transport_functions.transport_mb_electro_osmo_impl(list_Csp, Grid_to_transp_sp, pos_Mb, mesh_len, delta_t, delta_x, perm_sp, charge_sp)
    
#             elif Type_mb == "non permeable":
#                 # en cours
#                 (Grid_to_transp_sp) = (Grid_to_transp_sp)
    
#             # elif Type_mb == "electro-osmotique analytique":
#                 # en cours
                
#             else :
#                 print("no transport defined")
#                 sys.exit(1)
        
#         modif_c.append(Grid_to_transp_sp)
        
#     # re shape en gardant les liste : 
#     for i in range(0, len(modif_c[0])):
#         modif_c_space = []
#         for i1 in range(0, len(modif_c)):
#             modif_c_space.append(modif_c[i1][i])
#         Gtranp_L.append(modif_c_space)
    
#     Grid_final = Gtranp_L
#     # Gtranp_L = np.stack(modif_c, axis=1)
        
#     return (Grid_final)



""" Definition des paramètres """

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

(Final_Grid_L, integrals_grid_L, times) = Initiation_functions.calcul(SPECIES, spacesL, IC_dict, d_react, nbmailles_x, ti, tf, delta_t, delta_x, TFreaction, Mb_dict, Info_space)


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
    # Log = choix des coubres en log dans le temps
    # Equal : choix des courbe de facon equidistante dans le temps
    
change_log_time = "Equal"

images = Make_video.generate_images_new(Final_Grid_L, times, spacesL, Respace, Mb_pos, SPECIES, thickness, change_log_time, nb_frames=200)

# Make_video.create_video_from_images(images, "Results/video/IS_osmo_t5_P10.mp4", fps=8)



