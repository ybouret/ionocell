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

from Functions import Reaction_functions
from Functions import Diffusion_functions
from Functions import Transport_functions

ppl.rcParams['figure.figsize'] = [13, 6]


""" Definitions des fonctions """

'''CALCUL'''

def calcul(list_sp, x1, x2, dict_CI, react_dict, all_dict, nbmailles_x, t_init, t_fin, delta_t, delta_x, TFreaction, params_mod):
    """
    Simule l'évolution des concentrations chimiques dans un espace 2D, 
    en prenant en compte les reactions possibles, la diffusion d'especes et le transport à travers une membrane

    Entrées :
    ----------
    - list_sp : list
        Liste des noms des espèces étudiées 
    - x1, x2 : array
        Tableaux des valeurs des mailles spatiales pour les deux espaces
    - dict_CI : dict
        Dictionnaire des conditions initiales pour chaque espèce
        Format : {espèce: {"espace1": array, "espace2": array}}
    - react_dict : dict
        Dictionnaire des reactions
        Format : soit pas de reaction = {}
                soit reaction : {"reactionN": {"reactif" : espece(s), "produit" : espece(s)} }
    - all_dict : dict
        Dictionnaire des differentes caracteristiques des especes :
        Format : {"D": dict_coeffdiff, "P": dict_perm, "Z": dict_charge}
            {espèce: valeur} pour chaque dict
    - nbmailles_x: list
        Nombre de mailles par espace
    - t_init, t_fin : float
        Temps initial et final de la simulation
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - TFreaction : bool
        Indique si il y a reaction (True/False)
    - params_mod : list
        liste des parametres de notre simulation
        Format : [Type_mb]

    Sorties :
    ----------
    - CS_x1_dt, CS_x2_dt : array
        Tableaux des concentrations pour chaque espèce dans le temps et l'espace
        Format : array([nombre de pas de temps, nombre d'espèces, nombre de mailles d'espace])
    - integrals_TSx1, integrals_TSx2 : list
        Intégrales spatiales des concentrations dans espace1 et espace2 pour chaque pas de temps
        Format : list([nombre de pas de temps, nombre d'espèces])
    - times : list
        Liste contenant les valeurs des pas de temps

    Fonctionnement :
    ----------------
    1. Initialise les conditions de départ pour chaque espèce.
    2. Boucle sur les pas de temps :
        - Calcule les concentrations pour le pas suivant avec reaction, diffusion et transport à travers la membrane
        - Stocke les résultats dans les tableaux de concentrations
        - Calcule les intégrales spatiales des concentrations
    3. Retourne les résultats sous forme de tableaux et listes
    """
    (Grid_CI_x1, Grid_CI_x2) = condition_initiale(list_sp, dict_CI)  # cherche les conditions initiales des especes
    
    t = t_init
    times = []    
    
    print("La membrane modelisée est ", params_mod[0])
    
    integrals_Gridx1 = []
    integrals_Gridx2 = []
    
    (GridInit_x1, GridInit_x2) =  (Grid_CI_x1, Grid_CI_x2)
        
    CS_x1_dt = np.array([GridInit_x1])
    CS_x2_dt = np.array([GridInit_x2])

    while t < t_fin:  
        (GridFinal_x1, GridFinal_x2) = iteration(list_sp, GridInit_x1, GridInit_x2, nbmailles_x, delta_t, delta_x, TFreaction, react_dict, all_dict, params_mod)
                        # calcul les concentrations dans l'espace au temps d'apres
        
        CS_x1_dt = np.append(CS_x1_dt, [GridFinal_x1], axis=0)
              
        CS_x2_dt = np.append(CS_x2_dt, [GridFinal_x2], axis=0)
            
        # calcul des intégrales: methodes des trapezes
        integral_Gridx1 = np.trapz(GridFinal_x1, x1)
        integral_Gridx2 = np.trapz(GridFinal_x2, x2)
       
        integrals_Gridx1.append(integral_Gridx1)
        integrals_Gridx2.append(integral_Gridx2)
        
        # print(t)

        times.append(t)
        t += delta_t # passe au temps suivant (a la maille de temps d'apres)   
        
    return (CS_x1_dt, CS_x2_dt, integrals_Gridx1, integrals_Gridx2, times)



'''CONDITION INITIALES'''

def condition_initiale(list_sp, dict_CI):
    
    """
    Initialise les conditions de départ des espèces dans deux espaces.

    Entrées :
    ----------
    - list_sp : list
        Liste des noms des espèces étudiées 
    - dict_CI : dict
        Dictionnaire des conditions initiales pour chaque espèce
        Format : {espèce_x1: array, espèce_x2: array}.

    Sorties :
    ----------
    - TCI_x1, TCI_x2 : list
        Listes contenant les conditions initiales pour chaque espèce dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Parcourt la liste des espèces (`list_sp`)
    2. Récupère les conditions initiales de chaque espèce depuis `dict_CI` pour espace1 (`x1`) et espace2 (`x2`)
    3. Ajoute les conditions initiales aux tableaux `TCI_x1` et `TCI_x2`
    4. Retourne les tableaux contenant les conditions initiales pour les deux espaces
    """
    
    
    CI_x1 = []
    CI_x2 = []
    
    for specie in list_sp:
        
        Val_CI_x1 = dict_CI[f"{specie}_x1"]
        Val_CI_x2 = dict_CI[f"{specie}_x2"]
        
        CI_x1.append(Val_CI_x1)
        CI_x2.append(Val_CI_x2)
        
    return (CI_x1, CI_x2)


'''ITERATION'''

def iteration(list_sp, Grid_to_analyse_x1, Grid_to_analyse_x2, nbmailles_x, delta_t, delta_x, TFreaction, react_dict, all_dict, params_mod):
    """
    Calcule les concentrations dans l'espace pour le pas de temps suivant.

    Entrées :
    ----------
    - list_sp : list
        Liste des noms des espèces étudiées
    - TSx1, TSx2 : array
        Tableaux contenant les concentrations dans espace1 et espace2 au pas de temps précédent
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - TFreaction : bool
        Indique si il y a reaction (True/False)
    - react_dict : dict
        Dictionnaire décrivant les réactions chimiques impliquant les espèces
    - all_dict : dict
        Dictionnaire des differentes caracteristiques des especes :
        Format : {"D": dict_coeffdiff, "P": dict_perm, "Z": dict_charge}
            {espèce: valeur} pour chaque dict
    - params_mod : list
        liste des parametres de notre simulation
        Format : [Gradient, Type_mb]
    - t : int
        temps de la simulation

    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations calculées pour chaque espèce dans espace1 et espace2 respectivement, au pas de temps suivant
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Recupère les parametres necessaires
    2. Si il y a des réactions (`TFreaction` == True) :
        - Calcule les nouvelles concentrations apres réaction(s) via la méthode `reaction_euler`; donne TSRx1, TSRx2
    3. Pour chaque espèce dans `list_sp` :
        a. Calcule la diffusion dans espace1 et espace2 avec la fonction duffision(); donne  TF_x1, TF_x2
        b. Calcule le transport à travers la membrane (3 scenarios : osmose implicite, osmose analytique, electro-osmotique implicite)
            choix de la methode en fonction du params_mod[Gradient]; donne F_x1, F_x2 
    4. Retourne les tableaux `F_x1` et `F_x2` contenant les concentrations pour le pas de temps suivant
    """
    
    NS = len(list_sp)

    coeffdiff_dict = all_dict["D"]
    perm_dict = all_dict["P"]
    charge_dict = all_dict["Z"]
    
    # Gradient = params_mod[0]
    Type_mb = params_mod[0]
    
    # achangere pour automatisation
    nbmailles_x1 = nbmailles_x[0]
    nbmailles_x2 = nbmailles_x[1]
    
    # step analyse : 1. reaction 2. diffusion 3. transport

    # REACTION
    if TFreaction:
        (Greactx1) = Reaction_functions.reaction_euler(list_sp, Grid_to_analyse_x1, react_dict)
        (Greactx2) = Reaction_functions.reaction_euler(list_sp, Grid_to_analyse_x2, react_dict)
        
    else : # if no reaction described
        Greactx1 = Grid_to_analyse_x1
        Greactx2 = Grid_to_analyse_x2
    
    for species in range(0, NS):
        # DIFFUSION
        (Gdiff_x1) = Diffusion_functions.diffusion(list_sp, species, Greactx1, nbmailles_x1, delta_t, delta_x1, coeffdiff_dict)
        (Gdiff_x2) = Diffusion_functions.diffusion(list_sp, species, Greactx2, nbmailles_x2, delta_t, delta_x2, coeffdiff_dict)
        
        # TRANSPORT à la membrane
        if Type_mb == "osmotique implicite":
            (Gtransp_x1, Gtransp_x2) = Transport_functions.transport_mb_osmo(list_sp, species, Gdiff_x1, Gdiff_x2, nbmailles_x, delta_t, delta_x, perm_dict)
            
        elif Type_mb == "osmotique analytique":
            (Gtransp_x1, Gtransp_x2) = Transport_functions.transport_mb_osmo_analyt(list_sp, species, Gdiff_x1, Gdiff_x2, nbmailles_x, perm_dict, delta_t, delta_x)
                
        elif Type_mb == "electro-osmotique implicite":
            (Gtransp_x1, Gtransp_x2) = Transport_functions.transport_mb_electro_osmo_impl(list_sp, species, Gdiff_x1, Gdiff_x2, nbmailles_x, delta_t, delta_x, perm_dict, charge_dict)

        elif Type_mb == "non permeable":
            (Gtransp_x1, Gtransp_x2) = (Gdiff_x1, Gdiff_x2)

        # elif Type_mb == "electro-osmotique analytique":
        else :
            print("no transport defined")
            sys.exit(1)
            
    # no more step:
    GFinal_x1 = Gtransp_x1
    GFinal_x2 = Gtransp_x2
    
    return (GFinal_x1, GFinal_x2)



""" Definition des paramètres """

'''parametres modifiables'''

nbmailles_x1 = 2**5
nbmailles_x2 = 2**5
nbmailles_x = [nbmailles_x1, nbmailles_x2]

# ratio x1 vs x2: taille dans la cuve
Respace1 = 1/2
Respace2 = 1-(Respace1) # 1 car taille cuve = 1

# TEMPS: de 0 à 1
ti = 0
tf = 0.5

# CONSERVATION : seuil de tolérance
tolerance = 0.05

# CI GAUSSIENNE: especes decritent à t=0 par gaussienne
var_g = 0.005 # largeur de la gaussienne


'''Calcul pas de temps et d'espaces'''

# PAS D'ESPACE
delta_x1 = Respace1/(nbmailles_x1-1)
delta_x2 = Respace2/(nbmailles_x2-1)
delta_x = [delta_x1, delta_x2]

# espace
x1 = np.arange(nbmailles_x1)* delta_x1
x2 = np.arange(nbmailles_x2)* delta_x2 


# PAS DE TEMPS
delta_t1 = 0.9 * (delta_x1**2) / 2 # 0.9 * delta_t max 
delta_t2 = 0.9 * (delta_x2**2) / 2 # 0.9 * delta_t max 
delta_t = min(delta_t1,delta_t2) # recuperer le delta_t min = celui avec le plus de precisions


thickness = max(delta_x1, delta_x2)/2 #largeur de la membrane pour le plot

'''description d'especes'''

# liste des especes utilisées poru la simulation (decrites: "Na", "Cl", "OH", "NaCl")
# si pas decrite: rajouter les caracteristiques dans les dict + les CI
# SEULE LISTE A MODIFIER POUR AJOUTER/ENLEVER DES ESPECES POUR UNE SIMULATION
SPECIES = ["Na"] 

Marker = ["o", "d", "s", "h", "v"]


# COEFF DE DIFFUSION
d_coeffdiff = {
    "Na": 1,
    "Cl": 1,
    "NaCl": 1, 
    "OH": 1, 
    "H": 1, 
    "H2O": 1}

# PERMEABILITE
d_perm = {
    "Na": 10, # 1e-8,
    "Cl": 10, # 1e-7, 
    "OH": 1, # 1e-10,
    "NaCl": 10, 
    "H": 1, 
    "H2O": 1}

# CHARGE 
d_charge = {
    "Na": +1, # 1e-8,
    "Cl": -1, # 1e-7, 
    "H2O": 0, 
    "OH": -1, 
    "H": +1, # 1e-10,
    "NaCl": 0}

all_dict = {
    "D": d_coeffdiff, 
    "P": d_perm, 
    "Z": d_charge}

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

IC_dict = {
    "Na_x1": Gauss_1_4_x1, 
    "Na_x2": Null_x2, 
    "Cl_x1": cst_1_x1, 
    "Cl_x2": Null_x2, 
    "NaCl_x1": Null_x1, 
    "NaCl_x2": Null_x2,
    "OH_x1": Null_x1, 
    "OH_x2": Gauss_3_4_x2, 
    "H_x1": Gauss_1_4_x1, 
    "H_x2": Null_x2, 
    "H2O_x1": Null_x1, 
    "H2O_x2": Null_x2
    }

'''Params model pour calcul'''

# Calcul_gradient_avant = ["reaction", "diffusion", "transport"]
# Gradient = Calcul_gradient_avant[2]

Membrane_list = ["osmotique implicite", "osmotique analytique", "electro-osmotique implicite", "non permeable"]
Type_mb = Membrane_list[1]

Params_calcul = [Type_mb]

# Caract_space = {
#     "x1" : {
#         "BordG" : "W", 
#         "BordD" : "Mb"
#         }, 
#     "x2": {
#         "BordG" : "Mb",
#         "BordD" : "W"
#         }
#     }

# Caract_Mb = {
#     "MB1" : {
#         "BordG" : 0, 
#         "BordD" : 1
#         }
#     }


""" Calculs """

(TCF_x1, TCF_x2, integrals_TCx1, integrals_TCx2, times) = calcul(SPECIES, x1, x2, IC_dict, d_react, all_dict, nbmailles_x, ti, tf, delta_t, delta_x, TFreaction, Params_calcul)


NS = len(SPECIES)


""" Analyse de la conservation """    
integrals_TCF_tot = []
for i in range(0, np.shape(integrals_TCx1)[0]):
    integrals_TCF_tot.append(integrals_TCx1[i]+integrals_TCx2[i])

initial_TCF = integrals_TCF_tot[0]
variations_relative = [abs((s - initial_TCF) / initial_TCF) for s in integrals_TCF_tot]

for species in range(0,NS):    
    variations_relative_specie = [array[species] for array in variations_relative]
    non_conservative = any(var > tolerance for var in variations_relative_specie)
    
    if non_conservative:
        print("Le modèle avec", SPECIES[species], "n'est pas conservatif.")
    else:
        print("Le modèle", SPECIES[species], "est conservatif.")


""" Plot resultats """

#pas de temps:
nb_courbe_t = 8    # PEUT ETRE MODIFIÉ POUR PLOT PLUS DE COURBE = PLUS DE TEMPS


# description des couleurs en fct du nb de courbe
cmap = ppl.colormaps['Dark2']  
discrete_cmap = cmap.resampled(nb_courbe_t) 
colors = [discrete_cmap(i) for i in range(nb_courbe_t)]

'''parametres plot'''

#recuperation des elements espacés en log
log_indices = np.logspace(0, np.log10(len(times)-1), nb_courbe_t, dtype=int) 
times_courbelog = np.array(times)[log_indices.astype(int)]
                                  

TCF_x1_courbelog = TCF_x1[log_indices.astype(int)]
TCF_x2_courbelog = TCF_x2[log_indices.astype(int)]

TCF_courbelog = np.concatenate((TCF_x1_courbelog, TCF_x2_courbelog), axis = 2)


espace_2cuves = np.concatenate((x1, x2+Respace1), axis = None)
espace_2cuves_mb = espace_2cuves


espace_2cuves_mb[int(nbmailles_x1):] += thickness

    # au premier temps
TCF0 = np.concatenate((TCF_x1, TCF_x2), axis = 2)

    # au dernier temps
TCFF = np.concatenate((TCF_x1, TCF_x2), axis = 2)


'''Plot final'''

ppl.figure() 
for species in range(0,NS):
    ppl.plot(espace_2cuves, TCF0[0][species], color = 'black', marker=Marker[species], markersize=4, label=SPECIES[species])
    for i, t in enumerate(times_courbelog):
        ppl.plot(espace_2cuves, TCF_courbelog[i][species], color = colors[i], marker=Marker[species], markersize=4)

ppl.axvspan(Respace1, Respace1 + thickness, color='grey', alpha=0.5)  # griser la membrane
    
ppl.plot([], [], color='black', label='t=0')  
for i,t in enumerate(times_courbelog):
    ppl.plot([], [], color=colors[i], label=f't={t:.4f}')

ppl.title(f"Concentration des especes {SPECIES} dans l'espace au cours du temps")
ppl.xlabel("Position dans l'espace (x)")
ppl.ylabel("Concentration")
ppl.legend()
ppl.grid()

# savefig("electro_osmo_impl_Na_t1_Vm_p70.pdf")

ppl.show()

'''plot au premier temps (CI)'''
# ppl.figure()
# for species in range(0,NS):
#     ppl.plot(espace_2cuves, TCFF[0][species], color = 'black', marker=Marker[species], markersize=4, label=SPECIES[species])
    
# ppl.title(f"Concentration finale (times = {times[len(times)-1]} ) de A en fonction de l'espace")
# ppl.xlabel("Position dans l'espace (x)")
# ppl.ylabel("Concentration finale")
# ppl.legend()
# ppl.grid()



'''Plot au dernier temps'''
ppl.figure()
for species in range(0,NS):
    ppl.plot(espace_2cuves, TCFF[-1][species], marker=Marker[species], markersize=4, label=SPECIES[species])
    
ppl.title(f"Concentration finale (times = {times[len(times)-1]} ) de {SPECIES} dans l'espace")
ppl.xlabel("Position dans l'espace (x)")
ppl.ylabel("Concentration finale")
ppl.legend()
ppl.grid()

# ppl.savefig("electroosmo_n70_NaCl_t1_P10.pdf")



'''Plot dans l'espace 2'''
# ppl.figure() 
# for species in range(0,NS):
#     ppl.plot(x2, TCF_x2_courbelog[0][species], color = 'black', marker=Marker[species], markersize=4, label=SPECIES[species])
#     for i, t in enumerate(times_courbelog):
#         ppl.plot(x2, TCF_x2_courbelog[i][species], color = colors[i], marker=Marker[species], markersize=4)
        
# ppl.plot([], [], color='black', label='t=0')
# for i, t in enumerate(times_courbelog):
#     ppl.plot([], [], color=colors[i], label=f't={t:.4f}')

# ppl.title("Concentration de Aen fonction de l'espace x2 pour différents temps")
# ppl.xlabel("Position dans l'espace 2 (x2)")
# ppl.ylabel("Concentration")
# ppl.legend()
# ppl.grid()


'''Plot de la conservation'''
# integral
# ppl.figure()
# for species in range(0,NS):
#     ppl.plot(times, [array[species] for array in integrals_TCF_tot], label=SPECIES[species])
# ppl.xlabel("Temps")
# ppl.ylabel("Integrale sur l'espace")
# ppl.title("Vérification de la conservation")
# ppl.grid()
# ppl.legend()


# variations relatives
# ppl.figure()
# for species in range(0,NS):
#     ppl.plot(times, [array[species] for array in variations_relative], label=SPECIES[species])
# ppl.xlabel("Temps")
# ppl.ylabel("Variations relative")
# ppl.title("Vérification de la conservation")
# ppl.legend()
# ppl.grid()


ppl.show()


"""Video variation du temps params classique"""
nb_coubre_video = 50

#recuperation des elements espacés en log
log_indices_v = np.logspace(0, np.log10(len(times)-1), nb_coubre_video, dtype=int) 
times_courbelog_v = np.array(times)[log_indices_v.astype(int)]
                                  

TCF_x1_courbelog_v = TCF_x1[log_indices_v.astype(int)]
TCF_x2_courbelog_v = TCF_x2[log_indices_v.astype(int)]

TCF_courbelog_v = np.concatenate((TCF_x1_courbelog_v, TCF_x2_courbelog_v), axis = 2)

data = TCF_courbelog_v

# description des couleurs en fct du nb de courbe
cmap = ppl.get_cmap("tab20")  
colors_v = [cmap(i / 50) for i in range(50)]


def generate_images(data):
    num_frames, num_curves, num_points = data.shape  # (temps, courbes, points)
    x = espace_2cuves
    
    images = []
    for t in range(num_frames):
        fig, ax = ppl.subplots()
        
        if 0 < t < num_frames-1:
            for i in range(num_curves):
                ax.plot(x, data[t, i, :], color = colors_v[t], marker=Marker[i], markersize=4, label=f't={times_courbelog_v[t]:.4f}')
                ax.plot(espace_2cuves, TCF0[0][i], color = 'black', marker=Marker[i], markersize=4, label=f'{SPECIES[i]}, t=0')
            # ax.plot([], [], color=colors_v[t], label=f't={times_courbelog_v[t]:.4f}')    
            time_round = round(times_courbelog_v[t],5)
            
            ax.set_xlabel("Position dans l'espace (x)")
            ax.set_ylabel("Concentration")
            ax.set_title(f'times = {time_round}')
        
        elif t == 0 :
            for i in range(num_curves):
                ax.plot(espace_2cuves, TCF0[0][i], color = 'black', marker=Marker[i], markersize=4, label=f'{SPECIES[i]}, t=0')
                
            ax.set_xlabel("Position dans l'espace (x)")
            ax.set_ylabel("Concentration")
            ax.set_title(f'Concentration de {SPECIES[i]} à t=0')
            
        else :
            for i in range(num_curves):
                ax.plot(espace_2cuves, TCF0[0][i], color = 'black', marker=Marker[i], markersize=4, label=f'{SPECIES[i]}, t=0')
                
                for a in enumerate(times_courbelog_v[:-1]):
                    ax.plot(espace_2cuves, data[a][i], color = colors_v[a], marker=Marker[i], markersize=4)
                 
                ax.plot(espace_2cuves, TCFF[len(times)-1][i], color = 'black', marker=Marker[i], markersize=4, label=f'{SPECIES[i]}, {round(times_courbelog_v[-1],5)}')
                
            for c, d in enumerate(times_courbelog_v):
                if c %3 == 1:
                    ax.plot([], [], color=colors_v[c], label=f't={d:.4f}')
                    
                    
            time_round = round(times_courbelog_v[t],5)
                
            ax.set_xlabel("Position dans l'espace (x)")
            ax.set_ylabel("Concentration")
            ax.set_title(f'Evolution de la concentration de {SPECIES[i]} pour t entre [0,{tf}]')
            
        ax.axvspan(Respace1, Respace1 + thickness, color='grey', alpha=0.5)  # griser la membrane
        ax.legend()
        ax.set_ylim(0, 1)  # peut etre modifie
        ax.grid()
            
        
            
        # Convertir la figure en image NumPy
        fig.canvas.draw()
        image = np.array(fig.canvas.renderer.buffer_rgba())
        images.append(image)
        
        ppl.close(fig)

    return images

# Convertir les images en vidéo
def create_video_from_images(images, output_file, fps=8):
    clip = ImageSequenceClip(images, fps=fps)
    clip.write_videofile(output_file, codec='libx264')

# Générer les images et créer la vidéo
# images = generate_images(data)
# create_video_from_images(images, 'transport_electroosmo_impl_Na_n70_t1.mp4')


