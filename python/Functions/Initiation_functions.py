#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  9 10:38:56 2025

@author: jleclezio
"""
import numpy as np
import sys
import copy
import math

sys.path.insert(0, '/Users/jleclezio/Documents/ionocell/python/Functions/')

from Reaction_functions import reaction_euler
from Diffusion_functions import diffusion_two_wall, diffusion_one_wall_one_inject, diffusion_one_wall_one_supp
from Transport_functions import non_transport, transport_mb_osmo, transport_mb_osmo_analyt, transport_mb_electro_osmo_impl, transport_mb_electro_osmo_NaK, flux_electro_osmo_sp, flux_mb_osmo

"""Calcul"""

def calcul(list_Csp, list_space, dict_CI, react_dict, mesh_len, t_init, t_fin, delta_t, delta_x, TFreaction, Mb_dict, Space_info, L_transp, dict_transporter, DIFF, TRANSP):
    """
    discretisation dans le temps
    Simule l'évolution des concentrations chimiques dans un espace 2D, 
    en prenant en compte les reactions possibles, la diffusion d'especes et le transport à travers une membrane

    Entrées :
    ----------
    - list_Csp : list
        Liste des espèces étudiées : en class
    - list_space : list
        List contenant les arrays: Tableaux des valeurs des mailles spatiales pour les deux espaces
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
    - mesh_len: list
        Nombre de mailles par espace
    - t_init, t_fin : float
        Temps initial et final de la simulation
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - TFreaction : bool
        Indique si il y a reaction (True/False)
    # - params_mod : list
    #     liste des parametres de notre simulation
    #     Format : [Type_mb]

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

    Grid_CI_L = condition_initiale(list_Csp, dict_CI, list_space)  # cherche les conditions initiales des especes
    
    t = t_init
    times = [] 
    # t_fin = t_init + delta_t*5
    
    print("La membrane modelisée comprend ces transporteurs : ", L_transp)
    
    integrals_Grid_L = [] # analyse conservation
    
    CS_dt = []
    GridFinal_L = Grid_CI_L
    
    decimals_t = max(2, -int(math.floor(math.log10(t_fin))) + 2)
    last_print = None
    
    while t < t_fin:  
        
        CS_dt.append(GridFinal_L)
        
        GridFinal_L = iteration(list_Csp, GridFinal_L, mesh_len, delta_t, delta_x, TFreaction, react_dict, Mb_dict, Space_info, t, L_transp, dict_transporter, DIFF, TRANSP)
        
        # calcul des intégrales: methodes des trapezes
        integral_Grid_L = []
        for space in range(0, len(Grid_CI_L)): #nb d'espace
            
            integral_Grid = np.trapz(GridFinal_L[space], list_space[space])
            integral_Grid_L.append(integral_Grid)
        integrals_Grid_L.append(integral_Grid_L)
        
        
        times.append(t)
        t += delta_t # passe au temps suivant (a la maille de temps d'apres)   
        
        # print t :  if t - t+deltat > 1e-2
        val = round(t, decimals_t)
        if val != last_print:
            print(f"{val:.{decimals_t}f}")
            last_print = val
        
    CS_dt.append(Grid_CI_L)
    
    return (CS_dt, integrals_Grid_L, times)



'''CONDITION INITIALES'''

def condition_initiale(list_Csp, dict_CI, list_space):
    
    """
    Initialise les conditions de départ des espèces dans deux espaces.

    Entrées :
    ----------
    - list_Csp : list
        Liste des espèces étudiées , en Class
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
    CI_L = []
    for i in range (1, len(list_space)+1):
        CI_specie = []

        for specie in list_Csp:
            name_sp = specie.name[0]
            
            Val_CI = dict_CI[f"{name_sp}_x{i}"]

            CI_specie.append(Val_CI)
        CI_L.append(CI_specie)

    return (CI_L)


'''ITERATION'''

def iteration(list_Csp, Grid_to_analyse_L, mesh_len, delta_t, delta_x, TFreaction, react_dict, Mb_dict, Space_info, t, L_transp, dict_transporter, DIFF, TRANSP):
    """
    discretisation dans l'espace
    Calcule les concentrations dans l'espace pour le pas de temps suivant.

    Entrées :
    ----------
    - list_Csp : list
        Liste des espèces étudiées, en Class
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
    
    NS = len(list_Csp)
    
    # step analyse : 1. reaction 2. diffusion 3. transport
    Greact_L = []
    
    # REACTION
    for space1 in range (0, len(Grid_to_analyse_L)):
        Grid_to_analyse = Grid_to_analyse_L[space1]

        if TFreaction== True:
            (Greact_space) = reaction_euler(list_Csp, Grid_to_analyse , react_dict)

        else : # if no reaction described
            Greact_space = Grid_to_analyse
            
        Greact_L.append(Greact_space)
        
        
    if DIFF == True : 
        Gdiff_L = []
        # DIFFUSION
        for space2 in range (0, len(Greact_L)):
    
            Grid_to_diff = Greact_L[space2]
            
            nbmaille_space = mesh_len[space2]
            delta_x_space = delta_x[space2]
    
            Grid_diff = []
            for num_sp1 in range(0, NS):

                Grid_to_diff_2sp = Grid_to_diff[num_sp1]
                coeffdiff_sp = list_Csp[num_sp1].diff
                
                # if Space_info[space2].walls :
                (Gdiff_sp) = diffusion_two_wall(num_sp1, Grid_to_diff_2sp, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)

                # elif Space_info[space2].inje :
                #     if t < 0.5 : 
                #         (Gdiff_sp) = diffusion_one_wall_one_inject(num_sp1, Grid_to_diff, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)
                #     else : 
                #         (Gdiff_sp) = diffusion_one_wall_one_supp(num_sp1, Grid_to_diff, nbmaille_space, delta_t, delta_x_space, coeffdiff_sp)
                # else:
                #     Gdiff_sp = Grid_to_diff
                
                Grid_diff.append(Gdiff_sp)
    
            Gdiff_L.append(np.array(Grid_diff))
    else : 
        Gdiff_L = Greact_L
    
    
    # TRANSPORT à la membrane
    if TRANSP == True : 
        
        # position mb : 
        pos_Mb_L=[]
        for i in range (1, Mb_dict["nb_mb"]+1):
            key = f"Mb{i}"
            pos_Mb_L.append(Mb_dict[key])
            
        Em = Mb_dict["potential"]
    
        Grid_transp_sp = copy.deepcopy(Gdiff_L)
    
        for nb_membrane in range(0,Mb_dict["nb_mb"]) :
            pos_Mb = pos_Mb_L[nb_membrane]
            
            deltax1 = delta_x[pos_Mb[0]]
            deltax2 = delta_x[pos_Mb[1]]
            
            c_int_L = [row[-1] for row in Gdiff_L[pos_Mb[0]]] 
            c_ext_L = [row[0] for row in Gdiff_L[pos_Mb[1]]] 
            
            c_int_L_n = c_int_L.copy()
            c_ext_L_n = c_ext_L.copy()
            
            Flux_elec = {}
    
            # for each species
            for nb_sp, species in enumerate(list_Csp):  
                
                name_sp = species.name[0]
                
                C_int = c_int_L[nb_sp]
                C_ext = c_ext_L[nb_sp]
                
                # electro osmotic flux
                J_elect_sp = flux_electro_osmo_sp(species, C_int, C_ext, pos_Mb, Em, mesh_len, delta_t, delta_x)
    
                Flux_elec[name_sp] = J_elect_sp
    
                # transporter
                # if sp in multiple transporter : not working
                for transp in L_transp:
                    if transp == None :
                        flux_transp = 0
                    elif name_sp in dict_transporter[transp]:
                        
                        flux_transp = dict_transporter[transp][name_sp] * dict_transporter[transp]["Rho"]
                        
                    else:
                        flux_transp = 0
                      
                flux_sp_tot = J_elect_sp + flux_transp
                # print(flux_sp_tot, name_sp)
                
                if Mb_dict["Type"] == "non permeable":
                    flux_sp_tot = 0
                    
                if Mb_dict["Type"] == "osmotique":
                    flux_sp_tot = flux_mb_osmo(species, C_int, C_ext, pos_Mb, mesh_len, delta_t, delta_x)
                    
                # value after transport for border :        
                C_int_n = C_int - (flux_sp_tot*delta_t)/deltax1
                C_ext_n = C_ext + (flux_sp_tot*delta_t)/deltax2
                # print(deltax1, 'deltax1')
                # print(deltax2, "deltax2")
                
                c_int_L_n[nb_sp] = C_int_n
                c_ext_L_n[nb_sp] = C_ext_n      
                
                
    
            for i, (val_int, val_ext) in enumerate(zip(c_int_L_n, c_ext_L_n)):
    
                Grid_transp_sp[pos_Mb[0]][i][-1] = val_int
                Grid_transp_sp[pos_Mb[1]][i][0] = val_ext
    else : 
        Grid_transp_sp = Gdiff_L
        
    Grid_final = Grid_transp_sp

    # modif_c = Grid_transp_sp
   
    # # re shape en gardant les liste : 
    # for i in range(0, len(modif_c[0])):
    #     modif_c_space = []
    #     for i1 in range(0, len(modif_c)):
    #         modif_c_space.append(modif_c[i1][i])
    #     Gtranp_L.append(modif_c_space)
    
    # Grid_final = Gtranp_L
    # print(np.shape(Grid_final), "Grid_final")
    
    return (Grid_final)

