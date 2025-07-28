#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:52:28 2025

@author: jleclezio
"""

import numpy as np


'''REACTION EULER'''

def reaction_euler(list_SP, Grid_to_react, react_dict):
    """
    Calcule les concentrations dans l'espace apres une reaction

    Entrées :
    ----------
    - list_SP : list
        Liste des espèces étudiées, en Class : Specie
    - Grid_to_react : array
        Tableaux contenant les concentrations dans espace1 et espace2 au pas de temps précédent
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - react_dict : dict
        Dictionnaire décrivant les réactions chimiques impliquant les espèces
    - t : int
        temps de la simulation

    Sorties :
    ----------
    - TSx1, TSx1 : array
        Tableaux des concentrations calculées pour chaque espèce dans espace1 et espace2 respectivement apres les reactions
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    Methode : Euler
   
    """
    nb_react = len(react_dict)
    
    NS = len(list_SP)
    
    Grid_to_react_t = [np.zeros_like(arr) for arr in Grid_to_react]
    
    SPECIES = []
    for sp in list_SP :
        SPECIES.append(sp.name[0])

    for R in range (0, nb_react):
        reaction0 = react_dict[f"reaction{R}"]
        k_react = reaction0['cst_eq']
        
        REACTIF = []
        for reactif in reaction0['reactif']:
            REACTIF.append(Grid_to_react[SPECIES.index(reactif)])
            
        REACTIFS = np.prod(REACTIF, axis=0)
        
        for species in range(0, NS):
            sp = SPECIES[species]
            
            if sp in reaction0["reactif"]:
                Grid_to_react_t[species] = Grid_to_react[species] - k_react*REACTIFS
            elif sp in reaction0["produit"]:
                Grid_to_react_t[species] = Grid_to_react[species] + k_react*REACTIFS

    return (Grid_to_react_t)



