#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:49:52 2025

@author: jleclezio
"""
import numpy as np
def diffusion_two_wall(specie, Grid_to_diff, nbmailles_x, delta_t, delta_x, v_coeffdiff):
    """
    Étape de diffusion des espèces à travers l'espace pour un pas de temps.

    Entrées :
    ----------
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Grid_to_diff_x1, Grid_to_diff_x2 : arrays
        Tableaux des concentrations dans l'espace au pas de temps précédent
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - v_coeffdiff : int
        Value coefficient de diffusion pour specie

    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après la diffusion dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Récupère le coefficient de diffusion de l'espèce étudiée
    2. Calcule le facteur alpha pour la diffusion (lié à delta_t, delta_x, et D)
    3. Met à jour les concentrations pour chaque maille intérieure (1 à N-2 pour 
       espace1, 1 à M-2 pour espace2) à l'aide de l'équation de diffusion
    4. Applique les conditions aux bords :
        aux extrémités & aux frontières avec la membrane: pas de flux (valeurs extrapolées)
    5. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - TSx1, TSx2 : array
        Correspond aux concentrations au pas de temps précédent avec la reaction d'implementer (si TFreaction = True)
    """      
    # print(np.shape(Grid_to_diff), "fct")
    D = v_coeffdiff
    
    alpha = (D * delta_t) / delta_x**2 

    nbm = nbmailles_x
    
    Gfin_sp = Grid_to_diff.copy()

    # diffusion inside
    for i in range(1, nbm-1): # do not take into account the boundary
        Gfin_sp[i] = Grid_to_diff[i] + alpha *(Grid_to_diff[i+1] - 2*Grid_to_diff[i] + Grid_to_diff[i-1])

    # no flux: boundary
    Gfin_sp[0] = (4*Grid_to_diff[1] - Grid_to_diff[2])/3  # x = 0 
    Gfin_sp[nbm-1] = (4*Grid_to_diff[nbm-2] - Grid_to_diff[nbm-3])/3 # x = 1
    
    # Gfin_sp[0] = Grid_to_diff[1] # à x = 0 
    # Gfin_sp[nbm-1] = Grid_to_diff[nbm-2] # à x = 1

    return (Gfin_sp)

def diffusion_one_wall_one_inject(specie, Grid_to_diff, nbmailles_x, delta_t, delta_x, v_coeffdiff):
    """
    Étape de diffusion des espèces à travers l'espace pour un pas de temps.

    Entrées :
    ----------
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Grid_to_diff_x1, Grid_to_diff_x2 : arrays
        Tableaux des concentrations dans l'espace au pas de temps précédent
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - v_coeffdiff : int
        Value coefficient de diffusion pour specie

    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après la diffusion dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Récupère le coefficient de diffusion de l'espèce étudiée
    2. Calcule le facteur alpha pour la diffusion (lié à delta_t, delta_x, et D)
    3. Met à jour les concentrations pour chaque maille intérieure (1 à N-2 pour 
       espace1, 1 à M-2 pour espace2) à l'aide de l'équation de diffusion
    4. Applique les conditions aux bords :
        aux extrémités & aux frontières avec la membrane: pas de flux (valeurs extrapolées)
    5. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - TSx1, TSx2 : array
        Correspond aux concentrations au pas de temps précédent avec la reaction d'implementer (si TFreaction = True)
    """      
    
    D = v_coeffdiff
    
    alpha = (D * delta_t) / delta_x**2 
        
    nbm = nbmailles_x

    # diffusion inside
    for i in range(1, nbm-1): # do not take into account the boundary
        Grid_to_diff[specie][i] = Grid_to_diff[specie][i] + alpha *(Grid_to_diff[specie][i+1] - 2*Grid_to_diff[specie][i] + Grid_to_diff[specie][i-1])

    # no flux: boundary
    Grid_to_diff[specie][0] = 1
    Grid_to_diff[specie][nbm-1] = (4*Grid_to_diff[specie][nbm-2] - Grid_to_diff[specie][nbm-3])/3 # à x = 1

    Grid_diff = Grid_to_diff
    
    return (Grid_diff)

def diffusion_one_wall_one_supp(specie, Grid_to_diff, nbmailles_x, delta_t, delta_x, v_coeffdiff):
    """
    Étape de diffusion des espèces à travers l'espace pour un pas de temps.

    Entrées :
    ----------
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Grid_to_diff_x1, Grid_to_diff_x2 : arrays
        Tableaux des concentrations dans l'espace au pas de temps précédent
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - v_coeffdiff : int
        Value coefficient de diffusion pour specie

    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après la diffusion dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Récupère le coefficient de diffusion de l'espèce étudiée
    2. Calcule le facteur alpha pour la diffusion (lié à delta_t, delta_x, et D)
    3. Met à jour les concentrations pour chaque maille intérieure (1 à N-2 pour 
       espace1, 1 à M-2 pour espace2) à l'aide de l'équation de diffusion
    4. Applique les conditions aux bords :
        aux extrémités & aux frontières avec la membrane: pas de flux (valeurs extrapolées)
    5. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - TSx1, TSx2 : array
        Correspond aux concentrations au pas de temps précédent avec la reaction d'implementer (si TFreaction = True)
    """      
    
    D = v_coeffdiff
    
    alpha = (D * delta_t) / delta_x**2 
        
    nbm = nbmailles_x

    # diffusion inside
    for i in range(1, nbm-1): # do not take into account the boundary
        Grid_to_diff[specie][i] = Grid_to_diff[specie][i] + alpha *(Grid_to_diff[specie][i+1] - 2*Grid_to_diff[specie][i] + Grid_to_diff[specie][i-1])

    # no flux: boundary
    Grid_to_diff[specie][0] = 0
    Grid_to_diff[specie][nbm-1] = (4*Grid_to_diff[specie][nbm-2] - Grid_to_diff[specie][nbm-3])/3 # à x = 1

    Grid_diff = Grid_to_diff
    
    return (Grid_diff)


