#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:48:13 2025

@author: jleclezio
"""
import numpy as np

'''TRASNPORT A TRAVERS LA MEMBRANE GRADIENT OSMOTIQUE'''

'''flux osmotique, implicite'''

def transport_mb_osmo(SPECIES, species, Grid_1, Grid_2, nbmailles_x, delta_t, delta_x, dict_perm, Caract_space):
    """
    Étape de transport à travers la membrane via un gradient osmotique, de 1 vers 2

    Entrées :
    ----------
    - SPECIES : list
        Liste des noms des espèces étudiées
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Cgrad_x1, Cgrad_x1 : array
        Tableaux des concentrations dans espace1 et espace2 utilisé pour calculer les gradient
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - TF_x1, TF_x1 : array
        Tableaux des concentrations mises à jour après réactions et diffusion
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - dict_perm : dict
        Dictionnaire contenant les coefficients de perméabilité pour chaque espèce

    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après transport osmotique
        dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])

    Fonctionnement :
    ----------------
    1. Calcule le gradient osmotique (`grad_osmo`) entre espace1[N-1] et espace2[0] (frontiere gauche et frontiere droite de la mb)
    2. Met à jour les concentrations aux frontières avec la membrane en tenant compte des échanges osmotiques 
    3. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - Cgrad_x1, Cgrad_x2 : array
        Correspond aux concentrations au meme delta_t, apres reaction et diffusion
    """
    
    Cgrad_x1 = Grid_1
    Cgrad_x2 = Grid_2
    
    nbm_x1 = nbmailles_x[0]
    
    delta_x1 = delta_x[0]
    delta_x2 = delta_x[1]
    
    P = dict_perm.get(SPECIES[species])
    
    # gradient osmotique  : a verifier pour CG/dx)
    CG = Cgrad_x1[species][nbm_x1 - 1]
    CD = Cgrad_x2[species][0]

    grad_osmo = ((CG / delta_x1) - (CD / delta_x2))
    flux_osmo = P * grad_osmo
    
                   
    # exchange at the membrane : osmotique echange 
    # Euler application
    Grid_1[species][nbm_x1-1] = Grid_1[species][nbm_x1-1] - flux_osmo*delta_t
    Grid_2[species][0] = Grid_2[species][0] + flux_osmo*delta_t
    
    return (Grid_1, Grid_2)

'''flux osmotique, analytique'''

def transport_mb_osmo_analyt(SPECIES, species, Grid_1, Grid_2, nbmailles_x, dict_perm, delta_t, delta_x):
    """
    Étape de transport à travers la membrane via un gradient osmotique: 1 vers 2

    Entrées :
    ----------
    - SPECIES : list
        Liste des noms des espèces étudiées
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
   - Grid_1, Grid_2 : array
        Tableaux des concentrations mises à jour après réactions et diffusion pour notre espace temps
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - dict_perm : dict
        Dictionnaire contenant les coefficients de perméabilité pour chaque espèce
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    
    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après transport osmotique resolut sous forme analytique
        dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    
    Fonctionnement :
    ----------------
    1. Calcul la valeur propre a implementer dans la formule analytique
    2. Recupere les valeurs initiales
    3. Calcul le flux de transport
    4. Implement le flux dans les concentrations
    5. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - Cgrad_x1, Cgrad_x2 : array
        Correspond aux concentrations au meme delta_t, apres reaction et diffusion
    """
   
    Cgrad_x1 = Grid_1
    Cgrad_x2 = Grid_2
    
    nbm_x1 = nbmailles_x[0]
    nbm_x2 = nbmailles_x[1]
    
    #permeabilité
    P = dict_perm.get(SPECIES[species])
    
    # delta x 
    deltax1 = delta_x[0]
    deltax2 = delta_x[1]
    
    vp2 = - ( (P/deltax2) + (P/deltax1) )
    
    C0_x1 = Cgrad_x1[species][nbm_x1-1]
    C0_x2 = Cgrad_x2[species][0]
    
    # Calcul des constantes C1 et C2
    C1 = (C0_x1 + (C0_x2 * (deltax2 / deltax1)) ) / (1 + (deltax2 / deltax1))
    C2 = (C0_x2 - C0_x1) / (1 + (deltax2 / deltax1))
   
    # Calcul de CE_x1(t) et CE_x2(t)
    Ct_x1 = C1 - C2 * (deltax2 / deltax1) * np.exp(vp2*delta_t)
    Ct_x2 = C1 + C2 * np.exp(vp2*delta_t)
    
    # on a les concentrations qui sont transportés : on les implemente
    Grid_1[species][nbm_x1-1] = Ct_x1
    Grid_2[species][0] =  Ct_x2
    
    return (Grid_1, Grid_2)
       
'''flux electroosmotique, implicite'''

def PSI_fct(u):
    """
    Fonction regulatrice

    Entrées :
    ----------
    - u : int
    
    Sorties :
    ----------
    - res : int
        resultat de la fonction reg
        
    Fonctionnement : 
    ----------------
        deux choix: depend de la valeur de u, permet de ne pas diviser par 0 
    """

    if abs(u) > 0.01 :
        res = u/(np.exp(u)-1)

    else :
        # dvt limité de u/(e(u-1)) = 1 - (u/2) + ((u**2)/12)
        # schema de horner : 1 - (u/2) + ((u**2)/12) = 1  + u((-1/2) + u/12)
        res = 1  + u*((-1/2) + u/12) -(u**4)/720

    return (res)

def transport_mb_electro_osmo_impl(SPECIES, species, Grid_1, Grid_2, nbmailles_x, delta_t, delta_x, dict_perm, dict_charge):
    """
    Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014 : 1 vers 2 
    
    Entrées :
    ----------
    - SPECIES : list
        Liste des noms des espèces étudiées
    - species : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Cgrad_x1, Cgrad_x2 : array
        Tableaux des concentrations utilisé pour le calcul des gradiemt
    - TF_x1, TF_x1 : array
        Tableaux des concentrations mises à jour après réaction et diffusion pour notre espace temps
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - N, M : int
        Nombre de mailles dans espace1 (N) et espace2 (M)
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : float
        Pas d'espace pour la simulation
    - dict_perm : dict
        Dictionnaire contenant les coefficients de perméabilité pour chaque espèce
    - dict_charge : dict
        Dictionnaire contenant les charges pour chaque espèce
    
    Sorties :
    ----------
    - F_x1, F_x2 : array
        Tableaux des concentrations mises à jour après transport osmotique resolut sous forme analytique
        dans espace1 et espace2 respectivement
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    
    Fonctionnement :
    ----------------
    1. Calcul la variable zeta
    2. Appelle la fonction Psi(u)
    3. Recupere les concentrations à la membrane (`Xin`, `Xout`)
    $. Claule le flux de tranport (`J`)
    3. Implemente le flux dans les concentrations
    5. Retourne les tableaux mis à jour avec les nouvelles concentrations
    
    Notre simulation:
    ----------------
    - Cgrad_x1, Cgrad_x2 : array
        Correspond aux concentrations au meme delta_t, apres reaction et diffusion
    """
    
    Cgrad_x1 = Grid_1
    Cgrad_x2 = Grid_2
    
    # print("x1 = interieur, x2 = exterieur")
    P = dict_perm.get(SPECIES[species])
    Z = dict_charge.get(SPECIES[species])
    
    nbm_x1 = nbmailles_x[0]
    nbm_x2 = nbmailles_x[1]
    
    delta_x1 = delta_x[0]
    delta_x2 = delta_x[1]
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 #temperature en K : degrees Celsius + 273
    
    Vm = +70e-3 # A MODIFIER
    
    Zeta = (F*Vm)/(R*T)
    
    Psi = PSI_fct(Z*Zeta)
    
    Xint = (Cgrad_x1[species][nbm_x1-1])# /delta_x1
    Xout = (Cgrad_x2[species][0]) # /delta_x2
    
    # calcul le flux
    J = -P*Psi*(Xout - Xint * np.exp(Z*Zeta))
    
    # exchange at the membrane 
    Grid_1[species][nbm_x1-1] = Grid_1[species][nbm_x1-1] - (J*delta_t) # *delta_x1
    Grid_2[species][0] = Grid_2[species][0] + (J*delta_t) # *delta_x2

    return (Grid_1, Grid_2)


'''flux electro osmotique, analytique'''
# EN COURS

def transport_mb_electro_osmo_analy(SPECIES, species, Cgrad_x1, Cgrad_x2, TF_x1, TF_x2, N, M, delta_t, delta_x, dict_perm, dict_charge):
    """
    Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014
    """
    # print("x1 = interieur, x2 = exterieur")
    P = dict_perm.get(SPECIES[species])
    Z = dict_charge.get(SPECIES[species])
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 #temperature en K : degrees Celsius + 273
    
    Vm = -30e-3
    
    Zeta = (F*Vm)/(R*T)
    
    Psi = PSI_fct(Z*Zeta)
    
    E = np.exp(Z*Zeta)
    
    Xg0 = Cgrad_x1[species][N-1]
    Xd0 = Cgrad_x2[species][0]
    
    # calcul le flux
    
# full analytique : pas bon + trop de proba d'erreur de signe
# on veut : semi analytique
    Xg = (-1/(E-1)) * (-E*Xd0 + np.exp(Z*Zeta + ((E-1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + E*Xg0 - np.exp((E-1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)
    Xd = (1/(E+1)) * (Xd0 + np.exp(Z*Zeta + ((E+1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + Xg0 - np.exp((E+1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)

    
    # # exchange at the membrane 
    TF_x1[species][N-1] = Xg
    TF_x2[species][0] = Xd
   
    
    # if np.array_equal(TF_x1, TFx1_avt):
    #     print("cool x1")
    # else:
    #     print("pas egal x1")
    
    return (TF_x1, TF_x2)
