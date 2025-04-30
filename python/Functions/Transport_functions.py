#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:48:13 2025

@author: jleclezio
"""
import numpy as np

'''TRASNPORT A TRAVERS LA MEMBRANE GRADIENT OSMOTIQUE'''

'''flux osmotique, implicite'''

def transport_mb_osmo(list_SP, Grid_L, Mb_1, nbmailles_x, delta_t, delta_x, v_perm):
    """
    Étape de transport à travers la membrane via un gradient osmotique, de 1 vers 2

    Entrées :
    ----------
    - list_SP : list
        Liste des espèces étudiées, en Class: SPECIE
    - specie : int
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
    - v_perm : int
        Value coefficient de perméabilité 

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
    px1 = Mb_1[0]
    px2 = Mb_1[1]
    
    # Cgrad_x1 = Grid_L[px1]
    # Cgrad_x2 = Grid_L[px2]
    
    nbm_x1 = nbmailles_x[px1]
    
    delta_x1 = delta_x[px1]
    delta_x2 = delta_x[px2]
    
    P = v_perm
    
    # gradient osmotique  : a verifier pour CG/dx)
    CG = Grid_L[px1][nbm_x1 - 1]
    CD = Grid_L[px2][0]

    # # gradient osmotique  : a verifier pour CG/dx)
    # CG = Cgrad_x1[nbm_x1 - 1]
    # CD = Cgrad_x2[0]

    grad_osmo = (CG - CD)
    flux_osmo = P * grad_osmo
    
    # exchange at the membrane : osmotique echange 
    # Euler application
    # Grid_1 = Grid_L[px1]
    # Grid_2 = Grid_L[px2]
    
    Grid_L[px1][nbm_x1-1] = Grid_L[px1][nbm_x1-1] - (flux_osmo/delta_x2)*delta_t # prise en compte du delta_x
    Grid_L[px2][0] = Grid_L[px2][0] + (flux_osmo/delta_x1)*delta_t
    
    
    
    # Grid_1[nbm_x1-1] = Grid_1[nbm_x1-1] - (flux_osmo/delta_x1)*delta_t # prise en compte du delta_x
    # Grid_2[0] = Grid_2[0] + (flux_osmo/delta_x2)*delta_t
    
    # Grid_L[px1] = Grid_1
    # Grid_L[px2] = Grid_2
    return (Grid_L)

'''flux osmotique, analytique'''

def transport_mb_osmo_analyt(list_SP, Grid_L, Mb_1, nbmailles_x, v_perm, delta_t, delta_x):
    """
    Étape de transport à travers la membrane via un gradient osmotique: 1 vers 2

    Entrées :
    ----------
    - list_SP : list
        Liste des espèces étudiéesm en Class: Specie
    - specie : int
        Index de l'espèce étudiée dans la liste SPECIES
   - Grid_1, Grid_2 : array
        Tableaux des concentrations mises à jour après réactions et diffusion pour notre espace temps
        Format : array([nombre d'espèces, nombre de mailles d'espace])
    - v_perm : int
        Value coefficient de perméabilité 
    - delta_t : float
        Pas de temps pour la simulation
    - delta_x : list
        Pas d'espace pour la simulation :[dx1, dx2]
    
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
    px1 = Mb_1[0]
    px2 = Mb_1[1]
    
    Cgrad_x1 = Grid_L[px1]
    Cgrad_x2 = Grid_L[px2]
    
    nbm_x1 = nbmailles_x[px1]
    
    deltax1 = delta_x[px1]
    deltax2 = delta_x[px2]
    
    #permeabilité
    P = v_perm
    
    vp2 = - ( (P/deltax1) + (P/deltax2) )
    
    C0_x1 = Cgrad_x1[nbm_x1-1]
    C0_x2 = Cgrad_x2[0]
    
    # Calcul des constantes C1 et C2
    C1 = (C0_x1 + (C0_x2 * (deltax2 / deltax1)) ) / (1 + (deltax2 / deltax1))
    C2 = (C0_x2 - C0_x1) / (1 + (deltax2 / deltax1))
   
    # Calcul de CE_x1(t) et CE_x2(t)
    Ct_x1 = C1 - C2 * (deltax2 / deltax1) * np.exp(vp2*delta_t)
    Ct_x2 = C1 + C2 * np.exp(vp2*delta_t)
    
    # on a les concentrations qui sont transportés : on les implemente
    Grid_1 = Grid_L[px1]
    Grid_2 = Grid_L[px2]
    
    Grid_1[nbm_x1-1] = Ct_x1
    Grid_2[0] =  Ct_x2
    
    Grid_L[px1] = Grid_1
    Grid_L[px2] = Grid_2
    
    return (Grid_L)
       
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

def transport_mb_electro_osmo_impl(SPECIES, Grid_L, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge):
    # delta x a implementer 
    """
    Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014 : 1 vers 2 
    
    Entrées :
    ----------
    - SPECIES : list
        Liste des noms des espèces étudiées
    - specie : int
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
    - v_perm : int
        Value coefficient de perméabilité
    - v_charge : int
        Value charge
    
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
    
    px1 = Mb_1[0]
    px2 = Mb_1[1]
    
    Cgrad_x1 = Grid_L[px1]
    Cgrad_x2 = Grid_L[px2]
    
    nbm_x1 = nbmailles_x[px1]
    
    deltax1 = delta_x[px1]
    deltax2 = delta_x[px2]
    
    # print("x1 = interieur, x2 = exterieur")
    P = v_perm
    Z = v_charge
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 #temperature en K : degrees Celsius + 273
    
    Vm = +70e-3 # A MODIFIER
    
    Zeta = (F*Vm)/(R*T)
    
    Psi = PSI_fct(Z*Zeta)
    
    Xint = (Cgrad_x1[nbm_x1-1])# /delta_x1
    Xout = (Cgrad_x2[0]) # /delta_x2
    
    # calcul le flux
    J = -P*Psi*(Xout - Xint * np.exp(Z*Zeta))
    
    # exchange at the membrane 
    Grid_1 = Grid_L[px1]
    Grid_2 = Grid_L[px2]
    
    Grid_1[nbm_x1-1] = Grid_1[nbm_x1-1] - (J*delta_t) # *delta_x1
    Grid_2[0] = Grid_2[0] + (J*delta_t) # *delta_x2
    
    Grid_L[px1] = Grid_1
    Grid_L[px2] = Grid_2

    return (Grid_L)


'''flux electro osmotique, analytique'''
# EN COURS

def transport_mb_electro_osmo_analy(SPECIES, specie, Cgrad_x1, Cgrad_x2, TF_x1, TF_x2, N, M, delta_t, delta_x, v_perm, v_charge):
    """
    Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014
    """
    # print("x1 = interieur, x2 = exterieur")
    P = v_perm
    Z = v_charge
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 #temperature en K : degrees Celsius + 273
    
    Vm = -30e-3
    
    Zeta = (F*Vm)/(R*T)
    
    Psi = PSI_fct(Z*Zeta)
    
    E = np.exp(Z*Zeta)
    
    Xg0 = Cgrad_x1[specie][N-1]
    Xd0 = Cgrad_x2[specie][0]
    
    # calcul le flux
    
# full analytique : pas bon + trop de proba d'erreur de signe
# on veut : semi analytique
    Xg = (-1/(E-1)) * (-E*Xd0 + np.exp(Z*Zeta + ((E-1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + E*Xg0 - np.exp((E-1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)
    Xd = (1/(E+1)) * (Xd0 + np.exp(Z*Zeta + ((E+1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + Xg0 - np.exp((E+1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)

    
    # # exchange at the membrane 
    TF_x1[specie][N-1] = Xg
    TF_x2[specie][0] = Xd
   
    
    # if np.array_equal(TF_x1, TFx1_avt):
    #     print("cool x1")
    # else:
    #     print("pas egal x1")
    
    return (TF_x1, TF_x2)
