#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:48:13 2025

@author: jleclezio
"""
import math 
import numpy as np

'''aucun flux'''

def non_transport(Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge, Em):
    
    return (Grid_L)
    
'''flux osmotique, explicite'''

def transport_mb_osmo(Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge, Em):
    """
    Étape de transport à travers la membrane via un gradient osmotique

    Entrées :
    ----------
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
    
    nbm_x1 = nbmailles_x[px1]
    
    delta_x1 = delta_x[px1]
    delta_x2 = delta_x[px2]
    
    P = v_perm
    
    # gradient osmotique  
    CG = Grid_L[px1][nbm_x1 - 1].copy() # gauche de la membrane
    CD = Grid_L[px2][0].copy() # droite de la membrane

    grad_osmo = (CG - CD)
    flux_osmo = P * grad_osmo
    
    Grid_L[px1][nbm_x1-1] = Grid_L[px1][nbm_x1-1] - (flux_osmo * delta_t)/delta_x2
    Grid_L[px2][0] = Grid_L[px2][0] + (flux_osmo * delta_t)/delta_x1
    
    return (Grid_L)


def flux_mb_osmo(Species, C_int, C_ext, Mb_1, nbmailles_x, delta_t, delta_x):
        
        #Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge, Em):
   
    px1 = Mb_1[0]
    px2 = Mb_1[1]
    
    P = Species.perm
    
    # gradient osmotique  
    grad_osmo = (C_int - C_ext)
    flux_osmo = P * grad_osmo
    
    return (flux_osmo)


'''flux osmotique, analytique'''

def transport_mb_osmo_analyt(Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge, Em):
    """
    Étape de transport à travers la membrane via un gradient osmotique
    
    Entrées :
    ----------
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
    px1 = Mb_1[0] # espace gauche mb
    px2 = Mb_1[1] # espace droite mb
    
    Cgrad_x1 = Grid_L[px1]
    Cgrad_x2 = Grid_L[px2]
    
    nbm_x1 = nbmailles_x[px1]

    deltax1 = delta_x[px2]
    deltax2 = delta_x[px1]
    
    #permeabilité
    P = v_perm
    
    vp2 = - ( (P/deltax1) + (P/deltax2) )
    
    C0_x1 = Cgrad_x1[nbm_x1-1] # a gauceh de la mb
    C0_x2 = Cgrad_x2[0] # a droite de la mb
    
    # Calcul des constantes C1 et C2
    C1 = (C0_x1 + (C0_x2 * (deltax2 / deltax1)) ) / (1 + (deltax2 / deltax1))
    C2 = (C0_x2 - C0_x1) / (1 + (deltax2 / deltax1))
   
    # Calcul de CE_x1(t) et CE_x2(t)
    Ct_x1 = C1 - C2 * (deltax2 / deltax1) * math.exp(vp2*t)
    Ct_x2 = C1 + C2 * math.exp(vp2*t)
    
    # on a les concentrations qui sont transportées : on les implemente
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

def flux_electro_osmo_sp(Species, C_int, C_ext, Mb_1, Em, nbmailles_x, delta_t, delta_x):
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 # K  =  37 degrees C 
    
    Zeta = (F*Em)/(R*T) # dimensionless
    
    P = Species.perm
    Z = Species.charge
    
    Psi = PSI_fct(Z*Zeta)
    
    # calcul le flux
    J = -P*Psi*(C_ext - (C_int * math.exp(Z*Zeta)))

    return (J)
    
    
    

def transport_mb_electro_osmo_impl(Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, v_perm, v_charge, Em):
    """
    Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014 : 1 vers 2 
    
    Entrées :
    ----------
    - specie : int
        Index de l'espèce étudiée dans la liste SPECIES
    - Cgrad_x1, Cgrad_x2 : array
        Tableaux des concentrations utilisé pour le calcul des gradients
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
    
    Cgrad_x1 = Grid_L[px1].copy()
    Cgrad_x2 = Grid_L[px2].copy()
    
    nbm_x1 = nbmailles_x[px1]

    deltax1 = delta_x[px2]
    deltax2 = delta_x[px1]
    
    # print("x1 = interieur, x2 = exterieur")
    P = v_perm
    Z = v_charge
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 # K  =  37 degrees C 
    
    Vm = Em  # V
    
    Zeta = (F*Vm)/(R*T) # dimensionless
    
    Psi = PSI_fct(Z*Zeta)
    
    Xint = (Cgrad_x1[nbm_x1-1])
    Xout = (Cgrad_x2[0])
    
    # calcul le flux
    J = -P*Psi*(Xout - (Xint * math.exp(Z*Zeta)))
    
    # exchange at the membrane 
    Grid_1 = Grid_L[px1].copy()
    Grid_2 = Grid_L[px2].copy()
    
    Grid_1[nbm_x1-1] = Grid_1[nbm_x1-1] - (J*delta_t)/deltax1
    Grid_2[0] = Grid_2[0] + (J*delta_t)/deltax2  
    
    Grid_L[px1] = Grid_1
    Grid_L[px2] = Grid_2

    return (Grid_L)


def transport_mb_electro_osmo_NaK(SPECIES, Grid_L, t, Mb_1, nbmailles_x, delta_t, delta_x, Em, Rho_NaK):

    px1 = Mb_1[0]
    px2 = Mb_1[1]
    nbm_x1 = nbmailles_x[px1]

    deltax1 = delta_x[px2]
    deltax2 = delta_x[px1]
    
    # valeur fixe :
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 # K  =  37 degrees C 
    Vm = Em  # V
    
    Zeta = (F*Vm)/(R*T) # dimensionless
    
    
    sp_num = 0
    for sp in SPECIES:
        name_sp = sp.name[0]
        v_perm = sp.perm
        v_charge = sp.charge
        
        Grid_to_transp_sp = [arr[sp_num] for arr in Grid_L]
        
        if name_sp == "K":
            
            K_pos = sp_num
            
            # def J_elect_osmo_K :
            Cgrad_x1 = Grid_to_transp_sp[px1].copy()
            Cgrad_x2 = Grid_to_transp_sp[px2].copy()
            
            P = v_perm
            Z = v_charge
            
            Psi = PSI_fct(Z*Zeta)
            
            Xint = (Cgrad_x1[nbm_x1-1])
            Xout = (Cgrad_x2[0])
            
            # K_e = Xout
            
            # flux
            J_elect_osmo_K = -P*Psi*(Xout - (Xint * math.exp(Z*Zeta)))
            print(J_elect_osmo_K, "J_elect_osmo_K")
            
        elif name_sp == "Na":
            
            Na_pos = sp_num
            
            # def J_elect_osmo_Na :
            Cgrad_x1 = Grid_to_transp_sp[px1].copy()
            Cgrad_x2 = Grid_to_transp_sp[px2].copy()

            P = v_perm
            Z = v_charge
            
            Psi = PSI_fct(Z*Zeta)
            
            Xint = (Cgrad_x1[nbm_x1-1])
            Xout = (Cgrad_x2[0])
            
            # Na_i = Xint
            
            # flux
            J_elect_osmo_Na = -P*Psi*(Xout - (Xint * math.exp(Z*Zeta)))
            print(J_elect_osmo_Na, "J_elect_osmo_Na")
            
        sp_num = sp_num+1
    
    
    # Vmax_NaK = 1.64
    # K_Nai = 10e-3
    # K_Ke = 1e-3
    
    # Na_sat = (Na_i)/(K_Nai + Na_i)
    # K_sat = (K_e)/(K_Ke + K_e)
    
    # Rho_NaK = Vmax_NaK * Na_sat**3 * K_sat**2
    
    J_tot_K = J_elect_osmo_K - 2 * Rho_NaK
    # print(J_tot_K, "J_tot_K")
    
    J_tot_Na = J_elect_osmo_Na + 3 * Rho_NaK
    # print(J_tot_Na, "J_tot_Na") 
    
    
    # pour K : 
    Grid_K1 = Grid_L[px1][K_pos]
    Grid_K2 = Grid_L[px2][K_pos]
    
    # Grid_K1[nbm_x1-1] = Grid_K1[nbm_x1-1] - (J_tot_K*delta_t)/deltax1 
    Grid_K2[0] = Grid_K2[0] + (J_tot_K*delta_t)/deltax2
    
    Grid_L[px1][K_pos] = Grid_K1
    Grid_L[px2][K_pos] = Grid_K2
    
    # pour Na : 
    Grid_Na1 = Grid_L[px1][Na_pos]
    Grid_Na2 = Grid_L[px2][Na_pos]
    
    Grid_Na1[nbm_x1-1] = Grid_Na1[nbm_x1-1] - (J_tot_Na*delta_t)/deltax1 
    Grid_Na2[0] = Grid_Na2[0] + (J_tot_Na*delta_t)/deltax2
    
    Grid_L[px1][Na_pos] = Grid_Na1
    Grid_L[px2][Na_pos] = Grid_Na2
    
    return Grid_L

# '''flux electro osmotique, analytique'''
# # EN COURS

# def transport_mb_electro_osmo_analy(SPECIES, specie, Cgrad_x1, Cgrad_x2, TF_x1, TF_x2, N, M, delta_t, delta_x, v_perm, v_charge):
#     """
#     Étape de transport à travers la membrane avec equation GHK, recuperation du modele de 2014
#     """
#     # print("x1 = interieur, x2 = exterieur")
#     P = v_perm
#     Z = v_charge
    
#     # valeur fixe :
#     F = 96485 #cst de Faraday, J·V−1·mol−1
#     R = 8.314 #cst des gaz parfait, J·K−1·mol−1
#     T = 310 #temperature en K : degrees Celsius + 273
    
#     Vm = -30e-3
    
#     Zeta = (F*Vm)/(R*T)
    
#     Psi = PSI_fct(Z*Zeta)
    
#     E = np.exp(Z*Zeta)
    
#     Xg0 = Cgrad_x1[specie][N-1]
#     Xd0 = Cgrad_x2[specie][0]
    
#     # calcul le flux
    
# # full analytique : pas bon + trop de proba d'erreur de signe
# # on veut : semi analytique
#     Xg = (-1/(E-1)) * (-E*Xd0 + np.exp(Z*Zeta + ((E-1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + E*Xg0 - np.exp((E-1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)
#     Xd = (1/(E+1)) * (Xd0 + np.exp(Z*Zeta + ((E+1)*P*delta_t*Zeta*Psi/delta_x))*Xd0 + Xg0 - np.exp((E+1)*P*delta_t*Zeta*Psi/delta_x)*Xg0)

    
#     # # exchange at the membrane 
#     TF_x1[specie][N-1] = Xg
#     TF_x2[specie][0] = Xd
   
    
#     # if np.array_equal(TF_x1, TFx1_avt):
#     #     print("cool x1")
#     # else:
#     #     print("pas egal x1")
    
#     return (TF_x1, TF_x2)
