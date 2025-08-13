#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 11:17:01 2025

@author: jleclezio
"""
from Initiation_functions import condition_initiale
from Transport_functions import PSI_fct
import math

import numpy as np

def permeability_from_Na(Species_L, dict_CI, list_space, nbmailles, dict_Mb):
    
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 # K  =  37 degrees C 
    Em = dict_Mb["potential"]
    
    Zeta_values_NaK = (F*Em)/(R*T)
    
    for i in range(1, dict_Mb["nb_mb"]+1): # pour chaque mb
        key = f"Mb{i}"
        Mb_pos = dict_Mb[key]

        vars_dict = {} 
        for species in Species_L:
            name_sp = species.name[0] 

            key_i = f"{name_sp}_x{Mb_pos[0] + 1}"
            key_e = f"{name_sp}_x{Mb_pos[1] + 1}"
            
            i_array = dict_CI.get(key_i)
            e_array = dict_CI.get(key_e)
            
            if i_array is not None and e_array is not None:
                vars_dict[f"{name_sp}_i"] = i_array[-1]
                vars_dict[f"{name_sp}_e"] = e_array[0]
            else:
                print(f"Clé manquante pour {name_sp} : {key_i} ou {key_e}")
            
            
            if name_sp == "Na" :
                P_Na = species.perm

        nominateur = (vars_dict["Na_e"] - vars_dict["Na_i"] * math.exp(Zeta_values_NaK))
        denominateur = (vars_dict["K_e"] - vars_dict["K_i"] * math.exp(Zeta_values_NaK))

        fact_P = -(2/3) * nominateur / denominateur
    
    P_K = fact_P*P_Na
    
    return P_K

def Rho_from_Flux(Species_L, dict_CI, list_space, nbmailles, dict_Mb):
    
    # calcul de flux electro osmotique : 
        
    F = 96485 #cst de Faraday, J·V−1·mol−1
    R = 8.314 #cst des gaz parfait, J·K−1·mol−1
    T = 310 # K  =  37 degrees C 
    Em = dict_Mb["potential"]
    
    Zeta = (F*Em)/(R*T)
    
    for i in range(1, dict_Mb["nb_mb"]+1): # pour chaque mb
        key = f"Mb{i}"
        Mb_pos = dict_Mb[key]

        vars_dict = {} 
        for species in Species_L:
            
            name_sp = species.name[0]
                
            SP_i = f"{name_sp}_x{Mb_pos[0] + 1}"
            SP_e = f"{name_sp}_x{Mb_pos[1] + 1}"
            
            C_i = dict_CI.get(SP_i)[nbmailles[0]-1]
            C_e = dict_CI.get(SP_e)[0]
                
            P = species.perm
            Z = species.charge
                
                
            Psi = PSI_fct(Z*Zeta)
                
            # flux
            J_elect_osmo = -P*Psi*(C_e - (C_i * math.exp(Z*Zeta)))
                
            vars_dict[f"{name_sp}"] = J_elect_osmo
    
    # calcul de rho:
    
    for sp in vars_dict.keys():
        if sp == "Na":
            Rho_Na = (1/3)*(-vars_dict[sp])
        
        elif sp == "K":
            Rho_K = (-1/2)*(-vars_dict[sp])
    
    if (Rho_Na - Rho_K) < 1e-10 :
        Rho_NaK = Rho_Na
        print(Rho_NaK, "Rho_NaK")
        return Rho_NaK
    
    
    
    else :
        print("calcul Rho_NaK mauvais.")
        Rho_NaK = None
        return Rho_NaK
