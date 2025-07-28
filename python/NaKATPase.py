#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 09:21:02 2025

@author: jleclezio

def upper & lower limite membrane potential
"""

import math 
import numpy as np
import matplotlib.pyplot as plt


X = (10**(-7.2))/(1.8e-8) 

C = (1.8e-8)/(3.6e-6)

L = 10**3

def sigma_NHE(x, L, C):
    
    num1 = x*(1+x) + L*C*x*(1+C*x)
    num2 = L*(1+C*x)**2 + (1+x)**2
    
    res = num1/num2
    
    return res


'''Initial condition : '''


Na_i = 10e-3
Na_e = 140e-3
z_Na = +1

K_i = 140e-3
K_e = 4e-3
z_K = +1

H_i = 6.31e-8  #pH = 7.2
H_e = 3.98e-8 #pH = 7.4
H_i_acid = 1e-6 #pH = 6
z_H = +1

# valeur fixe :
F = 96485 #cst de Faraday, J·V−1·mol−1
R = 8.314 #cst des gaz parfait, J·K−1·mol−1
T = 310 #temperature en K : degrees Celsius + 273
F_Nernst = (R*T)/F


"""Study on NaK ATPase"""

'''Membrane potential limits: '''

#upper for transporter
UL_Em = (math.log(Na_e/Na_i))*F_Nernst/z_Na
#lower for transporter
LL_Em_NaK = (math.log(K_e)-math.log(K_i))*F_Nernst/z_K # physiologic concentration (pH = 7.2)


# membrane potentiel Na:
Vm_Na = -((R*T)/(z_Na*F))*math.log(Na_i/Na_e)

# membrane potentiel K:
Vm_K = -((R*T)/(z_K*F))*math.log(K_i/K_e)


'''Permeability dependance : '''

## LIMITS :

#upper
UZeta = (F*UL_Em)/(R*T)
Upper_fact_P = (-2/3) * (Na_e-Na_i*math.exp(UZeta))/(K_e-K_i*math.exp(UZeta)) 

#lower
LZeta_NaK = (F*LL_Em_NaK)/(R*T)
Lower_fact_P_NaK = (-2/3) * (Na_e-Na_i*math.exp(LZeta_NaK))/(K_e-K_i*math.exp(LZeta_NaK))

## VARIATIONS :
    
# Em :
Em_values_NaK = np.linspace(LL_Em_NaK, UL_Em, 5000)

# Zeta fact_P :
Zeta_values_NaK = (F * Em_values_NaK) / (R * T)
fact_P_values_NaK = -(2/3) * (Na_e - Na_i * np.exp(Zeta_values_NaK)) / (K_e - K_i * np.exp(Zeta_values_NaK))


fact_P_values_NaK_cut = fact_P_values_NaK[1::]
Em_values_NaK_cut = Em_values_NaK[1::]



# value of Em for fact_p = 1 (tolerence : 1e-2)
POS_fact_P_1 = np.where(np.isclose(fact_P_values_NaK, 1.0, atol=1e-2))[0]
Em_values_NaK[POS_fact_P_1]


'''Plot factor permeability vs membrane potential'''

plt.figure()

plt.plot(Em_values_NaK * 1000, np.log10(fact_P_values_NaK))

plt.axhline(y=Upper_fact_P, color='red', linestyle=':', label='Upper_fact_P')
plt.axhline(y=Lower_fact_P_NaK, color='blue', linestyle=':', label='Lower_fact_P')

plt.xlim([-70,+70])
plt.ylim([-1,1])

plt.xlabel("Em (mV)")
plt.ylabel("fact_P")
plt.title("Value of log10 Fact_Perm according to membrane potential (Em) variation, for NaKATPasw")

plt.grid(True)
plt.axhline(y=0, color='grey', linestyle='--')

# plt.savefig("log_fact_P_fct_Em_NaKATPase.pdf")
plt.show()




"""Study on NHE"""

'''Membrane potential limits: '''
#upper for transporter
UL_Em = (math.log(Na_e/Na_i))*F_Nernst/z_Na
#lower for transporter
LL_Em_NHE = (math.log(H_e)-math.log(H_i))*F_Nernst/z_H # physiologic concentration (pH = 7.2)

LL_Em_NHE_ac = (math.log(H_e)-math.log(H_i_acid))*F_Nernst/z_H # acidic concentration (pH = 6)


# membrane potentiel Na:
Vm_Na = -((R*T)/(z_Na*F))*math.log(Na_i/Na_e)
# membrane potentiel H:
Vm_H = -((R*T)/(z_H*F))*math.log(H_i/H_e)


'''Permeability dependance : '''

## LIMITS :

#upper
UZeta = (F*UL_Em)/(R*T)
Upper_fact_P = (-2/3) * (Na_e-Na_i*math.exp(UZeta))/(K_e-K_i*math.exp(UZeta)) 

#lower
LZeta_NHE = (F*LL_Em_NHE)/(R*T)
Lower_fact_P_NHE = -(Na_e-Na_i*math.exp(LZeta_NHE))/(H_e-H_i*math.exp(LZeta_NHE))

LZeta_NHE_ac = (F*LL_Em_NHE_ac)/(R*T)
Lower_fact_P_NHE_ac = -(Na_e-Na_i*math.exp(LZeta_NHE_ac))/(H_e-H_i_acid*math.exp(LZeta_NHE_ac))

## VARIATIONS :
    
# Em :
Em_values_NHE = np.linspace(LL_Em_NHE, UL_Em, 5000)
Em_values_NHE_ac = np.linspace(LL_Em_NHE_ac, UL_Em, 5000)

# Zeta fact_P :
Zeta_values_NHE = (F * Em_values_NHE) / (R * T)
fact_P_values_NHE = - (Na_e - Na_i * np.exp(Zeta_values_NHE)) / (H_e - H_i * np.exp(Zeta_values_NHE))

Zeta_values_NHE_ac = (F * Em_values_NHE_ac) / (R * T)
fact_P_values_NHE_ac = - (Na_e - Na_i * np.exp(Zeta_values_NHE_ac)) / (H_e - H_i_acid * np.exp(Zeta_values_NHE_ac))


fact_P_values_NHE_cut = fact_P_values_NHE[1::]
Em_values_NHE_cut = Em_values_NHE[1::]

fact_P_values_NHE_ac_cut = fact_P_values_NHE_ac[1::]
Em_values_NHE_ac_cut = Em_values_NHE_ac[1::]


# value of Em for fact_p = 1 (tolerence : 1e-2)
# POS_fact_P_1 = np.where(np.isclose(fact_P_values_NHE, 1.0, atol=1e-2))[0]
# Em_values[POS_fact_P_1]

# POS_fact_P_1 = np.where(np.isclose(fact_P_values_NHE_ac, 1.0, atol=1e-2))[0]
# Em_values[POS_fact_P_1]

'''Plot factor permeability vs membrane potential'''

## NHE
plt.figure()

plt.plot(Em_values_NHE * 1000, np.log10(fact_P_values_NHE))

plt.axhline(y=Upper_fact_P, color='red', linestyle=':', label='Upper_fact_P')
# plt.axhline(y=Lower_fact_P, color='blue', linestyle=':', label='Lower_fact_P')

# plt.xlim([-70,+70])
# plt.ylim([-1,1])

plt.xlabel("Em (mV)")
plt.ylabel("fact_P")
plt.title("Value of log10 Fact_Perm according to membrane potential (Em) variation, NHE")

plt.grid(True)
plt.axhline(y=0, color='grey', linestyle='--')

plt.savefig("log_fact_P_fct_Em_NHE.pdf")
plt.show()


## NHE acidic
plt.figure()

plt.plot(Em_values_NHE_ac * 1000, np.log10(fact_P_values_NHE_ac))

plt.axhline(y=Upper_fact_P, color='red', linestyle=':', label='Upper_fact_P')
# plt.axhline(y=Lower_fact_P, color='blue', linestyle=':', label='Lower_fact_P')

# plt.xlim([-70,+70])
# plt.ylim([-1,1])

plt.xlabel("Em (mV)")
plt.ylabel("fact_P")
plt.title("Value of log10 Fact_Perm according to membrane potential (Em) variation, NHE for acidic cell")

plt.grid(True)
plt.axhline(y=0, color='grey', linestyle='--')

# plt.savefig("log_fact_P_fct_Em_NHE.pdf")
plt.show()





