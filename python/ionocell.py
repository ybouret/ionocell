""" IonoCell """
import numpy as np


class Grid:
    """ Grid for one space """

    def __init__(self, _space, _species_list, _transport_meth, _pas_temps, _pas_espace, _concentration ):
        
        self.space = _space #xi
        self.length = len(_space) #xi
        self.species = _species_list #SPECIES
        self.meth = _transport_meth #Type_mb
        self.dt = _pas_temps
        self.dx = _pas_espace
        self.value = _concentration #TS_x1

    def show(self):
        """ show Grid parameters """
        print("length=", self.length)
        print("species name =", self.species)
        
    # def add_time(self, t_value):
    #     self.times.append(t_value)
        
# gridx1 = Grid(x1, ["Na","Cl"], "osmose implicite", delta_t, delta_x,  TCF_x1)

class Specie:
    """ params for each species """
    
    def __init__(self, _specie_name, _diff_coeff, _perm_coeff, _charge, _marker):
        
        self.name = _specie_name,
        self.diff = _diff_coeff
        self.perm = _perm_coeff
        self.charge = _charge
        self.marker = _marker
    
    def show(self):
        print("specie name =" , self.name)


# test = Specie("Na", 1, 10, +1, '-d')

# Specie(_specie_name = "Na", _diff_coeff = 1, _perm_coeff = 10, _charge = +1, _marker= '-d')

class Reaction:
    """ params for reactions """
    
    def __init__(self, _TF, _reactif, _produit, _cst_eq):
        
        self.TF = _TF,
        self.reac = _reactif
        self.pdt = _produit
        self.cst_eq = _cst_eq
    
    def show(self):
        print("reaction ? " , self.TF)

class Space :
    
    def __init__(self, _number, _type, _two_wall, _one_inje):
        
        self.nb = _number
        self.type = _type
        self.walls = _two_wall
        self.inje = _one_inje
    
    def show(self):
        print("which space ?" , self.nb)
    




