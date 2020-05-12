"""
Parameters relevant for the aerodynamics department. Update in case necessary

@author: Casper, Ihab, Eduardo

"""
# =============================================================================
# Import relevant modules 
import numpy as np 

# =============================================================================
# Main 
class Wing_A(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['rootchord', 'tipchord', 'span','sweepdegLErad','wingloading','taper','MTOM','surfacearea']
    def __init__(self):
        self.rootchord = 0.76
        self.tipchord = 0.113
        self.span = 3
        self.sweepdegLErad = 20*np.pi/180
        self.wingloading = 122.26
        self.taper = self.tipchord/self.rootchord
        self.MTOM = 16.217
        self.surfacearea = self.wingloading/(self.MTOM*9.81)
        
class Concept_A(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['MTOM', 'Vcruise', 'Vtransition', 'htransition']
    def __init__(self):
        self.MTOM = 16.217
        self.Vcruise = 28
        self.Vtransition = 14
        self.htransition = 30
        
class Wing_D(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['rootchord', 'tipchord', 'span','sweepdegLErad','wingloading','taper','MTOM','surfacearea']
    def __init__(self):
        self.rootchord = 0.76
        self.tipchord = 0.113
        self.span = 3
        self.sweepdegLErad = 20*np.pi/180
        self.wingloading = 122.26
        self.taper = self.tipchord/self.rootchord
        self.MTOM = 16.217
        self.surfacearea = self.wingloading/(self.MTOM*9.81)
    
class Concept_D(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['MTOM', 'Vcruise', 'Vtransition', 'htransition']
    def __init__(self):
        self.MTOM = 16.217
        self.Vcruise = 28
        self.Vtransition = 14
        self.htransition = 30

        
class Wing_C(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['wingloading', 'span', 'MTOM', 'surfaceaera','averagechord','sweepdegLErad']
    def __init__(self):
        self.wingloading = 122.26
        self.span = 3
        self.MTOM = 16.217
        self.surfaceaera = self.wingloading/(self.MTOM*9.81)
        self.averagechord = self.surfaceaera/self.span
        self.sweepdegLErad = 0 
        
class Concept_C(object):
     #NOTE use __slots__ to avoid unintended mistakes
    __slots__ = ['MTOM', 'Vcruise', 'Vtransition','htransition','taper','hcruise']
    def __init__(self):
        self.MTOM = 16.217
        self.Vcruise = 28
        self.Vtransition = 14
        self.htransition = 30
        self.hcruise = 500
        
class Parameters(object):
    __slots__ = ['Wing_A','Concept_A','Wing_C','Concept_C','Wing_D','Concept_D']
    def __init__(self):
        self.Wing_A = Wing_A()
        self.Concept_A = Concept_A()
        self.Wing_C = Wing_C()
        self.Concept_C = Concept_C()
        self.Wing_D = Wing_D()
        self.Concept_D = Concept_D()
        
parameters = Parameters()
        

        
     
    
    
