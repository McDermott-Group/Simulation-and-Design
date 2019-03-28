from scipy.constants import mu_0, epsilon_0, pi, hbar, e, c
from scipy.constants import physical_constants
from scipy.optimize import fsolve, root
import numpy as np

class Qubit(object):
    name = None
    omega_q = None
    omega_r = None
    C_r = None
    C_g = None
    C_q = None
    
    def __init__(self, name=None):
        if name:
            self.name = name
    
    def EkdivEc(self, ng=[0.5], Ec=1.0/4.0, Ej=2.0, N=50):
        # from https://github.com/thomasaref/TA_software/blob/3fea8a37d86239353dccd7a1b71dcedb41b7e64c/test_code/IDTQubitdesignpy.py#L62
        # calculates transmon energy level with N states (more states is better approximation)
        # effectively solves the mathieu equation but for fractional inputs (which doesn't work in scipy.special.mathieu_a
        d1=[]
        d2=[]
        d3=[]
        for a in ng:
            NL=2*N+1
            A=np.zeros((NL, NL))
            for b in range(0,NL):
                A[b, b]=4.0*Ec*(b-N-a)**2
                if b!=NL-1:
                    A[b, b+1]= -Ej/2.0
                if b!=0:
                    A[b, b-1]= -Ej/2.0
            w,v=np.linalg.eig(A)
            d1.append(min(w))#/h*1e-9)
            w=np.delete(w, w.argmin())
            d2.append(min(w))#/h*1e-9)
            w=np.delete(w, w.argmin())
            d3.append(min(w))#/h*1e-9)
        return np.array(d1), np.array(d2), np.array(d3)
    
    def E_c(self, c=None):
        if c is None:
            c = self.C_q + self.C_g
        return e**2/2/c

    def cap_g(self, g):
        return -((self.C_q+self.C_r)*(2*g)**2 + np.sqrt((self.C_q-self.C_r)**2*(2*g)**4 + 4*self.C_q*self.C_r*(2*g)**2*self.omega_q*self.omega_r)) / (2*((2*g)**2-self.omega_q*self.omega_r))

    def g(self):
        return 1./2*self.C_g/np.sqrt((self.C_q+self.C_g)*(self.C_r+self.C_g))*np.sqrt(self.omega_r*self.omega_q)
    
    def Chi_0(self):
        g = self.g()
        return -g**2/(self.omega_q-self.omega_r)
    
    def Chi(self):
        alpha = self.alpha(self.E_c(), self.E_j())
        Delta = self.omega_q-self.omega_r
        return -self.g()**2/Delta*(-alpha)/(Delta+alpha)
    
    def Q_r(self):
        return -self.omega_r/(2*self.Chi())
    
    def E01(self, E_c, E_j, ng=[0.0]):
        d1, d2, d3 = self.EkdivEc(ng, E_c, E_j)
        return d2[0]-d1[0]
    
    def alpha(self, E_c, E_j, ng=[0.0]):
        d1, d2, d3 = self.EkdivEc(ng, E_c, E_j)
        E12 = d3[0]-d2[0]
        E01 = d2[0]-d1[0]
        return (E12 - E01)/hbar
    
    def E_j(self):
        f = lambda ej: self.E01(self.E_c(), ej) - self.omega_q*hbar
        return fsolve(f, 2*pi*15e9*hbar)[0]
    
    def I_c(self):
        return 2*pi*self.E_j()/(physical_constants['mag. flux quantum'][0])