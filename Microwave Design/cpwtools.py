"""
Chris Wilen
McDermott Group, UW Madison

This file provides objects for a variety of different types of cpw resonators.
"""

from scipy.constants import c, epsilon_0, mu_0, h, hbar, e, pi
from scipy.special import ellipk, ellipkm1
import numpy as np
import capacitance as cap

def dBm_to_watts(dB):
    return 10**(dB/10.) * 0.001

def watts_to_dBm(watts):
    return np.log10( 1000 * watts ) * 10

class CPW(object):
    """
    w = width of center trace [um]
    s = width of gap [um]
    t = thickness [um]
    h = height of substrate [um]
    l = length of CPW [um]
    e1 = relative permitivity of substrate
    material = ['nb','al'], assigns Tc and rho
    tand =
    physics from goppl unless otherwise noted.
    """
    def __init__(self, w=10., s=6., t=.1, h=380., e1=11.8, material="nb", tand=.15e-6):
        self.w = w*1e-6
        self.s = s*1e-6
        self.t = t*1e-6
        self.h = h*1e-6
        self.e0 = 1.  # air above CPW
        self.e1 = e1  # dielectric
        self.tand = tand
        self.material = material
        self.Delta = (1.25*self.t/pi)*( 1 + np.log(4*pi*self.w/self.t) ) # correction for thickness of metal
        self.s_eff = self.s - self.Delta  # &&& make these into functions so if they are changed, everything works out nicely
        self.w_eff = self.w + self.Delta
        self.kineticInductanceCorrectionOn = True

        if material.lower() == "al":
            self.Tc = 1.23
            self.rho = 4e-9
        elif material.lower() == 'nb':
            self.Tc = 8
            self.rho = 4e-9

        self.l0 = 1.05e-3*np.sqrt(self.rho/self.Tc)

#     Turn on/off corrections due to metal thickness, kinetic inductance
    def setMetalThicknessCorrection(self, on):
        if on:
            self.s_eff = self.s - self.Delta
            self.w_eff = self.w + self.Delta
        else:
            self.s_eff = self.s
            self.w_eff = self.w

    def setKineticInductanceCorrection(self, on):
        self.kineticInductanceCorrectionOn = on

#     Effective Dielectric Constant from Silicon-Air Interface

    def _k0(self):
        return self.w_eff/(self.w_eff+2*self.s_eff)

    def _kp0(self):
        return np.sqrt(1-self._k0()**2)

    def _k1(self):
        return np.sinh(pi*self.w_eff/(4*self.h))/np.sinh(pi*(2*self.s_eff+self.w_eff)/(4*self.h))

    def _kp1(self):
        return np.sqrt(1-self._k1()**2)

    def Eeff(self):
        return 1 + ((self.e1-1)*ellipk(self._k1()**2)*ellipk(self._kp0()**2))/(2*ellipk(self._kp1()**2)*ellipk(self._k0()**2))

#     Kinetic Inductance Calculation
#     from Clem, "Inductances and attenuation constant for a thin-film SC CPW res"

    def _g(self, k, eps):
        return 1/( 2*(1-k)*ellipk(k**2)**2 ) * np.log( 2*(1-k)/eps/(1+k) )

    def _q(self, u):
        return (np.sinh(u) + u)/8/np.sinh(u/2)**2

    def Llk(self):
        """Kinetic inductance per unit length."""
        return mu_0*self.l0/self.w*self._q(self.t/self.l0)*self._g( (self.w/2)/(self.w/2+self.s), 1 )

#     Circuit Parameters

    def Llg(self):
        """Geometric inductance per unit length."""
        return mu_0*ellipk(self._kp0()**2)/(4*ellipk(self._k0()**2))

    def Ll(self):
        """Total Inductance per unit length.  Kinetic inductance corrections can
           be turned on or off."""
        if self.kineticInductanceCorrectionOn:
            return self.Llg() + self.Llk()
        else:
            return self.Llg()

    def Cl(self):
        """Capacitance per unit length"""
        return 4*epsilon_0*self.Eeff()*ellipk(self._k0()**2)/ellipk(self._kp0()**2)

    def Gl(self, w):
        return w * self.tand * self.Cl() * self.e1/(self.e1 + self.e0)

    def vph(self):
        """Phase velocity."""
        return 1/np.sqrt(self.Ll()*self.Cl())

    def z0(self):
        """Characteristic Impedance."""
        return np.sqrt(self.Ll()/self.Cl())

#     Loss
#     from Pozar p.50

    def alpha(self, w):
        return np.sqrt( 1j*w*self.Ll()*( self.Gl(w) + 1j*w*self.Cl() ) ).real

    def beta(self, w):
        return np.sqrt( 1j*w*self.Ll()*( self.Gl() + 1j*w*self.Cl() ) ).imag

#      String

    def __str__(self):
        return "Cl = {} pF/um\nLl = {} nH/um\nvph = {}e8 m/s\nz0 = {} Ohms".format( self.Cl()*1e12/1e6, self.Ll()*1e9/1e6, self.vph()/1e8, self.z0())


class CPWWithBridges(CPW):
    """Currently this does not support kinetic inductance."""
    def __init__(self, bridgeSpacing = 100, bridgeWidth = 2, oxideLossTan = 3e-3, t_oxide = 0.1, e_oxide = 3.9, **kwargs):
        self.bridgeSpacing = bridgeSpacing
        self.bridgeWidth = bridgeWidth
        self.oxideLossTan = oxideLossTan
        self.t_oxide = t_oxide
        self.e_oxide = e_oxide
        super(CPWWithBridges,self).__init__(**kwargs)

    def Cl_dielectric(self):
        """Capacitance per unit length through the substrate."""
        return self.Cl*self.e1/(1+self.e1)

    def Cl_air(self):
        """Capacitance per unit length through the air above the substrate."""
        return self.Cl*1/(1+self.e1)

    def Cl(self):
        """Note that this does not take into account the magnitude of the E-field,
        which changes along a standing wave."""
        cpwCap = super(CPWWithBridges,self).Cl()
        bridgeCap = cap.OverlapCapacitor(self.bridgeWidth*self.w, t = self.t_oxide, eps = self.e_oxide)
        return cpwCap + bridgeCap.cap()/self.l


class HalfLResonator:
    def __init__(self, cpw, l, cin, cout):
        self.l = l*1e-6
        self.cpw = cpw
        self.cki = cin
        self.cko = cout

    def f0(self):
        # return c/(np.sqrt(self.cpw.Eeff())*2*self.l) # Goppl, doesn't include kinetic inductance?
        return self.cpw.vph()/2/self.l

    def w0(self):
        return 2*pi*self.f0() # Goppl

    #     Circuit Parameters with Loss

    def L(self):
        # return 2*self.cpw.Ll()*self.l/(pi**2) # Goppl
        return 1/self.C()/self.w0()**2 # Pozar p.283

    def C(self):
        # these are equiv. (I checked)
        # return self.cpw.Cl()*self.l/2 # Goppl
        return np.pi/2/self.w0()/self.cpw.z0() # Pozar p.283

    def R(self):
        return self.cpw.z0()/(self.cpw.alpha(self.w0())*self.l) # Pozar p.283, Goppl

    def Qint(self):
        return self.R()*np.sqrt(self.C()/self.L()) # Goppl, Pozar p.283, all are equiv.
        # return self.w0()*self.R()*self.C()
        # return np.pi/2/cpw.alpha(self.w0())/self.l
        # return cpw.beta(self.w0())/2/cpw.alpha(self.w0())

    def wn(self):
        return 1/np.sqrt(self.L()*self.C()) # Goppl
        # return self.Qint()/(self.R()*self.C())

    #   Loading

    def Rin(self):
        '''Effective input resistance to ground'''
        return (1. + (self.wn()*self.cki*50.)**2)/(self.wn()*self.cki)**2/50. # Goppl
        # return 50*(1+self.Qs(50, self.cki)**2) # Ted's notes

    def Rout(self):
        '''Effective output resistance to ground'''
        return (1. + (self.wn()*self.cko*50.)**2)/(self.wn()*self.cko)**2/50. # Goppl
        # return 50*(1+self.Qs(50, self.cko)**2) # Ted's notes

    def Cin(self):
        '''Effective input capacitance to ground'''
        return self.cki/(1. + (self.wn()*self.cki*50.)**2)

    def Cout(self):
        '''Effective output capacitance to ground'''
        return self.cko/(1. + (self.wn()*self.cko*50.)**2)

    def wl(self):
        '''Loaded frequency in rad/s'''
        return 1./np.sqrt(self.L()*(self.C() + self.Cin() + self.Cout()))

    def fl(self):
        '''Loaded frequency'''
        return self.wl()/2/pi

    def Qc(self):
        '''Total coupling Q'''
        return 1/(1/self.Qin() + 1/self.Qout())

    def Qin(self):
        '''Output Q'''
        # return self.wn() * (self.C() + self.Cout())/(1./self.R() + 1./self.Rout())
        return self.Rin() * np.sqrt( (self.C() + self.Cout()) / self.L() )

    def Qout(self):
        '''Output Q'''
        # return self.wn() * (self.C() + self.Cout())/(1./self.R() + 1./self.Rout())
        return self.Rout() * np.sqrt( (self.C() + self.Cout()) / self.L() )

    def Ql(self):
        '''Loaded Q, including internal loss'''
        return self.wn() * (self.C() + self.Cin() + self.Cout())/(1./self.R() + 1./self.Rin() + 1./self.Rout())
        # return 1/(1/self.Qint() + 1/self.Qin() + 1/self.Qout())

    def kappa(self):
        '''Photon loss rate'''
        return self.wl()/self.Ql()

    def kin(self):
        return self.wl()/self.Qin()

    def kout(self):
        return self.wl()/self.Qout()

    def photons(self, Pin):
        '''Number of photons in the steady state cavity, based on input power [dBm]'''
        # From Matt Beck's notes on the wiki
        return dBm_to_watts(Pin) * self.Ql()**2 / self.Qc() / hbar / self.wl()**2

    def __str__(self):
        #return "l = {} um\nf = {} GHz\nQ = {}\nk = {} MHz".format(self.cpw.l*1e6,self.fl(), self.Ql(), self.kappa()/2e6/pi)
        return "l = {} um\nf = {} GHz\nQ = {}\nk = {} MHz\nC = {} pF\nL = {} nH\nR = {} Ohms".format(self.l*1e6, self.fl()/1e9, self.Ql(), self.kappa()/2e6/pi, self.C()*1e12, self.L()*1e9, self.R() )
