"""
Chris Wilen
McDermott Group, UW Madison

This file provides objects for a variety of different types of capacitors,
allowing one to easily calculate the capacitance of various geometries.
"""

from scipy.constants import c, epsilon_0, mu_0, h, hbar, e, pi
from scipy.special import ellipk, ellipkm1
import numpy as np

class OverlapCapacitor(object):

    def __init__(self, area, t=0.1, eps_r=3.9):
        """
        area = area of overlap in um^2.
        t = thickness of dielectric in um.
        """
        self.area = area*1e-12
        self.t = t*1e-6
        self.eps_r = eps_r

    def cap(self):
        return self.eps_r*epsilon_0*self.area/self.t


class MicrostripGap(object):
    """ Source: http://qucs.sourceforge.net/tech/node79.html """

    def __init__(self, w1, w2, s, eps_r=3.9, h=.1):
        """
        w1 = width of first microstrip in um
        w2 = width of second microstrip in um
        s  = distance between the two ends of the microstrip in um
        eps_r = relative permittivity of dielectric
        h  = ?? height of dielectric (I think) in um
        """
        self.w1 = w1*1e-6
        self.w2 = w2*1e-6
        self.s  = s *1e-6
        self.h  = h *1e-6
        self.eps_r = eps_r

    def Q1(self):
        return 0.04598*(0.03 + (self.w1/self.h)**self.Q5())*(0.272+0.07*self.eps_r*epsilon_0)

    def Q2(self):
        return 0.107*(self.w1/self.h + 9)*(self.s/self.h)**1.05 * (1.5+0.3*self.w1/self.h)/(1+0.6*self.w1/self.h)

    def Q3(self):
        return np.exp(-0.5978*(self.w2/self.w1)**1.35)-0.55

    def Q4(self):
        return np.exp(-0.5978*(self.w1/self.w2)**1.35)-0.55

    def Q5(self):
        return 1.23/( 1+0.12*(self.w2/self.w1-1)**0.9 )

    def cap(self):
        """Capacitance between the ends of the microstrips in Farads"""
        return 500*self.h*np.exp(-1.86*self.s/self.h)*self.Q1()*(1+4.19*(1-np.exp(-0.785*np.sqrt(self.h/self.w1)*self.w2/self.w1)))

    def capToGround1(self, openEndCap):
        """Capacitance to ground from the end of each microstrip
        in Farads.  The open end capacitance can be found from
        http://qucs.sourceforge.net/tech/node78.html#eq:Cend.  We
        probably want to implement this here, but I don't think we
        ever use this part, so until someone does, it will probably
        not be implemented."""
        return openEndCap*(self.Q2()+self.Q3())/(self.Q2()+1)

    def capToGround2(self, openEndCap):
        """Capacitance to ground from the end of each microstrip
        in Farads.  The open end capacitance can be found from
        http://qucs.sourceforge.net/tech/node78.html#eq:Cend.  We
        probably want to implement this here, but I don't think we
        ever use this part, so until someone does, it will probably
        not be implemented."""
        return openEndCap*(self.Q2()+self.Q4())/(self.Q2()+1)

    def Y(self,w,openEndCap1=0,openEndCap2=0):
        """Admittance parameters at frequency w=2*pi*f.
        Open ended capacitances default to 0, which eliminates
        the correction from the matrix below."""
        cp1 = self.capToGround1(openEndCap1)
        cp2 = self.capToGround2(openEndCap2)
        return 1j*w*np.matrix( [[cp1+self.cap(), -self.cap()   ],
                                [-self.cap()   , cp1+self.cap()]] )


class Interdigitated(object):
    """Source: Microstrip Lines and Slotlines by Gupta, p.133.
    Note that Gupta aproximated eps_0 by 10^-3/36pi.  I made it exact."""

    def __init__(self, w, s, l, n, eps_r=3.9):
        """
        w = width of finger in um
        s = gap between fingers in um (same on end and side)
        l = length of each finger in um (only the overlapping part)
        n = number of fingers.  must be same on both sides
        eps_r = relative permittivity of dielectric
        """
        self.w = w*1e-6
        self.s = s*1e-6
        self.l = l*1e-6
        self.n = n
        self.eps_r = eps_r

    def _a(self):
        return self.w/2

    def _b(self):
        return (self.w + self.s)/2

    def _k(self):
        return np.tan(self._a()*np.pi/4/self._b())**2

    def cap(self):
        return 2*self.eps_r*epsilon_0*self.l*(self.n - 1)*ellipk(self._k()**2)/ellipk(1-self._k()**2)
