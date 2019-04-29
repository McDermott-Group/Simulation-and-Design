__author__ = 'matthewbeck'
__email__ = 'mabeck2@wisc.edu'

import numpy as np
from numpy import *
import scipy as sp
from scipy import *
from scipy import constants

# Permeability of free space
mu0 = sp.constants.mu_0

# __________          ________           _______          ________          _________
#      g    |<- s1 ->|<--w1-->|<--s1--> |<--d-->|<--s2-->|<--w2-->|<--s2-->|   g

# Enter the dimensions of the coupled CPW's below in units of microns
w1 = 6.
w2 = 6.
s1 = 3.
s2 = 3.
d = 10.

# Ground planes are mult times longer than center trace
mult = 10.
g = mult * w1

# Function to calculate inductive coupling per unit meter
def CalcMutual(w1Width, w2Width, s1Width, s2Width, dWidth, gWidth):
    # Convert coupled CPW dimensions into nanometers
    w1Width = w1Width * 1000.
    w2Width = w2Width * 1000.
    s1Width = s1Width * 1000.
    s2Width = s2Width * 1000.
    dWidth = dWidth * 1000.
    gWidth = gWidth * 1000.

    # Number of points in trace 1
    points = 100

    # current return loop radius (Just needs to be large compared to cpw dimensions)
    returnLoop = 1.0e6

    # Zero temperature Penetration depth of superconductor in nm (100 nm thick Nb ~ 88 nm)
    penDepth = 88.0

    # number of points in W1
    w1Points = points

    # number of points in W2
    w2Points = np.round(w2Width / w1Width * points)

    # number of points in ground plane
    groundNumPoints = mult * points

    # number of points in D
    dNumPoints = np.round(dWidth / w1Width * points)

    # Left side ground plane vector
    ground1 = sp.linspace(
        -(dWidth / 2. + s1Width + w1Width + s1Width + gWidth),
        -(dWidth / 2. + s1Width + w1Width + s1Width),
        groundNumPoints
    )

    # W1 vector
    trace1 = sp.linspace(
        -(dWidth / 2. + s1Width + w1Width),
        -(dWidth / 2. + s1Width),
        w1Points
    )

    # D Vector
    separation = sp.linspace(
        -dWidth / 2.,
        dWidth / 2.,
        dNumPoints
    )

    # W2 Vector
    trace2 = sp.linspace(
        dWidth / 2. + s2Width,
        dWidth / 2. + s2Width + w2Width,
        w2Points
    )

    # right side ground plane vector
    ground2 = sp.linspace(
        dWidth / 2. + s2Width + w2Width + s2Width,
        dWidth / 2. + s2Width + w2Width + s2Width + gWidth,
        groundNumPoints
    )

    # Combine all vectors into master vector
    xVec = sp.concatenate((ground1, trace1, separation, trace2, ground2))

    # Length of y vector (Set to 0 for zero thickness approx)
    yVec = [0]

    # dummy value for KRON
    yDum = 1.

    # matrix consisting of CPW sites
    rMat = sp.kron(xVec, yDum)

    # total number of sites
    arrayLength = len(xVec) * len(yVec)

    # reshape rMat into vector
    rVec = sp.reshape(rMat.T, (arrayLength, 1))

    # dummy mutual matrix for KRON
    mDum = sp.ones((1, len(rVec)))

    # site positions KRON'd
    arg1 = sp.kron(mDum, rVec)

    # transpose of site positions KRON'd
    arg2 = sp.kron(rVec.T, mDum.T)

    # r_ij is the absolute distance between sites of the coupled CPWs. Obviouslly then the diagonal is 0's
    r_ij = abs(arg1 - arg2)

    # x-step in gournd plane 1
    dxG = gWidth / groundNumPoints

    # x-step in D
    dxD = dWidth / dNumPoints

    # x-step in trace 1
    dxW1 = w1Width / w1Points

    # x-step in trace 2
    dxW2 = w2Width / w2Points

    # self inductance of ground plane 1
    selfGround1 = dxG / 2. * sp.ones(len(ground1))

    # self inductance of trace 1
    selfW1 = dxW1 / 2. * sp.ones(len(trace1))

    # self inductance of D
    selfD = dxD / 2. * sp.ones(len(separation))

    # self inductance of trace 2
    selfW2 = dxW2 / 2. * sp.ones(len(trace2))

    # self inductance of ground plane 2
    selfGround2 = dxG / 2. * sp.ones(len(ground2))

    # combine self inductances into single array
    self = sp.concatenate((selfGround1, selfW1, selfD, selfW2, selfGround2))

    # put values of self of diagonal of 2D array of size r_ij
    selfDiag = sp.diag(self)

    r_ij = r_ij + selfDiag

    # Add self inductance to diagonal of r_ij
    r_ij = r_ij + selfDiag

    # Divide return current loop by position matrix (exponentiated M matrix)
    m_exp = returnLoop / r_ij

    # geometric contribution
    mGeo = mu0 / 2. / sp.pi * sp.log(m_exp)

    # kinetic inductance contribution
    lKinVec = mu0 * penDepth ** 2. / dxW1 ** 2. * np.ones(len(rVec))

    # put kinetic inductance values on diagonal of matrix size mGeo
    lKin = np.diag(lKinVec)

    # total inductance matrix
    m = mGeo + lKin

    # set phi of ground 1
    phiG1 = -1. * np.ones(len(ground1))

    # set phi of trace 1
    phiW1 = 1. * np.ones(len(trace1))

    # set phi of D
    phiD = -1. * np.ones(len(separation))

    # set phi fo ground 2
    phiG2 = -1. * np.ones(len(ground2))

    # initialize iSum
    iSum = np.zeros(len(ground1))
    iSum = iSum.reshape(len(iSum), 1)

    # initialize iCurrent
    iCurrent = sp.ones(len(ground1))
    iCurrent = iCurrent.reshape(len(iCurrent), 1)

    # Initial value for phi 2. range is from [-1,1]
    phiW2Val = -1.

    # initial j which is a place holder / counter for the while loop
    j = 0

    # initial step in phi
    dPhi = .05

    while phiW2Val < 1:
        phiW2 = phiW2Val * sp.ones(len(trace2))
        phi = sp.concatenate((phiG1, phiW1, phiD, phiW2, phiG2))
        phi = phi.reshape(len(phi), 1)
        current = np.linalg.solve(m, phi)
        current = current.reshape(len(current), 1)
        iSum[j] = sp.sum(
            current[int(groundNumPoints + points + dNumPoints):int(groundNumPoints + points + dNumPoints + points)])

        if j > 0:
            curValue = iSum[j]
            pastValue = iSum[j - 1]
            sign = curValue * pastValue
            if sign < 0:
                phiW2Val = phiW2Val - dPhi
                dPhi = dPhi / 10.
                j = j - 1
            if np.abs(curValue) < 1.0e-6:
                iCurrent[j] = phiW2Val
                iCurrent[j + 1] = phiW2Val
                break
        iCurrent[j] = phiW2Val
        j = j + 1
        phiW2Val = phiW2Val + dPhi

    phiW2I = iCurrent[j]
    i1Sum = sp.sum(current[int(groundNumPoints):int(groundNumPoints + points)])
    indLength = (phiW2I + 1.0) / i1Sum
    print('The inductance per unit length is %.2e H/m' % indLength)
    return indLength


CalcMutual(w1,w2,s1,s2,d,g)


# CalcMutual can easily be called from a for/while loop to explore many different geometries.