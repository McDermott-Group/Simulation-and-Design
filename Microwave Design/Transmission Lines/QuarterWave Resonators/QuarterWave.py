############################################################################
### Function that returns the quarter-wavelength given set of parameters ###
############################################################################
def calcQuarterWavelength(CPW_w, CPW_s, LEA_f0, SUB_er, Cp, Lp):
    # This function calculates the length of a resonator given a CPW geometry
    # and amount of parasitic capacitance (any capacitance to ground outside
    # of the self-capacitance of the CPW itself) and inductance.
    
    # CPW_w = width of the center trace
    # CPW_s = width of gap from center trace to ground plane
    # LEA_f0 = target frequency of resonator after accounting for parasitic capacitances
    # SUB_er = relative permittivity of substrate
    # STR_er = relative permittivity of oxide used in grounding straps over the CPW
    # Cp = extra capacitances beside those from purely the geometry of the CPW
    # Lp = extra inductances beside those from purely the geometry of the CPW
    
    ####################################
    ### Calculate Z0 of the geometry ###
    ####################################
    
    # Elliptic integral nonsense
    CPW_k0 = CPW_w / float(CPW_w + 2*CPW_s)
    CPW_kp = np.sqrt(1-CPW_k0**2)

    # Note that the argument is squared. This is due to the difference in how Gupta
    # defines K(k) and how the function is defined in SciPy.
    CPW_K0 = special.ellipk(CPW_k0**2)
    CPW_Kp = special.ellipk(CPW_kp**2)

    # Capacitance density of CPW to the substrate only.
    # Note: this isn't required to perform the calculation at hand, but is useful for troubleshooting.
    CPW_Csub = 2*epsilon_0*(CPW_K0 / CPW_Kp)

    # Total capacitance per unit length of the CPW geometry (substrate and air included)
    CPW_Cl = 2*epsilon_0*(SUB_er + 1) * (CPW_K0 / CPW_Kp)

    # Total capacitance per unit length replacing substrate with air
    CPW_Cl_air = 4*epsilon_0 * (CPW_K0 / CPW_Kp)

    # Effective relative permittivity of CPW geometry
    CPW_er = CPW_Cl / CPW_Cl_air

    # Phase velocity of this geometry
    CPW_vph = speed_of_light / np.sqrt(CPW_er)

    # CPW Impedance
    CPW_z0 = 1 / (CPW_vph * CPW_Cl)

    ################################################
    ### Calculate f0 of resonator w/o Parasitics ###
    ################################################
    
    LEA_om = 2 * pi * LEA_f0
    
    # From the lumped element model outlined above
    #RES_om = CPW_z0 * Cparasitic * LEA_om**2 / pi + LEA_om*np.sqrt((CPW_z0*Cparasitic*LEA_om/pi)**2+1)
    invDelOm = Lp * pi / float(4 * CPW_z0) + 4 * Cp * CPW_z0 / (pi)
    RES_om = (1/float(2)*(-invDelOm + np.sqrt(-4*Cp*Lp + invDelOm**2 + 4/(LEA_om**2))))**(-1)
    
    # Calculate CPW equivalent lumped elements
    LEA_C = pi / (4 * CPW_z0 * RES_om)
    LEA_L = 1 / (RES_om**2 * LEA_C)
    
    # Go from angular frequency to wavelength
    RES_lambda = 2 * pi * CPW_vph / RES_om

    # Output result in meters as well as other parameters
    output = namedtuple('CPWOutput','quarterlength z0 vph f0 fp Cp Lp er C L')
    CPW = output(RES_lambda/float(4), CPW_z0, CPW_vph, RES_om/(2*pi), LEA_f0,Cp,Lp, CPW_er, LEA_C, LEA_L)
    return CPW

########################################################################
### Function returns the mutual required for desired coupling Q ###
########################################################################
def calcMutual(Z0, L, om, Qc):
    return np.sqrt(2 * Z0 * L / (om * Qc))

###############################################################
### Function returns kappa given the frequency and coupling ###
###############################################################
def calcKappa(om, Qc):
    return om / Qc