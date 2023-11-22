#
# units.py
#
# Different units for phonon calculations
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import math

pi = math.pi
clight = 2.99792458e8      # m/s
kb = 1.38064852e-23        # m2*kg/(s2K)
plank = 6.62607004e-34     # m2*kg/s
hbar = plank/2./pi         # m2*kg/s
eps0 = 8.85418782e-12      # m-3kg-1s4A2, vacuum permittivity
AMU = 1.660539040e-27      # kg
Nav = 6.02214086e23        # mol^-1

##
DToCM = (1e-21/clight) # C*m/Debye

## Distance [M]
AToBohr = 1.889725989        # Bohr/A
BohrToA = 0.529177249        # A/Bohr
BohrToM = BohrToA * 1e-10    # m/Bohr

## Energy
EvToJ = 1.60217662e-19       # J/eV
EvToJmol = EvToJ * Nav       # J/mol
EvToRy = (1./13.605698066)

RyToEv = 13.605698066                # eV/Ry
RyToJ = RyToEv * EvToJ               # J/Ry
RyToHz = RyToEv*(EvToJ/(2*pi*hbar))  # Hz/Ry

HaToEv = 2.*RyToEv         # Ev/Hatree

JToEv = 1./EvToJ           # eV/J

HzToJ = 2.*pi*hbar                   # J/Hz
HzToEv = 2.*pi*hbar / EvToJ          # eV/Hz
HzToCm = (1./100./clight)            # cm^-1/Hz

KToJ  = kb                        # J/K
KToEv = kb * JToEv                # eV/K
KToHz = kb / (2.*pi*hbar)         # Hz/K

CmToHz  = clight * 100.              # Hz/cm^-1
CmToTHz = clight * 1e-10             # THz/cm^-1
CmToJ   = CmToHz * HzToJ             # J/cm^-1
CmToEv  = CmToJ  * JToEv             # eV/cm^-1
CmToK   = CmToJ / KToJ               # K/cm^-1
CmToRy  = (CmToEv / RyToEv)          # Ry/cm^-1

RyToCm  = 1./CmToRy                  # cm^-1/Ry
RyToTHz = RyToCm*CmToTHz             # cm^-1/Ry

THzToCm = (1./CmToTHz)

JToCm = 1./CmToJ           # cm^-1/J

AlmCmToHz = clight*100.*2.*pi        # (rad/s)/cm^-1

