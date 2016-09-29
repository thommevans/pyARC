#import pyARC
from pyARC_dev import pyARC
import numpy as np
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import glob
import pdb, sys, os


"""
Module containing a routine for installing the namelist
parameters in an ATMO object.
"""


# Fiducial values for HAT-P-18:
LOGG = 2.69
TEQ = 850
RPLANET = 0.995
RSTAR = 0.749
AAU = 0.0559
MDH = 0.0 # solar
CORATIO = 0.56 # solar



def Main( ATMO ):
    """
    Routine for setting the namelist parameters.
    """
    
    ATMO.executable = 'atmo.x'
    ATMO.nice = None
    ATMO.infile_path = 'ptchem.in' # ATMO input file path

    # PARAM: Parameters for ATMO
    ATMO.Debug = 1
    ATMO.fout = 'pt.ncdf'
    #ATMO.fin = 'temp.ncdf'
    ATMO.fin = 'pt.ncdf'

    #  EOS: Equation of state parameters
    #ATMO.gamma = 0.

    # GRID: Grid parameters
    ATMO.pmin = 1e-6
    ATMO.pmax = 5
    ATMO.taumin = 1e-6
    ATMO.taumax = 10 #8e3 #2e2
    ATMO.logg = LOGG
    ATMO.teff = 100.
    ATMO.ndepth = 15 #20 # 50 # this seems to be the parameter that affects speed a lot
    ATMO.Rp = RPLANET
    ATMO.pp_Rp = 0.001
    ATMO.nfreq = 250
    ATMO.nkmix = 10 #30
    ATMO.nband = 250 #350
    ATMO.nband_std = 32
    ATMO.corr_k = True
    ATMO.numax = 5e6

    # CHEMISTRY: Chemistry parameters
    ATMO.chem = 'eq'
    ATMO.MdH = MDH
    ATMO.COratio = CORATIO
    ATMO.fAin = 'chem_dummy.ncdf'
    ATMO.fAeqout = 'chem_eq.ncdf'
    ATMO.fAneqout = 'chem_neq.ncdf'
    #ATMO.fcoeff = '/home/tevans/code/pacode/atmo/chem/coeff_NASA_sc.dat'
    ATMO.print_chem = False

    # CHEM_NEQ: Non-equilibrium chemistry parameters
    ATMO.mixing = False
    ATMO.photochem = False
    ATMO.kzzcst = 1e9
    ATMO.nmol_eq = 107
    ATMO.tmax = 1e12
    ATMO.dtmax = 1e10
    ATMO.rate_limiter = True
    ATMO.Nmin = 1e-100
    ATMO.atol = 1e-10

    # RADTRANS: Radiative transfer parameters
    ATMO.nrays = 16
    ATMO.scatter = True
    ATMO.irrad = True
    ATMO.firad = 'lte048-4.5-0.0a+0.0.BT-NextGen.7.ncdf' # stellar spectrum filepath
    ATMO.rstar = RSTAR
    ATMO.rorbit = AAU
    ATMO.murad = 0.5
    ATMO.fred = 0.5
    ATMO.ftrans_spec = 'TransmissionModel.ncdf' # output file for transmission spectrum
    ATMO.fspectrum = 'EmissionModel.ncdf' # output file for the emission spectrum
    ATMO.fcfout = 'ContribFunc.ncdf' # output file for the normalised contribution function    
    
    # OPACITY: Opacity parameters
    ATMO.nkap = 6
    ATMO.art_haze = 1
    ATMO.cloud = False
    ATMO.cloud_top = 1
    ATMO.cloud_bottom = 20
    ATMO.cloud_strength = 1
    ATMO.kap_smooth = True
    ATMO.kerkap_smooth = 2

    # SOLVER: ATMO solver parameters
    ATMO.solve_hydro = False
    ATMO.solve_energy = False
    ATMO.minstep = 1e-3
    ATMO.maxstep = 9e-1
    ATMO.accuracy = 1e-1
    ATMO.psurf = 1e-6
    ATMO.print_err = False
    ATMO.transmission_spectrum = True
    ATMO.surface_spectrum = False
    #ATMO.hydrostatic = True
    ATMO.calc_cf = False
    
    # CONVECTION: Convection parameters
    ATMO.alpha = 0.

    return ATMO
