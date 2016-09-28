from pyATMO_dev import pyATMO
import InstallNamelistParameters
import numpy as np
import matplotlib.pyplot as plt
import pdb, sys, os, time


"""
Module for generating a synthetic transmission spectrum. Execution is:
>> import PrepDataSet
>> PrepDataset.MakeTransmissionData( make_plot=False )
"""


LOGG = InstallNamelistParameters.LOGG
TEQ = InstallNamelistParameters.TEQ
RPLANET = InstallNamelistParameters.RPLANET
RSTAR = InstallNamelistParameters.RSTAR
AAU = InstallNamelistParameters.AAU
MDH = InstallNamelistParameters.MDH
CORATIO = InstallNamelistParameters.CORATIO


TRANSMISSIONMODEL = 'TransmissionModelForFitting.ncdf'
TRANSMISSIONDATA = 'TransmissionDataForFitting.txt'


def MakeTransmissionData( make_plot=False ):

    MakeTransmissionModel()

    # Uncertainty on RpRs data:
    RpRsSig = 4e-4
    
    # Edges for the wavelength channels:
    chedges = np.linspace( 1.12, 1.6, 28 )

    # Read in the previously-generated transmission model:
    ATMO = pyATMO.ATMO()
    ATMO.ReadTransmissionModel( ncdf_fpath=TRANSMISSIONMODEL )
    wav_micron = ATMO.TransmissionModel[:,0]
    RpRsModel = ATMO.TransmissionModel[:,1]

    # Bin the model into the wavelength channels:
    nchannels = len( chedges )-1
    RpRsData = np.zeros( nchannels )
    for i in range( nchannels ):
        ixs = ( wav_micron>=chedges[i] )*( wav_micron<chedges[i+1] )
        RpRsData[i] = np.mean( RpRsModel[ixs] )

    # Add measurement uncertainties:
    RpRsData += RpRsSig*np.random.randn( nchannels )
    RpRsUncs = RpRsSig*np.ones( nchannels )

    # Save to file:
    output = np.column_stack( [ chedges[:-1], chedges[1:], RpRsData, RpRsUncs ] )
    np.savetxt( TRANSMISSIONDATA, output )

    if make_plot==True:
        x = 0.5*( chedges[:-1] + chedges[1:] )
        plt.ion()
        plt.figure()
        plt.plot( wav_micron, RpRsModel, '-c' )
        plt.plot( x, RpRsData, 'o', ms=10, mfc='r', mec='r' )
        plt.errorbar( x, RpRsData, yerr=4e-4, fmt='ok' )
        plt.xlim( [ 1, 2 ] )
        plt.xlabel( 'Wav (micron)' )
        plt.ylabel( 'Rp/Rs' )
    
    return None


def MakeTransmissionModel():
    ATMO = pyATMO.ATMO()
    InstallNamelistParameters.Main( ATMO )
    ATMO.Debug = 0
    ATMO.fin = None
    ATMO.nband = 5000
    ATMO.teff = TEQ
    ATMO.logg = LOGG
    ATMO.MdH = MDH
    ATMO.COratio = CORATIO
    ATMO.art_haze = 1.
    ATMO.cloud = False
    ATMO.solve_hydro = False
    ATMO.solve_energy = False
    ATMO.transmission_spectrum = True
    ATMO.surface_spectrum = False
    ATMO.ftrans_spec = TRANSMISSIONMODEL
    ATMO.RunATMO()
    return None


