import pdb, sys, os, time
import numpy as np
import scipy.interpolate

"""
This module contains routines for evaluating the log likelihood for some standard atmosphere models.
"""


def ClearChemEqTransmission( ATMO, pars, priors, DataArr, UncArr ):
    """
    Evaluates the log likelihood for a clear atmosphere in chemical equilibrium
    assuming an isothermal PT profile given an observed spectrum for the planet. 
    Free parameters are a vertical offset (dRpRs), planetary temperature (Teff), 
    atmospheric metallicity (MdH), and the C/O ratio (COratio).

    Inputs:
    ATMO - An ATMO object with namelist parameters set appropriately.
    pars - Dictionary containing values for each of the free parameters.
    priors - Dictionary containing prior functions for each of the free parameters.
    DataArr - Dictionary with separate named entries for each dataset comprised
              of Nix2 array where Ni is the number of datapoints in the ith dataset.
    UncArr - Same as DataArr except that it contains the associated measurement
             uncertainties.
    """
    
    t1 = time.time()

    # Evaluate the prior likelihood:
    logp_prior = 0
    for key in pars.keys:
        logp_prior = priors[key]( pars[key] )

    # Install values for the free parameters:
    ATMO.teff = pars['Teff']
    ATMO.MdH = pars['MdH']
    ATMO.COratio = pars['COratio']

    # Use tempfile to create input and output files so that there
    # will be no duplication e.g. if running many walkers: 
    tempfileobj1 = tempfile.NamedTemporaryFile( mode='w+b', delete=False )
    ATMO.ftrans_spec = tempfileobj1.name
    tempfileobj2 = tempfile.NamedTemporaryFile( mode='w+b', delete=False )
    ATMO.infile_path = tempfileobj2.name

    # Compute the model transmission spectrum:
    ATMO.RunATMO()
    ATMO.ReadTransmissionModel( ncdf_fpath=tempfileobj1.name )
    WavMicronModel = ATMO.TransmissionModel[:,0]
    RpRsModel = ATMO.TransmissionModel[:,1] - pars['dRpRs']

    # Bin the transmission spectrum into the data bandpasses:
    interpf = scipy.interpolate.interp1d( WavMicronModel, RpRsModel )
    datasets = ATMO.TransmissionData.keys()
    ndatasets = len( datasets )
    ModelArr = []
    for i in range( ndatasets ):
        dataseti = ATMO.TransmissionData[datasets[i]]
        ledges = dataseti[:,0]
        uedges = dataseti[:,1]
        nchannels = len( ledges )
        ModelArri = np.zeros( nchannels )
        for j in range( nchannels ):
            # Binning the model could be done more carefully with
            # actual instrument throughputs defined:
            WavChannel = np.r_[ ledges[j]:uedges[j]:1j*100 ]
            ModelArri[j] = np.mean( interpf( WavChannel ) )
        ModelArr += [ ModelArri ]
    ModelArr = np.concatenate( ModelArr )

    # Compute the residuals and data log likelihood:
    ResidsArr = DataArr - ModelArr
    ndat = len( ResidsArr )
    logp_data = logp_mvnormal_whitenoise( ResidsArr, UncArr, ndat )
    os.remove( tempfileobj1.name )
    os.remove( tempfileobj2.name )
    # TODO = COULD SAVE ALL OF THESE TRANSMISSION SPECTRA
    # FOR PLOTTING AT THE END

    logp = logp_prior + logp_data
    t2 = time.time()
    return logp


def logp_mvnormal_whitenoise( r, u, n  ):
    """
    Log likelihood of a multivariate normal distribution
    with diagonal covariance matrix.
    """
    term1 = -np.sum( np.log( u ) )
    term2 = -0.5*np.sum( ( r/u )**2. )
    return term1 + term2 - 0.5*n*np.log( 2*np.pi )

    
