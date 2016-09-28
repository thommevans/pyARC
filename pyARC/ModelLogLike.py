import pdb, sys, os, time
import numpy as np
import scipy.interpolate
import tempfile


"""
This module contains routines for evaluating the log likelihood for some standard atmosphere models.
"""


def ClearChemEqTransmission( parsarr, keys, ARC ):
    """
    Evaluates the log likelihood for a clear atmosphere in chemical equilibrium
    assuming an isothermal PT profile given an observed spectrum for the planet. 
    Free parameters are a vertical offset (dRpRs), planetary temperature (Teff), 
    atmospheric metallicity (MdH), and the C/O ratio (COratio).

    Inputs:
    parsarr - Array containing values for each of the model parameters.
    keys - String labels for each parameter in parsarr that can be used to 
           map to the prior functions.
    ARC - An ARC object.
    DataArr - Dictionary with separate named entries for each dataset comprised
              of Nix2 array where Ni is the number of datapoints in the ith dataset.
    UncArr - Same as DataArr except that it contains the associated measurement
             uncertainties.
    """
    
    t1 = time.time()

    datasets = ARC.TransmissionData.keys()
    ndatasets = len( datasets )
    DataArr = []
    UncArr = []    
    for i in range( ndatasets ):
        dataseti = ARC.TransmissionData[datasets[i]]
        DataArr += [ dataseti[:,2] ]
        UncArr += [ dataseti[:,3] ]        
    DataArr = np.concatenate( DataArr )
    UncArr = np.concatenate( UncArr )

    npar = len( keys )
    pars = {}
    for i in range( npar ):
        pars[keys[i]] = parsarr[i]

    
    # Evaluate the prior likelihood:
    logp_prior = 0
    for key in pars.keys():
        logp_prior += ARC.Priors[key]( pars[key] )

    # Install values for the free parameters:
    ARC.ATMO.teff = pars['Teff']
    ARC.ATMO.MdH = pars['MdH']
    ARC.ATMO.COratio = pars['COratio']

    # Use tempfile to create input and output files so that there
    # will be no duplication e.g. if running many walkers: 
    tempfileobj1 = tempfile.NamedTemporaryFile( mode='w+b', delete=False )
    ARC.ATMO.ftrans_spec = tempfileobj1.name
    tempfileobj2 = tempfile.NamedTemporaryFile( mode='w+b', delete=False )
    ARC.ATMO.infile_path = tempfileobj2.name

    # Compute the model transmission spectrum:
    ARC.ATMO.RunATMO()
    ARC.ATMO.ReadTransmissionModel( ncdf_fpath=tempfileobj1.name )
    WavMicronModel = ARC.ATMO.TransmissionModel[:,0]
    RpRsModel = ARC.ATMO.TransmissionModel[:,1] - pars['dRpRs']

    # Bin the transmission spectrum into the data bandpasses:
    interpf = scipy.interpolate.interp1d( WavMicronModel, RpRsModel )
    datasets = ARC.TransmissionData.keys()
    ndatasets = len( datasets )
    ModelArr = []
    for i in range( ndatasets ):
        dataseti = ARC.TransmissionData[datasets[i]]
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

    
