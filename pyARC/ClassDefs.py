import pdb, os, sys
import numpy as np
#import RunATMODef, OptimiseDef, SampleDef, TransmissionDef, EmissionDef, Utils
import RunMLEDef, RunEmceeDef, RunNestedSamplingDef


class ARC():
    """
    """

    def __init__( self ):
        """
        ARC object.
        """

        # List containing names of chemical species that
        # will have variable abundances in the retrieval:
        self.FreeSpecies = [ None ] # todo = think more about whether or not there's
                                    # a better way to do this.
        
        # Dictionaries for holding Transmission data and throughputs:
        self.TransmissionData = {}
        self.TransmissionBandpass = {}
        
        # Dictionaries for holding Eransmission data and throughputs:
        self.EmissionData = {}
        self.EmissionBandpass = {}
        
        # The log likelihood function:
        self.LogLikeFunc = None
        self.Priors = {}

        # Controls for MLE optimisation:
        # todo

        # Controls for emcee optimisation:
        self.InitParSampleFuncs = {}
        
        return None


    def RunMLE( self ):
        """
        """
        RunMLEDef.Main( self )
        return None

    
    def RunEmcee( self, nchains=1, nwalkers=100, nsteps=100, threads=1, ncorr_burn=0 ):
        """
        """
        RunEmceeDef.Main( self, nchains=nchains, nwalkers=nwalkers, \
                          nsteps=nsteps, threads=threads, \
                          ncorr_burn=ncorr_burn )
        return None

    
    def RunNestedSampling( self ):
        """
        """
        RunNestedSamplingDef.Main( self )
        return None

    
