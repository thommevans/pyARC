import pdb, os, sys, shutil
import numpy as np
#import RunATMODef, OptimiseDef, SampleDef, TransmissionDef, EmissionDef, Utils
import RunMLEDef, RunEmceeDef, RunNestedSamplingDef


class ARC():
    """
    """

    def __init__( self ):
        """
        ATMO object.
        """

        # Dictionaries for holding Transmission data and throughputs:
        self.TransmissionData = {}
        self.TransmissionBandpass = {}
        
        # Dictionaries for holding Eransmission data and throughputs:
        self.EransmissionData = {}
        self.EransmissionBandpass = {}
        
        # The log likelihood function:
        self.LogLikeFunc = None

        # Controls for MLE optimisation:
        # todo

        # Controls for emcee optimisation:
        
        return None


    def RunMLE( self ):
        """
        """
        RunMLEDef.Main( self )
        return None

    
    def RunEmcee( self ):
        """
        """
        RunEmceeDef.Main( self )
        return None

    
    def RunNestedSampling( self ):
        """
        """
        RunNestedSamplingDef.Main( self )
        return None

    
