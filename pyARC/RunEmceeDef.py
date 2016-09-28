import numpy as np
import emcee
from ProgressBar import progressbar
import pdb, sys, os, time


def Main( ARC, nchains=1, nwalkers=100, nsteps=100, threads=1, ncorr_burn=0 ):

    # Sample for each chain, i.e. group of walkers:
    keys = ARC.InitParSampleFuncs.keys()
    npar = len( keys )
    walker_chains = []
    print '\nRunning the MCMC sampling:'
    for k in range( nchains ):
    
        # Set initial parameter values for each walker:
        p0 = np.zeros( [ nwalkers, npar ] )
        newpriors = {}
        for i in range( npar ):
            p0[:,i] = ARC.InitParSampleFuncs[keys[i]]( nwalkers )

        # Initialise the emcee sampler:
        z = emcee.EnsembleSampler( nwalkers, npar, ARC.LogLikeFunc, threads=threads, args=(keys,ARC) )

        # Run the emcee sampler with progress bar:
        pbar = progressbar( nsteps )
        nstep_increment = max( [ int( np.round( 0.01*nsteps ) ), 1 ] )
        t1 = time.time()
        for i, result in enumerate( z.sample( p0, iterations=nsteps ) ):
            pbar.animate( i+1 )
        t2 = time.time()
        ttot = float( t2-t1 )
        tperlogp = ttot/float( nwalkers*nsteps )
        print '\nTotal time taken = {0:.3f} seconds (for {1:.0f} walkers and {2:.0f} steps'\
              .format( ttot, nwalkers, nsteps )
        print '\nTime per likelihood evaluation = {0:.3f} seconds\n'.format( tperlogp )

        walker_chain_k = {}
        for i in range( npar ):
            walker_chain_k[keys[i]] = z.chain[:,:,i].T
        walker_chain_k['logp'] = z.lnprobability.T
        walker_chains += [ walker_chain_k ]

    return None


def get_chain_from_walkers( walker_chains, acor_integs, ncorr_burn=3 ):
    nchains = len( walker_chains )
    acor = np.zeros( nchains )
    nsteps = np.zeros( nchains )
    nwalkers = np.zeros( nchains )
    for i in range( nchains ):
        walker_chain = walker_chains[i]
        nsteps[i], nwalkers[i] = np.shape( walker_chain['logp'] )
        keys = walker_chain.keys()
        keys.remove( 'logp' )
        npar = len( keys )
        acor_vals = np.zeros( npar )
        for j in range( npar ):
            acor_vals[j] = acor_integs[i][keys[j]]
        acor[i] = np.max( np.abs( acor_vals ) )
    y = nsteps/acor
    if y.min()<ncorr_burn:
        print '\nChains only run for {0:.2f}x correlation times'.format( y.min() )
        pdb.set_trace()
    else:
        acor = acor.max()
        nburn = int( np.round( ncorr_burn*acor ) )
        chain_dicts = []
        chain_arrs = []
        for i in range( nchains ):
            chain_dicts += [ collapse_walker_chain( walker_chains[i], nburn=nburn ) ]
        grs = gelman_rubin( chain_dicts, nburn=0, thin=1 )
        chain = combine_chains( chain_dicts, nburn=nburn, thin=1 )        
    return chain, grs, nchains, nwalkers, nsteps, acor, nburn


def collapse_walker_chain( walker_chain, nburn=0 ):
    """
    Takes emcee walker chains and collapses them into single
    chains, one for each parameter. An optional burn-in range
    from can be discarded from the beginning of each walker
    chain prior to combining them together.
    """
    if nburn==None:
        nburn = 0
    keys = walker_chain.keys()
    npar = len( keys )
    chain = {}
    for key in keys:
        chain[key] = walker_chain[key][nburn:,:].flatten()
    return chain


def walker_chain_autocorr( walker_chain, nburn=None, minlag=1, maxlag=50 ):
    """
    Computes the autocorrelation function and integrated autocorrelation times
    for each parameter, using the routines from emcee.
    """

    if nburn==None:
        nburn = 0

    keys = walker_chain.keys()
    keys.remove( 'logp' )
    npar = len( keys )
    nsteps, nwalkers = np.shape( walker_chain['logp'] )
    y = np.zeros( [ nsteps-nburn, npar ] ) 
    for i in range( npar ):
        y[:,i] = np.mean( walker_chain[keys[i]][nburn:,:], axis=1 ) # average over walkers
    
    acor_func_arr = emcee.autocorr.function( y, axis=0, fast=False )
    acor_integ_arr = emcee.autocorr.integrated_time( y, axis=0, c=1, fast=False, \
                                                     low=minlag, high=maxlag )

    acor_func = {}
    acor_integ = {}
    for i in range( npar ):
        acor_func[keys[i]] = acor_func_arr[:,i]
        acor_integ[keys[i]] = acor_integ_arr[i]

    return acor_func, acor_integ
    

def combine_chains( chain_list, nburn=0, thin=1 ):
    """
    Combines multiple chains into a single chain.

    CALLING

      new_chain = pyhm.combine_chains( chain_list, nburn=1000, thin=5 )

    The above would take a list of chains, cut 1000 samples from the start
    of each and thin by a factor of 5 before combining into a single chain.
    This is done separately for each parameter in the chain.
    """
    
    m = len( chain_list )
    keys = chain_list[0].keys()
    combined = {}
    for key in keys:
        combined[key] = []
        for i in range( m ):
            chain = chain_list[i][key][nburn:]
            n = len( chain.flatten() )
            ixs = ( np.arange( n )%thin==0 )
            combined[key] += [ chain[ixs] ]
        combined[key] = np.concatenate( combined[key] )

    return combined


def gelman_rubin( chain_list, nburn=0, thin=1 ):
    """
    Calculates the convergence statistic of Gelman & Rubin (1992).

    CALLING
    
        pyhm.gelman_rubin( chain_list, nburn=10000 )
    
      The above takes a list of chains as input and calculates the
      Gelman-Rubin convergence statistic for each parameter, after
      discarding nburn steps from the start of each chain. Note that
      each element of the chain list should be a dictionary, and each
      labelled element of each dictionary should be an array containing
      the samples for each parameter, i.e. the standard format in which
      chains are stored when they are generated as attributes of Sampler
      objects.

    DESCRIPTION

      The Gelman-Rubin convergence statistic is calculated as:

         B = ( n/( m-1 ) )*sum_j[ ( theta_j - theta_mean )^2  ]

         W = ( 1/m )*sum_j[ ( 1/( n-1 )*sum_i[ ( theta_ji - theta_j )^2 ] ) ]

         GelmanRubin = sqrt[ ( n-1 )/n + ( 1/n )*( B/W ) ]

      where m is the number of chains, n is the number of samples in each
      chain, the j summation is done over the m chains, the i summation is
      done over the n samples within the jth chain, theta_j is the mean of
      the jth chain, theta_mean is the mean of all samples across all m chains,
      and theta_ji is the ith sample of the jth chain.

      If the GelmanRubin value is close to 1, it suggests that the chains are
      all well-mixed and have likely reached a stable state. If the value is
      more than a few percent from 1, it suggests that the chains have not yet
      converged and more samples are required.

      Note that the formula given above does not include the scaling factor of
      df/(df-2) that Gelman & Rubin originally included. Brooks & Gelman (1998)
      subsequently pointed out that this factor was incorrect anyway, and that
      the correct factor should be (d+3)/(d+1). However, if the chains are close
      to converged, this scaling factor is very close to 1 anyway, which is why
      it is ignored in the pyhm implementation.
    """
    m = len( chain_list )

    keys = []
    for key in chain_list[0].keys():
        if ( key=='logp' )+( key=='accepted' ):
            continue
        else:
            keys += [ key ]

    n = len( chain_list[0][keys[0]] ) - nburn
    ixs = ( np.arange( n )%thin==0 )
    npars = len( keys )
    grs = {}
    for i in range( npars ):
        chains = []
        for j in range(m):
            chains += [ chain_list[j][keys[i]][nburn:][ixs] ]
        chains = np.column_stack(chains)
        W = np.mean( np.var( chains, axis=0, ddof=1 ) )
        B_over_n = np.var( np.mean( chains, axis=0 ), ddof=1 )
        sigma2 = ( ( n-1. )/n )*W + B_over_n
        Vhat = sigma2 + B_over_n/float( m )
        grs[keys[i]] = np.sqrt( Vhat/float( W ) )
                 
    return grs

