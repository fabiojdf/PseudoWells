import numpy as np
import random


def mcmc(nfacies, depth, P, code, seed, sf=1):
    """Markov Chain Monte Carlo Simulation

    Parameters
    ----------

    nfacies : int
        Number of facies.

    depth : array_like
        Depth curve.

    P : array_like
        Transition matrix.

    code : list
        List of names of each facies from 1 to nfacies.

    sf : int
        Start facies (facies on top of interval).

    Returns
    -------

    facies : array_like
        Random genetared facies log (number).

    code_facies : array_like
        Facies log with code (ex: (['shale', 'sand'...)]).
    """

    facies = np.empty(len(depth))
    facies[0] = sf

    a = np.arange(0, nfacies)

    zipf = zip(a, code)
    dic = dict(zipf)
    code_facies = []
    code_facies.append(dic.get(sf))
    
    if seed==True:
        random.seed(14)

    for i in range(1, len(depth)):

        row = facies[i-1]
        dist=[]
        for j in range(nfacies):
            dist.append(P[j, int(row)])
        facies[i] = random.choices(a, dist)[0]
        code_facies.append(dic.get(facies[i]))


    return facies[::-1], code_facies[::-1]
    
    





