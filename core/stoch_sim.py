import numpy as np
import scipy.stats


def gaussian_gen(mux, muy, cov):
    """Create 2D Gaussian Distribution.
    
    Parameters
    ----------
    mux : float
        Average of Gaussian distribution for variable x.

    muy : float
        Average of Gaussian distribution for variable y.

    cov : list, array_like
        Covariance list.

    Returns
    -------
    rv : scipy.stats._multivariate.multivariate_normal_gen object
        Gaussian distribution object.

    rv_pdf : array_like
        PDF for bivariate Gaussian distribution.
    """

    gridx = np.linspace(0, 1, 100)
    gridy = np.copy(gridx)

    X, Y = np.meshgrid(gridx, gridy)

    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X; pos[:, :, 1] = Y
    
    rv = scipy.stats.multivariate_normal([mux, muy], cov)
    rv_pdf = rv.pdf(pos)

    return rv, rv_pdf


def sampling_from_gaussian(rv, depth, facies):
    """Random samples from bivariate Gaussian distribution

    Parameters
    ----------
    rv : scipy.stats._multivariate.multivariate_normal_gen object
        Gaussian distribution object
    
    depth : array_like
        Depth curve.

    facies : array_like
        Facies boolean.

    Returns
    -------
    avr_x : float
        Average of variable x.

    avr_y : float
        Average of variable y.

    cov : array_like
        Covariance matrix between variable x and y.
    """

    samples = rv.rvs(size=len(depth[facies]))

    avr_x = np.mean(samples[:, 0])
    avr_y = np.mean(samples[:, 1])
    cov = np.cov([samples[:, 0], samples[:, 1]])

    return avr_x, avr_y, cov


def syn_variogram(depth, C, a):
    """Generate a exponential variogram model.

    Parameters
    ----------
    depth : array_like
        Depth curve.

    C : float
        Sill of variogram.

    a : float
        Range of variogram.

    Returns

    model : array_like
        Covariance model generated using exponential variogram model.
    """

    h = np.linspace(0, len(depth), len(depth))

    exp = C*(1 - np.exp((-3 * h) / a))
    model = (1 - exp/C)

    return model


def simulations(model, cov, avr_x, avr_y):
    """Simulations workflow.
    
    Parameters
    ----------
    model : array_like
        Covariance model.

    cov : array_like
        Covariance matrix between variables x and y.

    avr_x : array_like
        Trend for variable x.

    avr_y : array_like
        Trend for variable y.

    Returns
    -------
    sim_x : array_like
        Simulation for variable x.

    sim_y : array_like
        Simulation for variable y.
    """

    c = np.zeros((len(model), len(model)))

    for i in range(len(model)):
        for j in range(len(model)):
            
            if i == j:
                
                c[i][j] = model[0]
                
            elif (j > i):
                
                c[i][j] = model[j - i]
                    
            elif (j < i):
                    
                c[i][j] = model[i - j]

    K = np.kron(cov, c)

    R = np.linalg.cholesky(K)

    u = np.random.normal(loc=0, scale=1, size=2*len(model))

    w = np.dot(R, u)

    sim_x = avr_x + w[0:len(model)]
    sim_y = avr_y + w[len(model):]

    sim_x = np.clip(sim_x, 10**-5, 1.0)
    sim_y = np.clip(sim_y, 10**-5, 1.0)

    return sim_x, sim_y
