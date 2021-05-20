import numpy as np


def soft_sand(K, G, rho, phi, phic, n, P):
    """
    Function for calculation of soft sand rock physics model.

    Parameters
    ----------
    K : float
        Mineral bulk modulus.

    G : float
        Mineral shear modulus.

    rho : float
        Mineral density.

    phi : array_like
        Porosity array.
    
    phic : float
        Critical porosity.

    n : float
        Coordination number.

    P : float
        Hidrostatic Stress.

    Returns
    -------
    Ksoft : array_like
        Bulk modulus of soft sand model.

    Gsoft : array_like
        Shear modulus of soft sand model.
    """
    
    v = (3.*K-2.*G)/(2.*(3.*K+G))
    Khm = ((n**2. * (1.-phic)**2. * G**2 * P)/( 18.*np.pi**2 * (1.-v)**2))**(1./3.)
    Ghm = ((5.-4.*v)/(5.*(2.-v))) * ((3.*n**2 * (1.-phic)**2. * G**2. * P)/( 2*np.pi**2. * (1.-v)**2))**(1./3.)
    
    zhm = (Ghm/6.) * (9.*Khm + 8.*Ghm)/(Khm + 2.*Ghm)
    
    Ksoft = ((phi/phic)/(Khm + 4./3.* Ghm) + (1. - phi/phic)/(K + 4./3.*Ghm))**-1. - 4./3. * Ghm
    Gsoft = ((phi/phic)/(Ghm + zhm) + (1. - phi/phic)/(G + zhm))**-1. - zhm
    
    return Ksoft, Gsoft


#TODO Add other rock physics models


def gassmann(Kdry, Ks, Kf, phi):
    """Fluid substitution using Gassmann's Equation

    Parameters
    ----------
    Kdry : array_like
        Bulk modulus of dry rock.

    Ks : float
        Bulk modulus of mineral.

    Kf : float
        Bulk modulus of fluid.

    phi : array_like
        Porosity log.

    Returns
    -------
    Ksat : array_like
        Bulk modulus of saturated rock.
    """

    Ksat = np.zeros(len(phi))
    Ksat[0] = Ks
    Ksat [1:]= Ks * (phi[1:] * Kdry[1:] - ((1 + phi[1:]) * Kf * Kdry[1:] /Ks) + Kf) / ( (1-phi[1:]) * Kf + phi[1:] * Ks - Kf * Kdry[1:] / Ks)

    return Ksat


def rhob_gen(depth, phi, ow, rho, rho_oil, rho_water):
    """Pseudo-log of bulk density.
    
    Parameters
    ----------
    depth : array_like
        Depth log.

    phi : array_like
        Porosity log.

    ow : float
        Oil-water contact.

    rho : float
        Mineral density.

    rho_oil : float
        Oil density.

    rho_water : float
        Water density.

    Returns
    -------
    rhor : array_like
        Pseudo-log of bulk density.
    """

    rhor = np.empty(len(phi))

    rhor[(depth > ow)] = (1-phi[(depth > ow)]) * rho + phi[(depth > ow)] * rho_water
    rhor[(depth <= ow)] = (1-phi[(depth <= ow)]) * rho + phi[(depth <= ow)] * rho_oil

    return rhor
