import numpy as np


def sw_gen(depth, ow, so):
    """Generate pseudo-log of water saturation.
    
    Parameters
    ----------
    depth : array_like
        Depth log.

    ow : float
        Oil-water contact depth.

    so : float
        Oil saturation.

    Returns
    -------
    sw: array_like
        Sw pseudo-log.
    """

    sw = np.empty(len(depth))

    oil = (depth <= ow)
    water = (depth > ow)

    sw[oil] = 1 - so
    sw[water] = 1.0

    return sw