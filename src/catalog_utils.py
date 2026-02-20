"""
Generation of test datas.

@Author         : Nathan de Porcaro
@Date           : 19/02/2026
@last update    : 19/02/2026

Functions :
    - generate_truth_catalog()
    - add_astrometric_noise()
"""

import os
import sys
import numpy as np
from astropy import units as u
from astropy.table import QTable
# from astropy.coordinates import sky_coordinate

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


def generate_truth_catalog(n_sources: int, ra_range: tuple, dec_range: tuple) -> QTable :
    """Generate a QTable containing sky coordinates in a certain range

    Args:
        n_sources (int): number of coordinates
        ra_range (tuple): Limites of the coordinates (ra_min, ra_max) - deg
        dec_range (tuple): Limites of the coordinates (dec_min, dec_max) - deg

    Returns:
        QTable: Catalog of coordinates containing : (index, ra, dec)
    """

    qt = QTable()
    qt['index'] = np.linspace(0, n_sources - 1, n_sources, dtype=int)

    ra = np.random.uniform(*ra_range, n_sources)
    qt['ra'] = ra * u.deg

    dec_temp = np.arcsin(2*np.random.uniform(size=n_sources) - 1) + (np.pi / 2.0)
    len_dec = dec_range[1] - dec_range[0]
    dec = (len_dec / np.pi) * dec_temp + dec_range[0]
    qt['dec'] = dec * u.deg

    return qt




# EXECUTION :

nb = 20
ra_ = (10, 15)
dec_ = (45, 50)

qt_ = generate_truth_catalog(nb, ra_, dec_)
print(qt_)
