"""
Generation of test datas.

@Author         : Doxy_316
@Date           : 19/02/2026
@last update    : 20/02/2026

Functions :
    - generate_truth_catalog()
    - add_astrometric_noise()
"""

import numpy as np
from astropy import units as u
from astropy.table import QTable

# from astropy.coordinates import sky_coordinate


def generate_truth_catalog(n_sources: int, ra_range: tuple, dec_range: tuple) -> QTable:
    """Generate a QTable containing sky coordinates in a certain range

    Args:
        n_sources (int): number of coordinates
        ra_range (tuple): Limites of the coordinates (ra_min, ra_max) - deg
        dec_range (tuple): Limites of the coordinates (dec_min, dec_max) - deg

    Returns:
        QTable: Catalog of coordinates containing : (index, ra, dec)
    """

    qt = QTable()
    qt["index"] = np.linspace(0, n_sources - 1, n_sources, dtype=int)

    ra = np.random.uniform(*ra_range, n_sources)
    qt["ra"] = ra * u.deg

    # Uniform in sin(dec) :
    smin, smax = np.sin(np.deg2rad(dec_range))
    dec_temp = np.random.uniform(smin, smax, n_sources)
    dec = np.rad2deg(np.arcsin(dec_temp))
    qt['dec'] = dec * u.deg

    return qt


def add_astrometric_noise(truth_table: QTable, sigma: u.arcsec) -> QTable:
    """Takes a table and simulate an observation (Gaussian)

    Args:
        truth_table (QTable): Contain sky coordinates
        sigma (u.arcsec): Arg for the gaussian noise

    Returns:
        QTable: Table of coordinates with Gaussian noise
    """

    qt = QTable()
    qt["index"] = truth_table["index"]

    # Convergence of meridians :
    sigma_ra = sigma / np.cos(truth_table['dec'])
    ra_noise = np.random.normal(truth_table["ra"], sigma_ra.to(u.deg))
    qt["ra"] = ra_noise * u.deg

    dec_noise = np.random.normal(truth_table["dec"], sigma.to(u.deg))
    qt["dec"] = dec_noise * u.deg

    return qt





# TEST :

NBR = 20
ra_ = (10, 15)
dec_ = (45, 50)
sigma_ = 100 * u.arcsec

qt_1 = generate_truth_catalog(NBR, ra_, dec_)
qt_2 = add_astrometric_noise(qt_1, sigma_)

print("qt1 : \n", qt_1)
print("qt2 : \n", qt_2)
