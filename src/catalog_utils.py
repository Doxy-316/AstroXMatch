"""
Generation of test datas.

@Author         : Doxy_316
@Date           : 19/02/2026
@last update    : 21/02/2026

Functions :
    - generate_truth_catalog()
    - add_astrometric_noise()
"""

import numpy as np
from astropy import units as u
from astropy.table import QTable

rng = np.random.default_rng(seed=42)

def generate_truth_catalog(n_sources: int, ra_range: u.Quantity, dec_range: u.Quantity) -> QTable:
    """Generate a QTable containing sky coordinates in a certain range

    Args:
        n_sources (int): number of coordinates
        ra_range (u.Quantity): Limites of the coordinates (ra_min, ra_max) - deg
        dec_range (u.Quantity): Limites of the coordinates (dec_min, dec_max) - deg

    Returns:
        QTable: Catalog of coordinates containing : (index, ra, dec)
    """

    qt = QTable()
    qt["index"] = np.arange(n_sources)

    # ra :
    ra_min = ra_range.to(u.deg).value[0]
    ra_max = ra_range.to(u.deg).value[1]
    ra = rng.uniform(ra_min, ra_max, n_sources)
    qt["ra"] = ra * u.deg

    # dec :
    dec_min = dec_range.to(u.deg).value[0]
    dec_max = dec_range.to(u.deg).value[1]
    # Uniform in sin(dec) :
    smin = np.sin(np.deg2rad(dec_min))
    smax = np.sin(np.deg2rad(dec_max))
    dec_temp = rng.uniform(smin, smax, n_sources)
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
    sigma_ra = sigma / np.cos(truth_table['dec'].to(u.rad))
    ra_noise = rng.normal(truth_table["ra"], sigma_ra.to(u.deg))
    qt["ra"] = ra_noise * u.deg

    dec_noise = rng.normal(truth_table["dec"], sigma.to(u.deg))
    qt["dec"] = dec_noise * u.deg

    return qt





# TEST :

if __name__ == '__main__':
    NBR = 20
    ra_ = np.array((10, 15)) * u.deg
    dec_ = np.array((45, 50)) * u.deg
    sigma_ = 100 * u.arcsec

    qt_1 = generate_truth_catalog(NBR, ra_, dec_)
    qt_2 = add_astrometric_noise(qt_1, sigma_)

    print("qt1 : \n", qt_1)
    print("qt2 : \n", qt_2)
