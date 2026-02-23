"""
Matching between two catalog of stars

@Author         : Doxy_316
@Date           : 22/02/2026
@last-update    : 23/02/2026

functions :
    - cross_match()
"""


import numpy as np
import astropy.units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord
from catalog_utils import generate_truth_catalog, add_astrometric_noise



##### -------------- #####
### ---- Functions --- ###
##### -------------- #####


def cross_match(catalog_a: QTable, catalog_b: QTable, radius: u.Quantity) -> QTable :
    """Matching stars between 2 catalogs with a maximum distance

    Args:
        catalog_a (QTable): Table with colomns 'id', 'ra', 'dec' (Truth)
        catalog_b (QTable): Table with colomns 'id', 'ra', 'dec' (Data)
        radius (u.Quantity): max gap between two matched stars

    Returns:
        QTable: Table with index of matched stars and their distance
    """

    ra_truth = catalog_a['ra'].to(u.deg)
    dec_truth = catalog_a['dec'].to(u.deg)
    c_truth = SkyCoord(ra_truth, dec_truth, frame='icrs')

    ra_data = catalog_b['ra'].to(u.deg)
    dec_data = catalog_b['dec'].to(u.deg)
    c_data = SkyCoord(ra_data, dec_data, frame='icrs')

    match_result = c_truth.match_to_catalog_sky(c_data, nthneighbor=1)
    angular_separation = match_result[1]

    dist_mask = (angular_separation < radius.to(u.deg))

    qt = QTable()
    idx = catalog_a['index']
    qt['index'] = idx[dist_mask]
    qt['angular_separation'] = angular_separation[dist_mask]

    return qt



##### -------------- #####
### --- Execussion --- ###
##### -------------- #####


if __name__ == '__main__':
    NBR = 20
    ra_ = np.array((10, 15)) * u.deg
    dec_ = np.array((45, 50)) * u.deg
    sigma_ = 100 * u.arcsec
    max_radius = 70 * u.arcsec
    rng = np.random.default_rng(seed=42)

    qt_truth = generate_truth_catalog(NBR, ra_, dec_)
    qt_data = add_astrometric_noise(qt_truth, sigma_)

    # qt_matched = 
    qt_matched = cross_match(qt_data, qt_truth, max_radius)

    # print("\nqt_truth : \n", qt_truth)
    # print("\nqt_data : \n", qt_data)
    print("\nqt_matched : \n", qt_matched)