"""
Matching between two catalog of stars

@Author         : Doxy_316
@Date           : 22/02/2026
@last-update    : 22/02/2026

functions :
    - cross_match()
"""


import numpy as np
from astropy import units as u
from astropy.table import QTable
from astropy.coordinates import sky_coordinate
from catalog_utils import generate_truth_catalog, add_astrometric_noise


##### -------------- #####
### ---- Functions --- ###
##### -------------- #####


def cross_match(catalog_a: QTable, catalog_b: QTable, radius: u.Quantity) -> QTable :
    """Matching stars between 2 catalogs with a maximum gap

    Args:
        catalog_a (QTable): Catalog of coordinates (Data)
        catalog_b (QTable): Catalog of coordinates (Truth)
        radius (u.Quantity): max gap between two matched stars

    Returns:
        QTable: Table with index of matched stars and their distance
    """
    ...
    pass



##### -------------- #####
### --- Execussion --- ###
##### -------------- #####


if __name__ == '__main__':
    NBR = 20
    ra_ = np.array((10, 15)) * u.deg
    dec_ = np.array((45, 50)) * u.deg
    sigma_ = 100 * u.arcsec
    match_radius = 50 * u.arcsec

    qt_truth = generate_truth_catalog(NBR, ra_, dec_)
    qt_data = add_astrometric_noise(qt_truth, sigma_)

    qt_matched = cross_match(qt_data, qt_truth, match_radius)

    print("\nqt_truth : \n", qt_truth)
    print("\nqt_data : \n", qt_data)
    print("\nqt_matched : \n", qt_matched)