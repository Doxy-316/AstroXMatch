"""
Matching between two catalog of stars

@Author         : Doxy_316
@Date           : 22/02/2026
@last-update    : 24/02/2026

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


def cross_match(ref_catalog: QTable, target_catalog: QTable, radius: u.Quantity) -> QTable :
    """Matching stars between 2 catalogs

    Args:
        ref_catalog     (QTable): Table with colomns 'id', 'ra', 'dec'
        target_catalog  (QTable): Table with colomns 'id', 'ra', 'dec'
        radius (u.Quantity): max gap between two matched stars

    Returns:
        QTable: Table with index of matched stars and their distance
    """

    coords_ref = SkyCoord(ra=ref_catalog['ra'], dec=ref_catalog['dec'])
    coords_target = SkyCoord(ra=target_catalog['ra'], dec=target_catalog['dec'])

    idx, sep_2d, _ = coords_target.match_to_catalog_sky(coords_ref)

    mask = sep_2d < radius.to(u.deg)

    results = target_catalog[mask].copy()
    results['ref_index'] = idx[mask]
    results['angular_sep'] = sep_2d[mask]

    return results[['index', 'ref_index', 'angular_sep']]



##### -------------- #####
### --- Execussion --- ###
##### -------------- #####


if __name__ == '__main__':
    n_stars = 20
    ra_range = np.array((10, 15)) * u.deg
    dec_range = np.array((45, 50)) * u.deg
    sigma_ = 10 * u.arcsec
    max_radius = 10 * u.arcsec
    rng = np.random.default_rng(seed=42)

    qt_truth = generate_truth_catalog(n_stars, ra_range, dec_range)
    qt_data = add_astrometric_noise(qt_truth, sigma_)

    matched_table = cross_match(qt_truth, qt_data, max_radius)

    # print("\nqt_truth : \n", qt_truth)
    # print("\nqt_data : \n", qt_data)
    print("\nmatched_table : \n", matched_table)
