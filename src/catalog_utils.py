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
import astropy.units as u 
import astropy.coordinates.sky_coordinate

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
