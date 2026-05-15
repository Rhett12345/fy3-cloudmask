"""Cloud mask test modules for different surface types and conditions."""

from .land_day import land_day_standard, land_day_coast, land_day_desert, land_day_desert_coast
from .land_nite import land_nite
from .ocean_day import ocean_day
from .ocean_nite import ocean_nite
from .polar_day import (
    polar_day_land, polar_day_coast, polar_day_desert,
    polar_day_desert_coast, polar_day_ocean, polar_day_snow,
)
from .polar_nite import polar_nite_land, polar_nite_ocean, polar_nite_snow
from .snow_tests import day_snow, nite_snow, antarctic_day
from .restoral import (
    chk_land_restoral, chk_land_nite_restoral, chk_coast_restoral,
    chk_sunglint_restoral, chk_shallow_water, chk_spatial_var,
    chk_cloud_adj, chk_thin_cirrus_ir, chk_shadow,
)

__all__ = [
    # Land tests
    'land_day_standard', 'land_day_coast', 'land_day_desert', 'land_day_desert_coast',
    'land_nite',
    # Ocean tests
    'ocean_day', 'ocean_nite',
    # Polar tests
    'polar_day_land', 'polar_day_coast', 'polar_day_desert', 'polar_day_desert_coast',
    'polar_day_ocean', 'polar_day_snow',
    'polar_nite_land', 'polar_nite_ocean', 'polar_nite_snow',
    # Snow tests
    'day_snow', 'nite_snow', 'antarctic_day',
    # Restoral tests
    'chk_land_restoral', 'chk_land_nite_restoral', 'chk_coast_restoral',
    'chk_sunglint_restoral', 'chk_shallow_water', 'chk_spatial_var',
    'chk_cloud_adj', 'chk_thin_cirrus_ir', 'chk_shadow',
]
