import geopandas as gpd
import numpy as np
import shapely.geometry as shp_geo
import shapely.ops as shp_ops


def format_poly(poly):
    """Format polygons that crosses sphere limit (0, 2 PI)

    Args:
        poly (_type_): _description_

    Returns:
        _type_: _description_
    """
    _SPHERE_LIMIT_HIGH_ = shp_geo.LineString(
        [[2 * np.pi, -np.pi / 2], [2 * np.pi, np.pi / 2]]
    )
    polydiv = gpd.GeoSeries(
        [p for p in shp_ops.polygonize(poly.boundary.union(_SPHERE_LIMIT_HIGH_))]
    )
    transl_mask = polydiv.boundary.bounds["maxx"] > 2 * np.pi
    polydiv[transl_mask] = polydiv[transl_mask].translate(-2 * np.pi)
    return shp_geo.MultiPolygon(polydiv.values)
    
