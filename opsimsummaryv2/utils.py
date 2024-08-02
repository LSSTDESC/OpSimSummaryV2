"""This module contains usefull functions."""

import numpy as np

try:
    import geopandas as gpd
    import shapely.geometry as shp_geo
    import shapely.ops as shp_ops
    use_geopandas = True
except ImportError:
    use_geopandas = False


def format_poly(poly):
    """Fromat polygon that cross the 2 PI edge.

    Parameters
    ----------
    poly : shapely.polygons
        Polygon that represent fields

    Returns
    -------
    shapely.MultiPolygons
        Polygon that represent fields cutted on 2 PI edge.
    """
    if not use_geopandas:
        ModuleNotFoundError("Install geopandas library to use format_poly.")

    _SPHERE_LIMIT_HIGH_ = shp_geo.LineString(
        [[2 * np.pi, -np.pi / 2], [2 * np.pi, np.pi / 2]]
    )
    polydiv = gpd.GeoSeries(
        [p for p in shp_ops.polygonize(poly.boundary.union(_SPHERE_LIMIT_HIGH_))]
    )
    transl_mask = polydiv.boundary.bounds["maxx"] > 2 * np.pi
    polydiv[transl_mask] = polydiv[transl_mask].translate(-2 * np.pi)
    return shp_geo.MultiPolygon(polydiv.values)


def df_subdiviser(df, Nsub):
    """Divide a dataframe into a list of sub-dataframne

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas Dataframe to divide
    Nsub : int
        Number of subdivision

    Returns
    -------
    list(pandas.DataFrame)
        A list of dataframme subdivisions.
    """
    lensdf = len(df) // Nsub
    r = len(df) % Nsub
    subidx = np.arange(0, len(df), lensdf)
    subidx[:r] += np.arange(r)
    subidx[r:] += r
    sdfs = [df.iloc[i1:i2] for i1, i2 in zip(subidx[:-1], subidx[1:])]
    return sdfs


def host_joiner(survey_fields, host):
    """Use geopandas to match host in survey fields.

    Parameters
    ----------
    survey_fields : geopandas.GeoDataFrame
        Geodataframe describing survey fields
    host : pandas.DataFrame
        Dataframe that contains host informations

    Returns
    -------
    pandas.DataFrame
        Datafrane containing host that are in field with their correspind GROUPID
    """
        if not use_geopandas:
        ModuleNotFoundError("Install geopandas library to use host_joiner.")

    # Select host in circle
    host_pos = gpd.GeoDataFrame(
        index=host.index, geometry=gpd.points_from_xy(host.ra.values, host.dec.values)
    )
    grped_host = host_pos.sjoin(survey_fields, how="inner", predicate="intersects")

    # Create grp id
    grped_host = (
        grped_host.index_right.sort_values()
        .groupby(level=0)
        .apply(lambda x: np.random.choice(list(x)))
    )

    # Keep only hosts in fields
    survey_host = host.loc[grped_host.index]
    survey_host["GROUPID"] = grped_host
    return survey_host


def dataline(expMJD, ObsID, BAND, CCDgain, CCDnoise, SKYSIG, PSF, ZPT, ZPTNoise):
    """Write a SIMLIB data line

    Parameters
    ----------
    expMJD : float
        Date of observation in MJD
    ObsID : int
        ID of the observation
    BAND : str
        Band used for the observation
    CCDgain : float
        CCd gain
    CCDnoise : float
        CCD noise
    SKYSIG : float
        Sky noise error
    PSF : float
        Point spread function sigma
    ZPT : float
        Zero point
    ZPTNoise : float
        Zero point calibration error

    Returns
    -------
    str
        A SIMLIB LIB entry line
    """
    l = (
        "S: "
        f"{expMJD:5.4f} "
        f"{ObsID:10d}*2 "
        f"{BAND} "
        f"{CCDgain:5.2f} "  # CCD Gain
        f"{CCDnoise:5.2f} "  # CCD Noise
        f"{SKYSIG:6.2f} "  # SKYSIG
        f"{PSF:4.2f} "  # PSF1
        f"{0:4.2f} "  # PSF2
        f"{0:4.3f} "  # PSFRatio
        f"{ZPT:6.2f} "  # ZPTAVG
        f"{ZPTNoise:6.3f} "  # ZPTNoise
        f"{-99.:+7.3f} "
    )
    return l
