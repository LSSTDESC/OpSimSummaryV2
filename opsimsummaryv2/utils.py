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

def read_SNANA_WGTMAP(file):
    """Read a SNANA HOSTLIB WGTMAP

    Parameters
    ----------
    file : str
        path to WGTMAP

    Returns
    -------
    list(str), dict
        list of varnames, WGTMAP as a dict

    Raises
    ------
    ValueError
        No WGT key
    """    
    f = open(file, "r")

    data_starts = False
    for l in f:
        if 'VARNAMES_WGTMAP' in l:
            var_names = l.split()[1:]
            data = {k: [] for k in var_names}
            data_starts = True
        elif data_starts and 'WGT:' in l:
            for k, val in zip(var_names, l.split()[1:]):
                data[k].append(float(val))
    f.close()
    for k in data:
        data[k] = np.array(data[k])
        
    if 'WGT' in var_names:
        var_names.remove('WGT')
    else:
        raise ValueError("WGTMAP require 'WGT' key")
    if 'SNMAGSHIFT' in var_names:
        var_names.remove('SNMAGSHIFT')
    return var_names, data

def host_resampler(wgt_map_VAR, wgt_map_WGT, index, values, cdf_cut=0.95):
    """Resample host according to a WGTMAP

    Parameters
    ----------
    wgt_map_VAR : numpy.array(float)
        Variable of the WGTMAP
    wgt_map_WGT : numpy.array(float)
        Weight values of the WGTMAP
    index : numpy.array(int)
        Index of hosts
    values : numpy.array(float)
        Values of the Variable for host
    cdf_cut : float, optional
        A cut on the cdf to adjust , by default 0.95

    Returns
    -------
    _type_
        _description_
    """    
    count, edges = np.histogram(values, bins='rice')
    medges = (edges[1:] + edges[:-1]) * 0.5
    
    prob_select = count * np.interp(medges, wgt_map_VAR, wgt_map_WGT)
    cdf = np.cumsum(prob_select) 
    cdf /= cdf[-1]

    argmax = np.argmax(prob_select)
    count_max = count[argmax]
    prob_select_max = prob_select[argmax]
    
    N_to_draw = np.rint(prob_select / prob_select_max * count_max)
    cdf_mask = (cdf > 1.0 - cdf_cut) & (cdf < cdf_cut)
    correction_coeff = np.min(count[cdf_mask] / N_to_draw[cdf_mask])
    
    N_to_draw = np.rint(N_to_draw * correction_coeff).astype('int')
    
    wgt_values = np.interp(values, wgt_map_VAR, wgt_map_WGT)

    keep_idx = []
    for i, N in enumerate(N_to_draw):
        edmin, edmax = edges[i], edges[i + 1]
        mask = (values >= edmin) & (values < edmax)
        
        sdf_index = index[mask]

        if len(sdf_index) <  N:
            keep_idx.extend(sdf_index)
        elif N > 0:    
            wgt = wgt_values[mask]
            wgt /= np.sum(wgt)
            keep_idx.extend(np.random.choice(sdf_index, size=N, replace=False, p=wgt))
    return keep_idx
