"""Main module to read OpSim output database."""

import time
from functools import partial
from multiprocessing import Pool
from pathlib import Path
import pandas as pd

try:
    import geopandas as gpd
    import shapely.geometry as shp_geo
    import shapely.affinity as shp_aff

    use_geopandas = True
except ImportError:
    use_geopandas = False

import healpy as hp
import numpy as np
import sqlalchemy as sqla

from sklearn.neighbors import BallTree
from astropy.time import Time
from . import utils as ut


class OpSimSurvey:
    """
    A class to manipulate OpSim db data and turn them into simulation inputs.

    Attributes
    ----------
    db_path : pathlib.Path
        Path to the database file
    sql_engine : sqlalchemy.engine.base.Engine
        The sqlalchemy engine link to the db.
    opsimdf : pandas.DataFrame
        Dataframe of the simulated observations.
    tree : sklearn.neighbors.BallTree
        A BallTree used to select observations on the healpy representation of the survey.
    host : pandas.DataFrame
        Dataframe of the hosts.
    hp_rep : pandas.DataFrame
        The healpy representation oif the survey
    survey : pandas.DataFrame
        The healpy representation oif the survey
    survey_hosts : pandas.DataFrame
        The hosts that are inside the survey.
    """

    __LSST_FIELD_RADIUS__ = np.radians(1.75)  # LSST Field Radius in radians
    __LSST_pixelSize__ = 0.2  # LSST Pixel size in arcsec^-1
    __LSST_DDF_TAGS__ = {
        'TAGS': np.array(['ELAISS1', 'XMM_LSS', 'ECDFS', 'COSMOS', 'EDFS_a', 'EDFS_b']), # Name of DDF fields
        'RA': np.array([9.45, 35.71, 53.12, 150.10, 58.90, 63.60]), # RA coords of DDF fields
        'Dec': np.array([-44.00, -4.75, -28.10, 2.18, -49.32, -47.60]) # Dec coords of DDF fields
        }
    
    def __init__(
        self,
        db_path,
        table_name="observations",
        MJDrange=None,
        host_file=None,
        host_config={},
    ):
        """Construct an OpSimSurvey object from an OpSim DataBase.

        Parameters
        ----------
        db_path : str
            The path to the Opsim db file
        table_name : str, optional
            Name of the observations table in the OpSIm db file, by default "observations"
        MJDrange : (int, int) or (str,str), optional
            Min and Max date to query if float assumed to be mjd, by default None
        host_file : str, optional
            Path to a parquet file containg hosts, by default None
        host_config : dict, optional
            Configuration for reading host file, by default {}
        """
        self.db_path = Path(db_path)
        self.sql_engine = self._get_sql_engine(db_path)
        self.opsimdf = self._get_df_from_sql(self.sql_engine, MJDrange=MJDrange)
        self.opsimdf.attrs["OpSimFile"] = self.db_path.name

        self.tree = BallTree(
            self.opsimdf[["_dec", "_ra"]].values, leaf_size=50, metric="haversine"
        )

        self.host = self._read_host_file(host_file, **host_config)

        self._hp_rep = None
        self._survey = None

    @staticmethod
    def _get_sql_engine(dbname):
        """Read a sql db Opsim output file.

        Parameters
        ----------
        dbname : str
            Path to sql file.

        Returns
        -------
        sqlalchemy.engine.base.Engine
            The sqlalchemy engine link to the db.
        """
        if not Path(dbname).exists():
            raise ValueError(f"{dbname} does not exists.")
        # Prepend the abs path with sqlite for use with sqlalchemy
        if not dbname.startswith("sqlite"):
            dbname = "sqlite:///" + dbname
        print("Reading from database {}".format(dbname))
        engine = sqla.create_engine(dbname, echo=False)
        return engine

    @staticmethod
    def _get_df_from_sql(sql_engine, MJDrange=None):
        """Load data from the db file.

        Parameters
        ----------
        sql_engine : sqlalchemy.engine.base.Engine
            The sqlalchemy engine link to the db.
        MJDrange :( , ) str or float, optional
            Min and Max date to query if float assumed to be mjd, by default None

        Returns
        -------
        pandas.DataFrame
            Dataframe of the simulated observations.
        """
        tstart = time.time()
        query = "SELECT * FROM observations"
        if MJDrange is not None:
            if isinstance(MJDrange, str):
                time_format = None
            else:
                time_format = "mjd"
            MJDrange = Time(MJDrange, format=time_format)
            query += f" WHERE observationStartMJD > {MJDrange.mjd[0]} AND observationStartMJD < {MJDrange[1].mjd}"

        df = pd.read_sql(query, con=sql_engine)
        df["_ra"] = np.radians(df.fieldRA)
        df["_dec"] = np.radians(df.fieldDec)
        df.set_index("observationId", inplace=True)
        print(f"Read N = {len(df)} observations in {time.time() - tstart:.2f} seconds.")
        return df.sort_values(by="observationStartMJD")

    def _read_host_file(
        self,
        host_file,
        col_ra="ra",
        col_dec="dec",
        ra_dec_unit="radians",
        wgt_map=None,
        add_SNMAGSHIFT=False,
    ):
        """Read a parquet file containing hosts.

        Parameters
        ----------
        host_file : str
            Path to the parquet file
        col_ra : str, optional
            Key of column containing RA, by default 'ra'
        col_dec : str, optional
            Key of column containing Dec, by default 'dec'
        ra_dec_unit : str, optional
            Unit of ra_dec (radians or degrees), by default 'radians'

        Returns
        -------
        pandas.DataFrame
            Dataframe of the hosts.
        """
        if host_file is None:
            print("No host file.")
            return None

        print("Reading host from {}".format(host_file))
        hostdf = pd.read_parquet(host_file)

        if wgt_map is not None:
            print(f"Reading and applying HOST WGT MAP from {wgt_map}")
            var_names, wgt_map = ut.read_SNANA_WGTMAP(wgt_map)
            if len(var_names) > 1:
                raise NotImplementedError("HOST RESAMPLING ONLY WORK FOR 1 VARIABLES")
            keep_index = ut.host_resampler(
                wgt_map[var_names[0]],
                wgt_map["WGT"],
                hostdf.index.values,
                hostdf[var_names[0]].values,
            )

            hostdf = hostdf.loc[keep_index]

            if add_SNMAGSHIFT and "SNMAGSHIFT" in wgt_map:
                snmagshift = np.zeros(len(hostdf))
                for i in range(len(wgt_map["WGT"]) - 1):
                    mask = np.ones(len(hostdf), dtype=bool)
                    for v in var_names:
                        mask &= hostdf[v].between(wgt_map[v][i], wgt_map[v][i + 1])
                    snmagshift[mask] = wgt_map["SNMAGSHIFT"][i]
                hostdf["SNMAGSHIFT"] = snmagshift

        if ra_dec_unit == "degrees":
            hostdf[col_ra] += 360 * (hostdf[col_ra] < 0)
            hostdf.rename(columns={col_ra: "RA_GAL", col_dec: "DEC_GAL"}, inplace=True)
            hostdf["ra"] = np.radians(hostdf["RA_GAL"])
            hostdf["dec"] = np.radians(hostdf["DEC_GAL"])
        elif ra_dec_unit == "radians":
            hostdf[col_ra] += 2 * np.pi * (hostdf[col_ra] < 0)
            hostdf.rename(columns={col_ra: "ra", col_dec: "dec"}, inplace=True)
            hostdf["RA_GAL"] = np.degrees(hostdf["ra"])
            hostdf["DEC_GAL"] = np.degrees(hostdf["dec"])
        hostdf.attrs["file"] = host_file
        return hostdf

    def compute_hp_rep(
        self,
        nside=256,
        minVisits=500,
        maxVisits=100_000,
        ddf_nobs_thresh=1100,
        add_ddf_tag=False,
        angle_sep_tol=10
    ):
        """Compute a healpy version of the survey.

        Parameters
        ----------
        nside : int, optional
            Healpy nside, by default 256
        minVisits : int, optional
            Minimum number of observations required in the survey, by default 500
        maxVisits : int, optional
            Maximum number of observations required in the survey, by default 100000
        """
        ipix = np.arange(hp.nside2npix(nside))
        hp_ra, hp_dec = np.radians(hp.pix2ang(nside, ipix, lonlat=True))

        self._hp_rep = pd.DataFrame(dict(ipix=ipix, hp_ra=hp_ra, hp_dec=hp_dec))

        # Compute number of visits
        self._hp_rep["n_visits"] = self.tree.query_radius(
            self.hp_rep[["hp_dec", "hp_ra"]].values,
            r=self.__LSST_FIELD_RADIUS__,
            count_only=True,
        )

        # Apply mask if needed
        visits_mask = np.ones(len(self._hp_rep), dtype="bool")
        if minVisits is not None:
            visits_mask &= self._hp_rep["n_visits"] >= minVisits

        if maxVisits is not None:
            visits_mask &= self._hp_rep["n_visits"] <= maxVisits

        self._hp_rep = self._hp_rep[visits_mask]
        self._hp_rep.set_index("ipix", inplace=True)

        self._hp_rep.attrs["nside"] = nside
        self._hp_rep.attrs["survey_area_deg"] = hp.nside2pixarea(
            self.hp_rep.attrs["nside"], degrees=True
        ) * len(self._hp_rep)
        self._hp_rep.attrs["survey_area_rad"] = hp.nside2pixarea(
            self.hp_rep.attrs["nside"], degrees=False
        ) * len(self._hp_rep)

        field_labels = np.empty_like(self._hp_rep["n_visits"], dtype="U20")
        ddf_mask = self._hp_rep["n_visits"] >= ddf_nobs_thresh
        
        # WFD naming
        field_labels[~ddf_mask] = 'WFD'
        
        # DDF naming
        if add_ddf_tag:
            angle_sep_to_field = np.empty((len(self.hp_rep["hp_ra"]), 
                                           len(self.__LSST_DDF_TAGS__['TAGS'])),
                                          dtype='float32')
            
            for i, (R, D) in enumerate(zip(self.__LSST_DDF_TAGS__['RA'], self.__LSST_DDF_TAGS__['Dec'])):
                angle_sep_to_field.T[i] = ut.compute_angle_sep(
                        self.hp_rep["hp_ra"].values, 
                        self.hp_rep["hp_dec"].values,
                        np.deg2rad(R), np.deg2rad(D)
                    )
            angle_sep_argmins = np.argmin(angle_sep_to_field, axis=1)
            angle_sep_min =  np.min(angle_sep_to_field, axis=1)
            ddf_tags = np.char.add('DDF_' , self.__LSST_DDF_TAGS__['TAGS'][angle_sep_argmins])
            ddf_tags[np.degrees(angle_sep_min) > angle_sep_tol] = 'DDF_NONE'
        else:
            ddf_tags = 'DDF'
        
        field_labels[ddf_mask] = ddf_tags

        self._hp_rep["field_label"] = field_labels
        
        print(
            f"Finished compute healpy representation, total number of fields : {len(self._hp_rep)}."
        )

    def plot_hp_rep(self, **hp_mollview_kwargs):
        """Plot the healpy mollview of the survey representation.

        Returns
        -------
        matplotlib.Figure
            The mollview figure
        """

        ipix_map = np.zeros(hp.nside2npix(self.hp_rep.attrs["nside"]))
        ipix_map[self.hp_rep.index.values] = self.hp_rep.n_visits.values
        return hp.mollview(map=ipix_map, unit="N visits", **hp_mollview_kwargs)

    def sample_survey(self, N_fields, random_seed=None, nworkers=10):
        """Sample Nfields inside the survey's healpy representation.

        Parameters
        ----------
        N_fields : int
            Number of fields to sample
        random_seed : int or numpy.random.SeedSquence, optional
            The random  seed used to sample the fields, by default None
        nworkers : int, optional
            Number of cores used to run multiprocessing on host matching, by default 10

        Notes
        ------
        Random seed only apply on field sampling.
        """
        if N_fields > len(self.hp_rep):
            raise ValueError(
                f"N_fields ({N_fields}) > survey fields ({len(self.hp_rep)})"
            )

        seed = np.random.SeedSequence(random_seed)
        rng = np.random.default_rng(seed)

        self._survey = self.hp_rep.sample(n=N_fields, replace=False, random_state=rng)
        self._survey.reset_index(inplace=True)
        self._survey.attrs = self.hp_rep.attrs.copy()
        self._survey.attrs["N_fields"] = N_fields

        if self.host is not None:
            print("Compute survey hosts")
            self._survey_hosts = self.get_survey_hosts(nworkers=nworkers)

    def get_obs_from_coords(self, ra, dec, is_deg=True, formatobs=False, keep_keys=[]):
        """Get observations at ra, dec coordinates.

        Parameters
        ----------
        ra : numpy.ndarray(float)
            RA coordinate
        dec : numpy.ndarray(float)
            Dec coordinate
        is_deg : bool, optional
            is RA, Dec given in degrees, by default True
        formatobs : bool, optional
            format obs for simulation, by default False
        keep_keys : list(str)
            List of keys to keep in addition to formatted quantities

        Yields
        ------
        pandas.DatFrame
            Dataframes of observations.
        """
        if is_deg:
            ra = np.radians(ra)
            dec = np.radians(dec)

        obs_idx = self.tree.query_radius(
            np.array([dec, ra]).T,
            r=self.__LSST_FIELD_RADIUS__,
            count_only=False,
            return_distance=False,
        )
        for idx in obs_idx:
            if formatobs:
                yield self.formatObs(self.opsimdf.iloc[idx], keep_keys=keep_keys)
            else:
                yield self.opsimdf.iloc[idx]

    def get_survey_obs(self, formatobs=True, keep_keys=[]):
        """Get survey observations.

        Parameters
        ----------
        formatobs : bool, optional
            Format the observation to get only quantities of interest for simulation, by default True
        keep_keys : list(str)
            List of keys to keep in addition to formatted quantities

        Yields
        ------
        pandas.DatFrame
            Dataframes of observations.
        """
        return self.get_obs_from_coords(
            *self.survey[["hp_ra", "hp_dec"]].values.T,
            is_deg=False,
            formatobs=formatobs,
            keep_keys=keep_keys,
        )

    def get_survey_hosts(self, nworkers=10):
        """Get survey hosts.

        Parameters
        ----------
        nworkers : int, optional
            Number of cores used to run multiprocessing, by default 10

        Returns
        -------
        pandas.DataFrame
            Dataframe of host inside the survey, matched to their field indicated by GROUPID.
        """
        if not use_geopandas:
            raise ModuleNotFoundError("Install geopandas library to use host matching.")

        if self.host is None:
            raise ValueError("No host file set.")

        # Cut in 2 circles that intersect edge limits (0, 2PI)
        _SPHERE_LIMIT_LOW_ = shp_geo.LineString([[0, -np.pi / 2], [0, np.pi / 2]])

        _SPHERE_LIMIT_HIGH_ = shp_geo.LineString(
            [[2 * np.pi, -np.pi / 2], [2 * np.pi, np.pi / 2]]
        )

        survey_fields = gpd.GeoDataFrame(
            index=self.survey.index,
            geometry=gpd.points_from_xy(
                self.survey.hp_ra.values, self.survey.hp_dec.values
            ).buffer(self.__LSST_FIELD_RADIUS__),
        )

        # scale for dec dependance
        survey_fields = survey_fields.map(
            lambda x: shp_aff.scale(
                x, xfact=np.sqrt(2 / (1 + np.cos(2 * x.centroid.xy[1][0])))
            )
        )

        # mask for edge effect
        mask = survey_fields.intersects(_SPHERE_LIMIT_LOW_)
        survey_fields[mask] = survey_fields[mask].translate(2 * np.pi)
        mask |= survey_fields.intersects(_SPHERE_LIMIT_HIGH_)
        survey_fields.loc[mask, "geometry"] = gpd.GeoSeries(
            data=[ut.format_poly(p) for p in survey_fields[mask].geometry],
            index=survey_fields[mask].index,
        )

        host_joiner = partial(ut.host_joiner, survey_fields)

        sdfs = ut.df_subdiviser(self.host, Nsub=nworkers)

        with Pool(nworkers) as p:
            res = p.map(host_joiner, sdfs)

        survey_host = pd.concat(res)

        return survey_host

    def formatObs(self, OpSimObs, keep_keys=[]):
        """Format function to get quantities of interest for simulation.

        Parameters
        ----------
        OpSimObs : pandas.DataFrame
            Dataframe of OpSim output
        keep_keys : list(str)
            List of keys to keep in addition to formatted quantities

        Returns
        -------
        pandas.DataFrame
            Dataframe that contains quantities for simulation.

        Notes
        -----
        Quantities are obtained following arxiv:1905.02887

        Details of the calculs:

        The :math:`\sigma_\mathrm{PSF}` in units of :math:`\mathrm{arcsec}^{-1}` is obtained as

        .. math::
            \sigma_\mathrm{PSF} = \\frac{\mathbf{seeingFwhmEff}}{2\sqrt{2\ln2}}.

        The :math:`\mathbf{PSF}` is computed in units of :math:`\mathrm{pixel}^{-1}` as

        .. math::
            \mathbf{PSF} = \\frac{\sigma_\mathrm{PSF}}{\mathrm{PixelSize}},

        where :math:`\mathrm{PixelSize}` is the pixel size in :math:`\mathrm{arcsec}^{-1}`.

        The zero-point :math:`\mathbf{ZPT}` is computed as

        .. math::
            \mathbf{ZPT} = 2 m_5 - m_\mathrm{sky} + 2.5 \log\left[25 A \left(1 + A^{-1}10^{-0.4(m_5-m_\mathrm{sky})}\\right)\\right],

        where :math:`m_5 = \mathbf{fiveSigmaDepth}`, :math:`m_\mathrm{sky} = \mathbf{skyBrightness}` and :math:`A` is the noise equivalent area given by
        :math:`A = 4 \pi \sigma_\mathrm{PSF}^2`.

        The sky noise :math:`\mathbf{SKYSIG}` in unit of :math:`\mathrm{ADU}.\mathrm{pixel}^{-1}` is computed such as

        .. math::
            \mathbf{SKYSIG}^2 = 10^{-0.4\left(m_\mathrm{sky} - \mathbf{ZPT}\\right)} \\times \mathrm{PixelSize}^2.
        """
        sigPSF = OpSimObs["seeingFwhmEff"] / (2 * np.sqrt(2 * np.log(2)))

        # PSF in pixels^{-1}
        PSF = sigPSF / self.__LSST_pixelSize__

        # Noise Effective Area in arcsec^{-2}
        NoiseArea = 4 * np.pi * sigPSF**2

        # ZPT
        dmag = OpSimObs["fiveSigmaDepth"] - OpSimObs["skyBrightness"]

        ZPT = 2 * OpSimObs["fiveSigmaDepth"] - OpSimObs["skyBrightness"]
        ZPT += 2.5 * np.log10(25 * NoiseArea)
        ZPT += 2.5 * np.log10(1 + 10 ** (-0.4 * dmag) / NoiseArea)

        # SKYSIG
        dskymag = OpSimObs["skyBrightness"] - ZPT
        SKYSIG = np.sqrt(10 ** (-0.4 * dskymag) * self.__LSST_pixelSize__**2)

        formatobs_df = pd.DataFrame(
            {
                "expMJD": OpSimObs["observationStartMJD"],
                "PSF": PSF,
                "ZPT": ZPT,
                "SKYSIG": SKYSIG,
                "BAND": OpSimObs["filter"],
                **{k: OpSimObs[k] for k in keep_keys},
            }
        ).reset_index(names="ObsID")

        return formatobs_df

    @property
    def hp_rep(self):
        if self._hp_rep is None:
            raise ValueError("hp_rep not set. Please run compute_hpix_survey before.")
        return self._hp_rep

    @property
    def survey(self):
        if self._survey is None:
            raise ValueError("survey not set. Please run sample_survey() before.")
        return self._survey

    @property
    def survey_hosts(self):
        if self.host is None:
            return None
        return self._survey_hosts
