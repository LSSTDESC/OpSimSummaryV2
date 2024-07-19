import pandas as pd
import geopandas as gpd
import healpy as hp
import numpy as np
import sqlalchemy as sqla
from pathlib import Path
from sklearn.neighbors import BallTree

class OpSimSurvey:
    __LSST_FIELD_RADIUS__ = np.radians(1.75)
    __LSST_pixelSize__ = 0.2

    def __init__(self, db_path, table_name="observations", host_file=None, host_config={}):
        self.db_path = Path(db_path)
        self.sql_engine = self._get_sql_engine(db_path)
        self.opsimdf = self._get_df_from_sql(self.sql_engine)
        self.opsimdf.attrs['OpSimFile'] = self.db_path.name

        self.tree = BallTree(self.opsimdf [['_dec', '_ra']].values,
                             leaf_size=40,
                             metric='haversine')
        
        self.host = self._read_host_file(host_file, **host_config)
        
        self._hp_rep = None
        self._survey = None
        
        

    @staticmethod
    def _get_sql_engine(dbname):
        """Read a sql db Opsim output file.

        Parameters
        ----------
        dbname : str
            _description_

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        """        
        if not Path(dbname).exists():
            raise ValueError(f"{dbname} does not exists.")
        # Prepend the abs path with sqlite for use with sqlalchemy
        if not dbname.startswith('sqlite'):
            dbname = 'sqlite:///' + dbname
        print('Reading from database {}'.format(dbname))
        engine = sqla.create_engine(dbname, echo=False)
        return engine
    
    @staticmethod
    def _get_df_from_sql(sql_engine):
        df = pd.read_sql('observations', con=sql_engine)
        df['_ra'] = np.radians(df.fieldRA)
        df['_dec'] = np.radians(df.fieldDec)
        df.set_index('observationId', inplace=True)
        return df
    
    def _read_host_file(self, host_file, col_ra='ra', col_dec='dec', ra_dec_unit='radians'):
        if host_file is None:
            print('No host file.')
            return 
        
        print('Reading host from {}'.format(host_file))

        hostdf = pd.read_parquet(host_file)
        
        if ra_dec_unit == 'degrees':
            hostdf[col_ra] += 360 * (hostdf[col_ra] < 0) 
            hostdf.rename(columns={col_ra: 'RA_GAL', col_dec: 'DEC_GAL'}, inplace=True)
            hostdf['ra'] = np.radians(hostdf['RA_GAL'])
            hostdf['dec'] = np.radians(hostdf['DEC_GAL'])
        elif ra_dec_unit == 'radians':
            hostdf[col_ra] += 2 * np.pi * (hostdf[col_ra] < 0) 
            hostdf.rename(columns={col_ra: 'ra', col_dec: 'dec'}, inplace=True)
            hostdf['RA_GAL'] = np.degrees(hostdf['ra'])
            hostdf['DEC_GAL'] = np.degrees(hostdf['dec'])
        
        return gpd.GeoDataFrame(
            data=hostdf,
            geometry=gpd.points_from_xy(hostdf.ra.values, hostdf.dec.values))
        
    def compute_hp_rep(self, nside=256, minVisits=None, maxVisits=None):
        ipix = np.arange(hp.nside2npix(nside))
        hp_ra, hp_dec = np.radians(hp.pix2ang(nside, ipix, lonlat=True))
        
        self._hp_rep = pd.DataFrame(dict(ipix=ipix, hp_ra=hp_ra, hp_dec=hp_dec))

        # Compute number of visits
        self._hp_rep['n_visits'] = self.tree.query_radius(self.hp_rep[['hp_dec', 'hp_ra']].values, 
                                                            r=self.__LSST_FIELD_RADIUS__, 
                                                            count_only=True)
        
        # Apply mask if needed
        visits_mask = np.ones(len(self._hp_rep), dtype='bool')
        if minVisits is not None:
            visits_mask &= self._hp_rep['n_visits'] >= minVisits
            
        if maxVisits is not None:
            visits_mask &= self._hp_rep['n_visits'] <= maxVisits
        
            
        self._hp_rep = self._hp_rep[visits_mask]
        self._hp_rep.set_index('ipix', inplace=True)
        
        self._hp_rep.attrs['nside'] = nside
        
        print(f'Finished compute healpy representation, total number of fields : {len(self._hp_rep)}.')

    def plot_hp_rep(self, **hp_mollview_kwargs):
        ipix_map = np.zeros(hp.nside2npix(self.hp_rep.attrs['nside']))
        ipix_map[self.hp_rep.index.values] = self.hp_rep.n_visits.values
        return hp.mollview(map=ipix_map, unit='N visits', **hp_mollview_kwargs)
    
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
        
    def sample_survey(self, N_fields, random_seed=None):
        if N_fields > len(self.hp_rep):
            raise ValueError(f"N_feilds ({N_fields}) > survey fields ({len(self.hp_rep)})")
        
        seed = np.random.SeedSequence(random_seed)
        rng = np.random.default_rng(seed)
        
        self._survey = self.hp_rep.sample(n=N_fields, replace=False, random_state=rng)
        self._survey.reset_index(inplace=True)
        self._survey.attrs['nside'] = self.hp_rep.attrs['nside']
        self._survey.attrs['survey_area_deg'] = hp.nside2pixarea(self.hp_rep.attrs['nside'], degrees=True) * N_fields
        self._survey.attrs['survey_area_rad'] = hp.nside2pixarea(self.hp_rep.attrs['nside'], degrees=False) * N_fields
        self._survey.attrs['N_fields'] = N_fields

    def get_survey_obs(self, formatobs=True):
        obs_idx = self.tree.query_radius(self.survey[['hp_dec', 'hp_ra']],
                                              r=self.__LSST_FIELD_RADIUS__,
                                              count_only=False,
                                              return_distance=False)
        
        for idx in obs_idx:
            if formatobs:
                yield self.formatObs(self.opsimdf.iloc[idx])
            else:
                yield self.opsimdf.iloc[idx]
            
    def get_survey_host(self):
        # Cut in 2 circles that intersect edge limits (0, 2PI)
        _SPHERE_LIMIT_LOW_ = shp_geo.LineString([[0, -np.pi / 2],
                                                 [0, np.pi / 2]])

        _SPHERE_LIMIT_HIGH_ = shp_geo.LineString([[2 * np.pi, -np.pi / 2],
                                                  [2 * np.pi, np.pi / 2]])

        
        survey_fields = gpd.points_from_xy(self.survey.hp_ra.values,
                                           self.survey.hp_dec.values).buffer(self.__LSST_FIELD_RADIUS__)
        
        mask = survey_fields.intersects(_SPHERE_LIMIT_LOW_)
        survey_fields[mask] = survey_fields[mask].translate(2 * np.pi)
        mask |= survey_fields.intersects(_SPHERE_LIMIT_HIGH_)
        survey_fields[mask] = gpd.GeoSeries([op.utils.format_poly(p) for p in survey_fields[mask]])

        # Select host in circle
        grped_host = OpSimSurv.host[['geometry']].sjoin(
            gpd.GeoDataFrame(geometry=FieldPoints), how="inner", predicate="intersects"
        )
        
        # Create grp id
        grped_host = grped_host.index_right.sort_values().groupby(level=0).apply(lambda x:
                                                                    np.random.choice(list(x)))

        # Keep only hosts in fields
        survey_host = OpSimSurv.host.loc[grped_host.index]
        survey_host["GROUPID"] = grped_host
        
        return survey_host
    
    def formatObs(self, OpSimObs):
        
        sigPSF = OpSimObs['seeingFwhmEff'] / (2 * np.sqrt(2 * np.log(2)))
        
        # PSF
        PSF = sigPSF / self.__LSST_pixelSize__
        
        # Noise Effective Area
        NoiseArea = 4 * np.pi * sigPSF**2
        
        # ZPT
        dmag = OpSimObs['fiveSigmaDepth'] - OpSimObs['skyBrightness']
        
        ZPT = 2 * OpSimObs['fiveSigmaDepth'] - OpSimObs['skyBrightness']
        ZPT += 2.5 * np.log10(25 * NoiseArea) 
        ZPT += 2.5 * np.log10(1 + 10**(-0.4 * dmag) / NoiseArea)
        
        # SKYSIG
        dskymag = OpSimObs['skyBrightness'] - ZPT
        SKYSIG = np.sqrt(10**(-0.4 * dskymag) * self.__LSST_pixelSize__**2)
        
        return pd.DataFrame({'expMJD': OpSimObs['observationStartMJD'],
                             'PSF': PSF,
                             'ZPT': ZPT,
                             'SKYSIG': SKYSIG,
                             'BAND': OpSimObs['filter']})
