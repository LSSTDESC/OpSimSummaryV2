import pandas as pd
import geopandas as gpd
import healpy as hp
import numpy as np
import sqlalchemy as sqla
from pathlib import Path
from sklearn.neighbors import BallTree

class OpSimSummarizer:
    def __init__(self, db_path, table_name="observations"):
        self.sql_engine = self._get_sql_engine(db_path)
        self.opsimdf = self._get_df_from_sql(self.sql_engine)
        self.tree = BallTree(self.opsimdf [['_dec', '_ra']].values,
                             leaf_size=40,
                             metric='haversine')
        
        self._hp_survey = None

    @staticmethod
    def _get_sql_engine(dbname):
        """read a sql Opsim output file.

        Parameters
        ----------
        dbname : _type_
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
        print(' Reading from database {}'.format(dbname))
        engine = sqla.create_engine(dbname, echo=False)
        return engine
    
    @staticmethod
    def _get_df_from_sql(sql_engine):
        df = pd.read_sql('observations', con=sql_engine)
        df['_ra'] = np.radians(df.fieldRA)
        df['_dec'] = np.radians(df.fieldDec)
        df.set_index('observationId', inplace=True)
        return gpd.GeoDataFrame(data=df, 
                                geometry=gpd.points_from_xy(df['_ra'], df['_dec']))
        
    def compute_hpix_survey(self, nside=256, minVisits=None, maxVisits=None, ang_radius=1.75):
        ipix = np.arange(hp.nside2npix(nside))
        hp_ra, hp_dec = np.radians(hp.pix2ang(nside, ipix, lonlat=True))
        
        self._hp_survey = pd.DataFrame(dict(ipix=ipix, hp_ra=hp_ra, hp_dec=hp_dec))
        self._hp_survey['n_visits'] = self.tree.query_radius(self.hp_survey[['hp_dec', 'hp_ra']].values, 
                                                            r=np.radians(ang_radius), 
                                                            count_only=True)
        
        
        # Apply mask if needed
        visits_mask = np.ones(len(self._hp_survey), dtype='bool')
        if minVisits is not None:
            visits_mask &= self._hp_survey['n_visits'] >= minVisits
        if maxVisits is not None:
            visits_mask &= self._hp_survey['n_visits'] <= maxVisits
            
        self._hp_survey = self._hp_survey[visits_mask]
        self._hp_survey.set_index('ipix', inplace=True)
        
    @property
    def hp_survey(self):
        if self._hp_survey is None:
            raise ValueError("hp_survey not set. Run compute_hpix_survey before.")
        return self._hp_survey
        

        
        