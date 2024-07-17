import pandas as pd
import sqlalchemy as sqa

class OpSimSummarizer:
    def __init__(self, dbname, table_name="observations"):
        self.sql_engine = self._get_sql_engine(dbname)
        self.opsimdf = pd.read_sql(table_name, self.sql_engine)
        

    @staticmethod
    def _get_sql_engine(dbname):
        # Prepend the abs path with sqlite for use with sqlalchemy
        if not dbname.startswith('sqlite'):
            dbname = 'sqlite:///' + dbname
        print(' Reading from database {}'.format(dbname))
        engine = sqa.create_engine(dbname, echo=False)
        return engine
    
    