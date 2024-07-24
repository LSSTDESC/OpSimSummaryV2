import os
import time
from pathlib import Path
import datetime
import opsimsummaryv2 as oss
import numpy as np
from . import utils as ut

class SNANA_Simlib():
    """
    A class to write SNANA simulation output from OpSimSurvey object.

    Attributes
    ----------
    OpSimSurvey : opsimsummaryv2.summary_opsim.OpSimSurvey
        The OpSimSurvey object
    author_name : str, optional
        The author of the SIMLIB
    out_path : str, optional
        The path of the output file, default './'
    date_time : str
        The current time, to be written in SIMLIB
    dataline: numpy.vectorize
        A vectorized function to write SIMLIB dataline
    """
    def __init__(self, OpSimSurvey, out_path=None, author_name=None, ZPTNoise=0.005, CCDgain=1., CCDnoise=0.25):
        self.OpSimSurvey = OpSimSurvey
        self.author_name = author_name
        self.out_path = self._init_out_path(out_path)
        self.date_time = datetime.datetime.now()
        
        self.ZPTNoise = ZPTNoise
        self.CCDnoise = 0.25
        self.CCDgain = 1.
        
        self.dataline = self._init_dataline()
        
    def _init_out_path(self, out_path):
        """Format output path for SIMLIB and HOSTLIB.

        Parameters
        ----------
        out_path : str
            Output directory or file.

        Returns
        -------
        pathlib.Path
            Path to the output SIMLIB.
        """        
        if out_path is None:
            out_path = './'
            
        out_path = Path(out_path)
        
        if out_path.is_dir():
            out_path /= Path(self.OpSimSurvey.opsimdf.attrs['OpSimFile']).stem  + '.SIMLIB'
        if self.OpSimSurvey.survey_hosts is not None:
            out_path = out_path.with_stem(f"{out_path.stem}_{Path(self.OpSimSurvey.host.attrs['file']).stem}")
        return out_path
    
    def get_SIMLIB_doc(self):
        """Give the DOCUMENTATION string for SIMLIB.

        Returns
        -------
        str
            DOCUMENTATION string
        """        
        minMJD = self.OpSimSurvey.opsimdf.observationStartMJD.min()
        maxMJD = self.OpSimSurvey.opsimdf.observationStartMJD.max()
        OpSimFile = self.OpSimSurvey.opsimdf.attrs['OpSimFile']
        doc = 'DOCUMENTATION:\n'
        doc += f'    PURPOSE: simulate LSST based on mock opsim version {OpSimFile}\n'
        doc += '    INTENT:   Nominal\n'
        doc += '    USAGE_KEY: SIMLIB_FILE\n'
        doc += '    USAGE_CODE: snlc_sim.exe\n'
        doc += '    VALIDATION_SCIENCE: \n'
        doc += '    NOTES: \n'
        doc += '        PARAMS MINMJD: {:.4f}\n'.format(minMJD)
        doc += '        PARAMS MAXMJD: {:.4f}\n'.format(maxMJD)
        doc += '        PARAMS TOTAL_AREA: {:.3f}\n'.format(self.OpSimSurvey.survey.attrs['survey_area_deg'])
        doc += '        PARAMS SOLID_ANGLE: {:.3f}\n'.format(self.OpSimSurvey.survey.attrs['survey_area_rad'])
        if self.OpSimSurvey.survey_hosts is not None:
            doc += '        ASSOCIATED HOSTLIB: {}\n'.format(self.out_path.with_suffix('.HOSTLIB').absolute())
        doc += '    VERSIONS:\n'
        doc += f'    - DATE : {self.date_time.strftime(format="%y-%m-%d")}\n'
        doc += f'    AUTHORS : {self.author_name}, OpSimSummaryV2 version {oss.__version__}\n'
        doc += 'DOCUMENTATION_END:\n'
        return doc
    
    def get_SIMLIB_header(self, saturation_flag=1024, comments='\n'):
        """Give the SIMLIB header string.

        Parameters
        ----------
        saturation_flag : int, optional
            The flag corresponding to saturated obs, by default 1024
        comments : str, optional
            Comments to add to the header, by default '\n'

        Returns
        -------
        str
            The SIMLIB header string.
        """        
        try:
            user = os.getlogin()
        except:
            user = 'NONE'
        try:
            host = os.getenv( 'HOSTNAME' )
        except:
            host = 'NONE'
        # comment: I would like to generalize ugrizY to a sort but am not sure
        # of the logic for other filter names. so ducking for now
        header = '\n\n\n'
        header += 'SURVEY: LSST   FILTERS: ugrizY  TELESCOPE: LSST\n'
        header += 'USER: {0:}     HOST: {1}\n'.format(user, host) 
        header += f'NLIBID: {len(self.OpSimSurvey.survey)}\n'
        header += 'NPE_PIXEL_SATURATE:   100000\n'
        header += f'PHOTFLAG_SATURATE:    {saturation_flag}\n'
        header += comments + '\n'
        header += 'BEGIN LIBGEN\n'
        return header
    
    def LIBheader(self, LIBID, ra, dec, opsimdf, mwebv=0.0, groupID=None):
        """Give the string of the header of a LIB entry.

        Parameters
        ----------
        LIBID : int
            The LIBID of the entry.
        ra : float
            RA [deg] coordinate of the entry 
        dec : float
            Dec [deg] coordinate of the entry 
        opsimdf : pandas.DataFrame
            LIB entry observations
        mwebv: float, optional
            MWEBV of the entry, default = 0.0
        groupID : int, optional
            GROUPID of the entry used to match with HOSTLIB hosts.

        Returns
        -------
        str
            The LIB entry header string.
        """     
        nobs = len(opsimdf)
        # String formatting
        s = '# --------------------------------------------' +'\n' 
        s += 'LIBID: {0:10d}'.format(LIBID) +'\n'
        tmp = 'RA: {0:+10.6f} DEC: {1:+10.6f}   NOBS: {2:10d} MWEBV: {3:5.2f}'
        tmp += ' PIXSIZE: {4:5.3f}'
        s += tmp.format(ra, dec, nobs, mwebv, self.OpSimSurvey.__LSST_pixelSize__) + '\n'
        if groupID is not None:
            s += f"HOSTLIB_GROUPID: {groupID}" + "\n"
        s += '#                               CCD  CCD         PSF1 PSF2 PSF2/1' +'\n'
        s += '#     MJD      ID*NEXPOSE  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n'
        return s
    
    def _init_dataline(self):
        f = lambda expMJD, ObsID, BAND, SKYSIG, PSF, ZPT: ut.dataline(expMJD, ObsID, BAND, 
                                                                      self.CCDgain,  self.CCDnoise, 
                                                                      SKYSIG, PSF, ZPT, self.ZPTNoise)
        
        return np.vectorize(f)

    def LIBdata(self, opsimdf):
        """Give the string of a LIB entry.

        Parameters
        ----------
        opsimdf : pandas.DataFrame
            LIB entry observations

        Returns
        -------
        str
            The str of the LIB entry.
        """        
        opsimdf['BAND'] = opsimdf['BAND'].map(lambda x: x.upper() if x=='y' else x)
        lib = '\n'.join(self.dataline(opsimdf['expMJD'].values, opsimdf['ObsID'].values, 
                                      opsimdf['BAND'].values,opsimdf['SKYSIG'].values,
                                      opsimdf['PSF'].values, opsimdf['ZPT'].values))
        return lib + '\n'
    
    def LIBfooter(self, LIBID):
        """Give the string of a LIB entry footer.

        Parameters
        ----------
        LIBID : int
            The LIBID of the entry

        Returns
        -------
        str
            The string of a LIB entry footer
        """        
        footer = 'END_LIBID: {0:10d}'.format(LIBID)
        footer += '\n\n'
        return footer
        
    def get_SIMLIB_footer(self):
        """Give SIMLIB footer.
        """
        s = 'END_OF_SIMLIB:    {0:10d} ENTRIES'.format(self.OpSimSurvey.survey.attrs['N_fields'])
        return s
    
    def write_SIMLIB(self, write_batch_size=10, buffer_size=8192):  
        """write the SIMLIB (and the HOSTLIB) file(s).
        
        Parameters
        ----------
        write_batch_size : int
            Number of LIBID to write at the same time
        buffer_size : int
            buffering option for open() function
        """              
        tstart = time.time()
        print(f'Writing SIMLIB in {self.out_path}')
        
        with open(self.out_path, 'w', buffering=buffer_size) as simlib_file:
            simlib_file.write(self.get_SIMLIB_doc())
            simlib_file.write(self.get_SIMLIB_header())
            
            bcount = 0
            simlibstr = ''
            for (i, field), obs in zip(self.OpSimSurvey.survey.iterrows(), self.OpSimSurvey.get_survey_obs()):
                LIBID = i
                RA = np.degrees(field.hp_ra)
                DEC = np.degrees(field.hp_dec)
                if self.OpSimSurvey.survey_hosts is not None:
                    groupID = LIBID
                else:
                    groupID = None
                    
                simlibstr += self.LIBheader(LIBID, RA, DEC, obs, groupID=groupID)
                simlibstr += self.LIBdata(obs)
                simlibstr += self.LIBfooter(LIBID)
                
                if not bcount % write_batch_size:
                    simlib_file.write(simlibstr)
                    simlibstr = ''
                    
                bcount +=1
                
            simlib_file.write(simlibstr)
            simlib_file.write(self.get_SIMLIB_footer())
    
        print(f'SIMLIB wrote in {time.time() - tstart:.2f} sec.\n')

        if self.OpSimSurvey.survey_hosts is not None:
            tstart = time.time()
            print(f"Writting {len(self.OpSimSurvey.survey_hosts)} hosts in {self.out_path.with_suffix('.HOSTLIB')}")
            self.write_hostlib(self.OpSimSurvey.survey_hosts.drop(columns=['ra', 'dec']), buffer_size=buffer_size)
            print(f"HOSTLIB file wrote in {time.time() - tstart:.2f} sec.")

    def get_HOSTLIB_doc(self):
        """Give docstring for HOSTLIB file.

        Returns
        -------
        str
            Docstring of the HOSTLIB file
        """        
        doc = ("DOCUMENTATION:\n"
                "PURPOSE: HOSTLIB for LSST based on mock opsim\n"
                "    VERSIONS:\n"
                f"    - DATE : {self.date_time}\n"
                f"    - ASSOCIATED SIMLIB : {self.out_path.absolute()}\n"
                f"    AUTHORS : {self.author_name}\n"
                "DOCUMENATION_END:\n"
                "# ========================\n")
        return  doc

    def get_HOSTLIB_header(self, hostdf):
        """Give HOSTLIB header.

        Parameters
        ----------
        hostdf : pandas.DataFrame
            Hosts dataframe

        Returns
        -------
        str
            Header of HSOTLIB file
        """        
        header = (f"# Z_MIN={hostdf.ZTRUE_CMB.min()} Z_MAX={hostdf.ZTRUE_CMB.max()}\n\n"
                   "VPECERR: 0\n\n")
        return header
    
    def write_HOSTLIB(self, hostdf, buffer_size=8192):
        """Write the HOSTLIB file. Called in write_SIMLIB.
        
        Parameters
        ----------
        hostdf : pandas.DataFrame
            Hosts dataframe
        buffer_size : int
            buffering option for open() function
        """  
        import csv          

        VARNAMES = "VARNAMES:"
        for k in hostdf.columns:
            VARNAMES += f" {k}"

        with open(self.out_path.with_suffix('.HOSTLIB'), 'w', buffering=buffer_size) as hostf:
            hostf.write(self.get_HOSTLIB_doc())
            hostf.write(self.get_HOSTLIB_header(hostdf))
            columns = hostdf.columns.values
            columns = np.insert(columns, 0, 'VARNAMES: ')
            hostdf['VARNAMES: '] = 'GAL: '
            hostdf[columns].to_csv(hostf, sep=' ', index=False, quoting=csv.QUOTE_NONE, escapechar=' ')
            