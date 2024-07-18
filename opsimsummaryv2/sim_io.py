import os
from pathlib import Path
import datetime
import opsimsummaryv2 as oss
import numpy as np

class SNANA_Simlib():
    _pixelSize = 0.2

    def __init__(self, OpSimSurvey, output_path, author_name=None):
        self.OpSimSurvey = OpSimSurvey
        self.author_name = author_name
        self.output_path = Path(output_path)
    
    def get_survey_doc(self):
        dt = datetime.datetime.now()
        ts = dt.isoformat()
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
        doc += '    VERSIONS:\n'
        doc += f'    - DATE : {dt.strftime(format="%y-%m-%d")}\n'
        doc += f'    AUTHORS : {self.author_name}, OpSimSummaryV2 version {oss.__version__}\n'
        doc += 'DOCUMENTATION_END:\n'
        return doc
    
    def get_survey_header(self, saturation_flag=1024, comments='\n'):
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
    
    def LIBheader(self, LIBID, ra, dec, opsimtable, mwebv=0.0):
        """
        Parameters
        ----------
        fieldID : int
            integer for the unique field ID
        ra : float, degrees
            ra of the field location
        dec : float, degrees
            dec of the field location
        opsimtable : `np.array` of `pd.DataFrame`
            sequence of OpSim observations in above format to find number of
            observations.
        mwebv : float, defaults to 0.0
            milky way E(B-v) value. This is usually recomputed in SNANA
            depending on flags, and hence can be left as 0.0
        fieldtype : string, defaults to None
            string used to construct `Field: fieldtype` line, if None this
            line is left out.
        """
        nobs = len(opsimtable)
        # String formatting
        s = '# --------------------------------------------' +'\n' 
        s += 'LIBID: {0:10d}'.format(LIBID) +'\n'
        tmp = 'RA: {0:+10.6f} DEC: {1:+10.6f}   NOBS: {2:10d} MWEBV: {3:5.2f}'
        tmp += ' PIXSIZE: {4:5.3f}'
        s += tmp.format(ra, dec, nobs, mwebv, self._pixelSize) + '\n'
        s += '#                           CCD  CCD         PSF1 PSF2 PSF2/1' +'\n'
        s += '#     MJD      ID*NEXPOSE  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n'
        return s
        
    def LIBdata(self, obs):
        obs = self.formatObs(obs)
        for ObsID, row in opsimtable.iterrows():
            lst = ['S:',
                   "{0:5.4f}".format(row.observationStartMJD),
                   "{0:10d}*2".format(ObsID),
                   data['filter'],
                   "{0:5.2f}".format(1.),                  # CCD Gain
                   "{0:5.2f}".format(0.25),                # CCD Noise
                   "{0:6.2f}".format(row.SKYSIG),          # SKYSIG
                   "{0:4.2f}".format(row.PSF),             # PSF1
                   "{0:4.2f}".format(0.),                  # PSF2
                   "{0:4.3f}".format(0.),                  # PSFRatio
                   "{0:6.2f}".format(row.ZPTAVG),          # ZPTAVG
                   "{0:6.3f}".format(0.005),               # ZPTNoise 
                   "{0:+7.3f}".format(-99.)]               # MAG
            lib = ' '.join(lst)
            lib += '\n'
        return lib
    
    def formatObs(self, obs):
        opsim_seeing = obs['seeingFwhmEff']
        
        # magsky is in units of mag/arcsec^2
        # opsim_maglim is in units of mag
        opsim_maglim = obs['fiveSigmaDepth']
        opsim_magsky = obs['skyBrightness']
                
        # COMPUTE SIMLIB PSF VALUE
        obs['PSF'] = opsim_seeing / 2.35 / self._pixelSize
        
        # COMPUTE SIMLIB ZPT
        opsim_snr = 5.
        area = (1.51 * opsim_seeing)**2
        arg = area * opsim_snr**2
        
        # Background dominated limit assuming counts with system transmission only
        # is approximately equal to counts with total transmission
        zpt_approx = 2.0 * opsim_maglim - opsim_magsky + 2.5 * np.log10(arg)
        
        # Additional term to account for photons from the source, again assuming
        # that counts with system transmission approximately equal counts with total
        # transmission.
        tmp = 10.0**(-0.4 * (opsim_magsky - opsim_maglim))
        zpt_cor = 2.5 * np.log10(1.0 + 1.0 / (area * tmp))
        
        obs['ZPTAVG'] =  zpt_approx + zpt_cor
        
        # COMPUTE SKYSIG 
        # npix_asec = 1. / self._pixelSize**2.
        obs['SKYSIG'] = self._pixelSize * np.sqrt(10.0 **(-0.4 * (opsim_magsky - obs['ZPTAVG'])))
        return obs
    
    def writeSimlib(self):
        with open(self.output_path, 'w') as simlib_file:
            simlib_file.write(self.get_survey_doc())
            simlib_file.write(self.get_survey_header())
            for (i, field), obs in zip(self.OpSimSurv.survey.iterrows(), self.OpSimSurv.get_survey_obs()):
                LIBID = i
                RA = np.degrees(field.hp_ra)
                DEC = np.degrees(field.hp_dec)
                
                simlib_file.write(self.LIBheader(LIBID, RA, DEC, obs))




        