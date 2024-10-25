import opsimsummaryv2 as op
import argparse
import sys
from pathlib import Path

def limit_numpy(threads=4):
    import os
    os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
    os.environ["OPENBLAS_NUM_THREADS"] = "4" # export OPENBLAS_NUM_THREADS=4 
    os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6
    os.environ["VECLIB_MAXIMUM_THREADS"] = "4" # export VECLIB_MAXIMUM_THREADS=4
    os.environ["NUMEXPR_NUM_THREADS"] = "4" # export NUMEXPR_NUM_THREADS=6

parser = argparse.ArgumentParser(
                    prog='OpSimSummaryV2 script',
                    description='Execute OpSimSummaryV2 to create a SIMLIB from a OpSim output database.',
                    )

parser.add_argument('db_file', 
                    help='absolute path to the opsim database.')

parser.add_argument('--Nfields', '-Nf',
                    help='Number of fields to sample',
                    default=50000, type=int)

parser.add_argument('--host_file', '-hf',
                    help='absolute path to a host file.',
                    default=None)

parser.add_argument('--author', '-auth',
                    help='Author of the file.',
                    default=None)

parser.add_argument('--hf_RA_col', '-hfra', 
                    help='RA column keys in host file',
                    default='RA_GAL')

parser.add_argument('--hf_DEC_col', '-hfdec', 
                    help='DEC column keys in host file',
                    default='DEC_GAL')

parser.add_argument('--hf_radec_unit',
                    help='DEC column keys in host file',
                    default='degrees')

parser.add_argument('--min_MJD',
                    help='Minimum date to query',
                    default=None, type=float)

parser.add_argument('--max_MJD',
                    help='Maximum date to query',
                    default=None, type=float)

parser.add_argument('--min_visits',
                    help='Minimum observation visits',
                    default=500, type=int)

parser.add_argument('--max_visits',
                    help='Maximum observation visits',
                    default=100000, type=int)

parser.add_argument('--output_dir',
                    help='Output dir for the SIMLIB',
                    default='./')

parser.add_argument('--random_seed', '-rs',
                    help='Random seed for survey sampling',
                    default=None, type=int)

parser.add_argument('--limit_numpy_threads', '-np_threads',
                    help='Limit the number of threads numpy could  use.',
                    default=4, type=int)

parser.add_argument('--n_cpu',
                    help='Number of cpu to use for matching survey and hosts.',
                    default=10, type=int)

parser.add_argument('--snana_wgtmap',
                    help='SNANA WGTMAP to apply to host.',
                    default=None, type=str)

parser.add_argument('--wfd_only',
                    help="Only write WFD",
                    action='store_true')

parser.add_argument('--ddf_only',
                    help="Only write WFD",
                    action='store_true')

parser.add_argument('--wfd_ddf_nobs_thresh', 
                    help='Number of observations threshold between WFD and DDF fields',
                    default=1100, type=int)

args = parser.parse_args()

limit_numpy(args.limit_numpy_threads)

command = sys.argv.copy()
command[0] = Path(command[0]).name
NOTES = {'COMMAND': str(command)}


# format MJD range
MJDrange = None
if args.min_MJD is not None:
    MJDrange = [args.min_MJD, 1e15]
if args.max_MJD is not None:
    if args.min_MJD is not None:
        MJDrange[1] = args.max_MJD
    else:
        MJDrange = [0, args.max_MJD]
        
# Set min and max visits
minVisits = args.min_visits
maxVisits = args.max_visits
NOTES['WARNING_NOTICE'] = f'DDF approximately determined by nobs>={args.wfd_ddf_nobs_thresh}'

if args.wfd_only and args.ddf_only:
    raise ValueError('wfd_only and ddf_only options can not be set at the same time')
elif args.wfd_only:
    print(f'Replace max_visits by WFD/DDF threshold: nobs < {args.wfd_ddf_nobs_thresh}')
    maxVisits = args.wfd_ddf_nobs_thresh
    file_suffix = '_wfd'
elif args.ddf_only:
    print(f'Replace min_visits by WFD/DDF threshold: nobs >= {args.wfd_ddf_nobs_thresh}')
    minVisits = args.wfd_ddf_nobs_thresh
    file_suffix = '_ddf'


host_config = {
    'col_ra': args.hf_RA_col, 
    'col_dec': args.hf_DEC_col,
    'ra_dec_unit':  args.hf_radec_unit
    }

if args.snana_wgtmap is not None:
    host_config['wgt_map'] = args.snana_wgtmap

OpSimSurv = op.OpSimSurvey(
    args.db_file, 
    MJDrange=MJDrange,
    host_file=args.host_file,
    host_config=host_config
    )

# Compute healpy rep
OpSimSurv.compute_hp_rep(nside=256, minVisits=minVisits, maxVisits=maxVisits)

# Sample survey
OpSimSurv.sample_survey(args.Nfields, random_seed=args.random_seed, nworkers=args.n_cpu)

# Write the SIMLIB
sim = op.sim_io.SNANA_Simlib(OpSimSurv, out_path=args.output_dir, author_name=args.author, NOTES=NOTES)
sim.write_SIMLIB(write_batch_size=100)