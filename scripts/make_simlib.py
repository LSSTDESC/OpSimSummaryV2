import opsimsummaryv2 as op
import argparse

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
                    default=50000)

parser.add_argument('--host_file', '-hf',
                    help='absolute path to a host file.',
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
                    default=None)

parser.add_argument('--max_MJD',
                    help='Maximum date to query',
                    default=None)

parser.add_argument('--output_dir',
                    help='Output dir or file for the SIMLIB',
                    default='./')

parser.add_argument('--random_seed', '-rs',
                    help='Random seed for survey sampling',
                    default=None)

parser.add_argument('--limit_numpy_threads', '-np_threads',
                    help='Limit the number of threads numpy could  use.',
                    default=4)

parser.add_argument('--n_cpu',
                    help='Number of cpu to use for matching survey and hosts.',
                    default=10)


args = parser.parse_args()

limit_numpy(args.limit_numpy_threads)

# format MJD range
MJDrange = None
if args.min_MJD is not None:
    MJDrange = [args.min_MJD, 1e15]
if args.max_MJD is not None:
    if args.min_MJD is not None:
        MJDrange[1] = args.max_MJD
    else:
        MJDrange = [0, args.max_MJD]

OpSimSurv = op.OpSimSurvey(args.db_file, 
                           MJDrange=MJDrange,
                           host_file=args.host_file,
                           host_config={'col_ra': args.hf_RA_col, 
                                        'col_dec': args.hf_DEC_col,
                                        'ra_dec_unit':  args.hf_radec_unit})

# Compute healpy rep
OpSimSurv.compute_hp_rep(nside=256, minVisits=500, maxVisits=10000)

# Sample survey
OpSimSurv.sample_survey(args.Nfields, random_seed=args.random_seed, nworkers=args.n_cpu)

# Write the SIMLIB
sim = op.sim_io.SNANA_Simlib(OpSimSurv, out_path=args.output_dir)
sim.write_SIMLIB(write_batch_size=100)


