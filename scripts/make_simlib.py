import opsimsummaryv2 as op
import argparse
import sys
import os
from pathlib import Path

def limit_numpy(threads=4):
    os.environ["OMP_NUM_THREADS"] = "4"  # export OMP_NUM_THREADS=4
    os.environ["OPENBLAS_NUM_THREADS"] = "4"  # export OPENBLAS_NUM_THREADS=4
    os.environ["MKL_NUM_THREADS"] = "4"  # export MKL_NUM_THREADS=6
    os.environ["VECLIB_MAXIMUM_THREADS"] = "4"  # export VECLIB_MAXIMUM_THREADS=4
    os.environ["NUMEXPR_NUM_THREADS"] = "4"  # export NUMEXPR_NUM_THREADS=6


parser = argparse.ArgumentParser(
    prog="OpSimSummaryV2 script",
    description="Execute OpSimSummaryV2 to create a SIMLIB from a OpSim output database.",
)


# Args DB file
parser.add_argument("db_file", help="Absolute path to the OpSim database. Only the version could be specified if used with download option '-d' flag.")
parser.add_argument(
    "--download",
    "-d",
    help="Dowload the baseline db file if it did not exist",
    action="store_true",
)


# Args survey sampling
parser.add_argument(
    "--Nfields", "-Nf", help="Number of fields to sample.", default=50_000, type=int
)

parser.add_argument(
    "--hp_nside", "-hpns", help="Nside resolution of healpy pixel.", default=256, type=int
)


parser.add_argument(
    "--min_visits", help="Minimum observation visits", default=500, type=int
)

parser.add_argument(
    "--max_visits", help="Maximum observation visits", default=1e15, type=int
)

parser.add_argument("--wfd_only", help="Only write WFD", action="store_true")

parser.add_argument("--ddf_only", help="Only write WFD", action="store_true")

parser.add_argument("--ddf_tags", help="Write DDF tags", action="store_true")


parser.add_argument(
    "--ddf_nobs_thresh",
    help="Number of observations threshold between WFD and DDF fields",
    default=1100,
    type=int,
)

parser.add_argument(
    "--random_seed",
    "-rs",
    help="Random seed for survey sampling",
    default=None,
    type=int,
)

parser.add_argument(
    "--simlib_coadd", help="Run SNANA simlib_coadd.exe. Consider using it if you are simulating DDF.", action="store_true"
)

# Args hosts
parser.add_argument(
    "--host_file", "-hf", help="absolute path to a host file.", default=None
)

parser.add_argument(
    "--snana_wgtmap", help="SNANA WGTMAP to apply to host.", default=None, type=str
)

parser.add_argument("--author", "-auth", help="Author of the file.", default=None)

parser.add_argument(
    "--hf_RA_col", "-hfra", help="RA column keys in host file", default="RA_GAL"
)

parser.add_argument(
    "--hf_DEC_col", "-hfdec", help="DEC column keys in host file", default="DEC_GAL"
)

parser.add_argument(
    "--hf_radec_unit", help="DEC column keys in host file", default="degrees"
)


# Args MJD range
parser.add_argument("--min_MJD", help="Minimum date to query", default=0, type=float)

parser.add_argument("--max_MJD", help="Maximum date to query", default=1e15, type=float)

parser.add_argument("--output_dir", help="Output dir for the SIMLIB. If not specified, the code will search for $", default=None)


# Args computation
parser.add_argument(
    "--limit_numpy_threads",
    "-np_threads",
    help="Limit the number of threads numpy could  use.",
    default=4,
    type=int,
)

parser.add_argument(
    "--n_cpu",
    help="Number of cpu to use for matching survey and hosts.",
    default=10,
    type=int,
)


###############
## MAIN CODE ##
###############

args = parser.parse_args()

# Limit numpy cpu usage
limit_numpy(args.limit_numpy_threads)

# COPY THE COMMAND USED TO STORE IT
command = sys.argv.copy()
command[0] = Path(command[0]).name

NOTES = {"COMMAND": str(command)}
SNANA_LSST_ROOT = os.getenv("SNANA_LSST_ROOT")

# CHECK & SET OUTPUT_DIR
if args.output_dir is None:
    if SNANA_LSST_ROOT is not None:
        output_dir = Path(SNANA_LSST_ROOT) / "simlibs"
    else:
        raise ValueError("Set output_dir or SNANA_LSST_ROOT env var.")
else:
    output_dir = Path(args.output_dir)

# CHECK & SET DB FILE
db_file = Path(args.db_file)
if not db_file.exists():
    if SNANA_LSST_ROOT is not None and (Path(SNANA_LSST_ROOT) / 'simlibs' / db_file.name).exists():
        db_file = Path(SNANA_LSST_ROOT) / 'simlibs' / db_file.name
    elif args.download:
        import re

        version = re.findall("v([0-9].[0-9])", db_file.name)[0]
        print(f'Downloading OpSim output v{version}')
        op.utils.download_rubinlsst_baseline_dbfile(version, output_dir=str(output_dir))
        db_file = output_dir / f"baseline_v{version}_10yrs.db"
    else:
        raise ValueError(
            "db file does not exist please verify path or use -d flag to download it"
        )

# MJD RANGE
MJDrange = [args.min_MJD, args.max_MJD]

# SET MIN AND MAX VISITS
minVisits = args.min_visits
maxVisits = args.max_visits
NOTES["WARNING_NOTICE"] = (
    f"DDF approximately determined by nobs>={args.ddf_nobs_thresh}"
)

file_suffix = "_all"
if args.wfd_only and args.ddf_only:
    raise ValueError("wfd_only and ddf_only options can not be set at the same time")
elif args.wfd_only:
    print(f"Replace max_visits by WFD/DDF threshold: nobs < {args.ddf_nobs_thresh}")
    maxVisits = args.ddf_nobs_thresh
    file_suffix = "_WFD"
elif args.ddf_only:
    print(
        f"Replace min_visits by WFD/DDF threshold: nobs >= {args.ddf_nobs_thresh}"
    )
    minVisits = args.ddf_nobs_thresh
    file_suffix = "_DDF"


host_config = {
    "col_ra": args.hf_RA_col,
    "col_dec": args.hf_DEC_col,
    "ra_dec_unit": args.hf_radec_unit,
}

if args.snana_wgtmap is not None:
    host_config["wgt_map"] = args.snana_wgtmap

OpSimSurv = op.OpSimSurvey(
    str(db_file), MJDrange=MJDrange, host_file=args.host_file, host_config=host_config
)

# Compute healpy rep
OpSimSurv.compute_hp_rep(
    nside=args.hp_nside, 
    minVisits=minVisits, 
    maxVisits=maxVisits,
    ddf_nobs_thresh=args.ddf_nobs_thresh,
    add_ddf_tag=args.ddf_tags
    )

# Sample survey
OpSimSurv.sample_survey(args.Nfields, random_seed=args.random_seed, nworkers=args.n_cpu)

# Write the SIMLIB
sim = op.sim_io.SNANA_Simlib(
    OpSimSurv,
    out_path=str(output_dir),
    author_name=args.author,
    NOTES=NOTES,
    file_suffix=file_suffix,
)
sim.write_SIMLIB(write_batch_size=100)

if args.simlib_coadd:
    import subprocess
    
    print('\n\n#################\n\n')
    print('Executing simlib_coadd.exe')
    
    SNANA_DIR = os.getenv("SNANA_DIR")
    if SNANA_DIR is None:
        raise ValueError("$SNANA_DIR need to be defined to use coadd")
    p = subprocess.run([f"{SNANA_DIR}/bin/simlib_coadd.exe", f"{sim.out_path}", "SORT_BAND"])
    print('Return code:', p.returncode)