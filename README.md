`OpSimSummaryV2` is a rework of the previous [`OpSimSummary` library](https://github.com/LSSTDESC/OpSimSummary).

`OpSimSummary` is a codebase developed to interact with the LSST Operations Simulator outputs. Currently they are used for catalog Time Domain Simulations. 
This includes a library that can be called by simulation codes to obtain the set of LSST pointings observing a particular point, as well as a script which uses
the library and precomputes such pointings and store them in an observation library. This storage is in a format specific to [`SNANA`](http://snana.uchicago.edu/)

## Documentation
See the [ReadTheDoc page](https://opsimsummaryv2.readthedocs.io/en/latest/installation.html) (under construction).

## Installation  and Software Requirements

```
git clone https://github.com/bastiencarreres/OpSimSummaryV2.git
cd opsimsummaryv2
pip install .
```
