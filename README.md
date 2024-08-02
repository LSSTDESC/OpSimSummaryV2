[![Documentation Status](https://readthedocs.org/projects/opsimsummaryv2/badge/?version=latest)](https://opsimsummaryv2.readthedocs.io/en/latest/?badge=latest)

## Presentation
`OpSimSummaryV2` is a rework of the previous [`OpSimSummary`](https://github.com/LSSTDESC/OpSimSummary) library. This code is developed to interact with the LSST Operations Simulator outputs available [here](http://astro-lsst-01.astro.washington.edu:8080/?runId=1). The OpSim output can be converted into simulation inputs: [SNANA](https://github.com/RickKessler/SNANA) SIMLIB files or [SNSim](https://github.com/bastiencarreres/snsim) observations files.

## Documentation
See the [ReadTheDoc page](https://opsimsummaryv2.readthedocs.io/en/latest/installation.html) (under construction).

## Installation  and Software Requirements

```
git clone https://github.com/LSSTDESC/OpSimSummaryV2.git
cd opsimsummaryv2
pip install .
```
If you want to use host matching feature install using
```
pip install .[hostmatch]
```
