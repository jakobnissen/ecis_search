# Ecis
Extracellular injection systems
Building a model to predict ecis packability from leading 20 aas

Author: Jakob Nybo Nissen
Date of creation: 2022-10-19

## Directory structure
* `raw`: Raw data, e.g. experimental data, or data from external research groups.
  Should not be modified at all.
* `src`: This directory contains code and scripts used to reproduce the results.
  The file `main.py` or `main.jl` should produce all results using only the data
  in directories `raw` and `choices`.
* `tmp`: Directory for throwaway analyses and intermediate results.
  Anything in this directory should be able to be deleted with no big loss
* `cache`: Also for intermediate results, but for content that is troublesome
  to recreate, e.g. results of long-running simuations or long-running computation
* `choices`: For files that are not raw files, but impossible to recreate automatically,
  because they rely on humans (you!) making judgement calls.
* `results`: For final analytic results. `main.jl/py` should write results to this
  directory, primarily
* `paper`: For results related to submission of any papers, e.g. manuscripts or
  publication-ready figures.

## Reproducibility
Hardware: Intel i7-1185G7
OS: Linux 6.2.16-2 Manjaro

Julia 1.9.1
Python 3.10.12
