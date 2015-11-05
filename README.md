# PI_Monkey

### (Probabilistic Inference Monkey data)

This repository contains code used to analyze neural recordings of monkeys while they learned an orientation-discrimination task.

<strong><span style="color:#990000">Only our analysis code is kept under version control, not data nor results</span></strong> (and so the organization of some of the code may seem strange, given that it is designed to analyze a particularly formatted structure array of data). The data folder is deliberately empty, and .mat and .fig files are ignored. 

---

### Files

The `analyze_` files are the top-level analysis scripts. They can serve as templates for loading and preprocessing data.
The `test_` files also contain some plots, though mostly 'sanity checks' to ensure that the other functions work properly.

Loosely speaking, the lower-case functions and scripts are helpers, and those formatted like Function_Name implement the
main analyses and preprocessing.

The files `+Vis/boundedlines.m` and `+Vis/outlinebounds.m` come from [kakearney's boundedline repository on Github](https://github.com/kakearney/boundedline-pkg)