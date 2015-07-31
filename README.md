# PI_Monkey

### (Probabilistic Inference Monkey data)

This repository contains code used to analyze neural recordings of monkeys while they learned an orientation-discrimination task.

Note that the data folder is deliberately empty, and .mat and .fig files are ignored. Only analysis code is kept here, not data nor results (and so the organization of some of the code may seem strange, given that it is designed to analyze a particularly formatted structure array of data)

---

### Files

The `analyze_` files are the top-level analysis scripts. They can serve as templates for loading and preprocessing data.
The `test_` files also contain some plots, though mostly 'sanity checks' to ensure that the other functions work properly.

Loosely speaking, the lower-case functions and scripts are helpers, and those formatted like Function_Name implement the
main analyses and preprocessing.