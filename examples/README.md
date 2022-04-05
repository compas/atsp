# ATSP2K README

For an atsp2K calculation, the odd and even states should be divided into
groups of interacting terms (in the Breit-Pauli approximation).  It
is convenient to perform calculations for each group (all Z)
in a separate directory.  For example, in Si-like, the terms of
3s(2)3p(2) 1D, 1S, 3P all interact for some J and are the lowest group of even
parity, so they could be defined as e1.

It has also been found convenient to generate expansions in a directory
files_c prior to actual calculations to make sure that expansions, using
the same rules for all groups, are not too large.  In general then, a
calculation consist of the directory structure:

```
- files_c      contains the .c expansions
- e1, e2, ..   directories for calculations of even groups e1, e2, etc
- o1, o2, ..   directories for calculations of even groups e1, e2, etc
- tr           directory containing the LS transitions
- z14, z15, .. directory containing LSJ transitions for Z=14, 15, etc
```

### Naming conventions

```
  Variables
  ---------_
    parity = o or e
    grp    = 1, 2, 3, ..
    LS     = 2P, 4D, ..
    n      = 3, 4, 5, ..
    At     = Na, P, Cl, ..

  files_c:
  --------
    ${LS}${parity}${grp}.$n.c       for mchf
    ${parity}${grp}.$n.c            for LSJ

  	 ${At}.${LS}${parity}${grp}.$n.c for special cases:

  mchf output for each group
  --------------------------
    ${LS}${parity}${grp}.${z}_$n.l
    ${parity}${grp}.${z}_$n.w

  bp output for each group
  ------------------------
    ${parity}${grp}.${z}_$n.j
```

### Example runs
This directory contains two complete scripts for performing calculations
for forbidden transitions between levels of 3s(2)3p(2) in LS and LSJ.
```
- sh_example_ls: script for performing an LS calculation for the forbidden
                  transitions between the levels of Si-like 3s(2)3p(2)

- sh_example_lsj: shows how the calculations can be extended to Breit-Pauli
```

A new example case as of 2022 is provided for F-like systems - which includes also hyperfine calculations.
