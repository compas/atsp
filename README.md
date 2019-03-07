### ATSP - A Multiconfiguration Hartree-Fock (MCHF) Atomic Structure Package 
This is the official repository for the development of ATSP maintained by the [CompAS](https://github.com/compas) team. Fundamental insight into the approach can be obtained from e.g. Froese Fischer et al.'s book "*Computational Atomic Structure: An MCHF Approach*" [[1]](https://www.crcpress.com/Computational-Atomic-Structure-An-MCHF-Approach/Froese-Fischer-Brage-Johnsson/p/book/9780750304665) or the 2016 review article [[2]](https://doi.org/10.1088/0953-4075/49/18/182004) as well as the latest published version from 2007 [[3]](http://dx.doi.org/10.1016/j.cpc.2007.01.006).

### Brief Description
The ATSP suite can be used to determine energy levels and associated wave functions for states of atoms and ions in the MCHF (LS) or Breit–Pauli (LSJ) approximation. Given the wave function, various atomic properties can be computed such as electric (Ek) and magnetic (Mk) multipole radiative transition probabilities (k <= 10) between LS or LSJ states, isotope shift constants, hyperfine parameters, and  factors.

The code is based on dynamic memory allocation, sparse matrix methods, and a recently developed angular library. It is meant for large-scale calculations in a basis of orthogonal orbitals for groups of LS terms of arbitrary parity. For Breit–Pauli calculations, all operators — spin–orbit, spin–other orbit, spin–spin, and orbit–orbit - may be included. For transition probabilities the orbitals of the initial and final state need not be orthogonal. A bi-orthogonal transformation is used for the evaluation of matrix elements in such cases. In addition to transition rates of all types, isotope shifts and hyperfine constants can be computed as well as  factors.

### Published versions
- ATSP2K - *An MCHF atomic-structure package for large-scale calculations* (Froese Fischer et al. CPC, 2007) [[3]](http://dx.doi.org/10.1016/j.cpc.2007.01.006)

- ATSP - *The MCHF atomic-structure package* (Froese Fischer CPC, 2000) [[4]](https://doi.org/10.1016/S0010-4655(00)00009-6)

### References
[1] C. Froese Fischer, T. Brage, and P. Jönsson, 
*Computational Atomic Structure: An MCHF Approach*,
[Institute of Physics (Bristol, 1997)](https://www.crcpress.com/Computational-Atomic-Structure-An-MCHF-Approach/Froese-Fischer-Brage-Johnsson/p/book/9780750304665) (ISBN 9780750304665 - CAT# IP192)

[2] C. Froese Fischer, M. R. Godefroid, T. Brage, P. Jönsson and G. Gaigalas,
*Advanced multiconfiguration methods for complex atoms: I. Energies and wave functions*,
[J. Phys. B: At. Mol. Opt. Phys. **49** 182004 (2016)](https://doi.org/10.1088/0953-4075/49/18/182004)

[3] C. Froese Fischer, G. Tachiev, G. Gaigalas, & M. R. Godefroid,
*An MCHF atomic-structure package for large-scale calculations*, 
[Comput. Phys. Comm., **176(8)** 559 (2007)](http://dx.doi.org/10.1016/j.cpc.2007.01.006)

[4] C. Froese Fischer,
*The MCHF atomic-structure package*, 
[Comput. Phys. Comm. **128** 635 (2000)](https://doi.org/10.1016/S0010-4655(00)00009-6)
