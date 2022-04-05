Example case provided by Messaoud Nemouchi and Fatima Z. Boualili (Spring, 2022)
Hfs and E1, E2 and M1 calculations for F-like Ti, V, Cr See A&A A25 (2014) 563 for details.

States: ground state (2s(2)2p(5) 2^P^o) = gs
        excited state (2s(1)2p(6) 2^S^e = es

The scripts in this example calculate hfs constants for the two states,
E1 transition between the  two states and E2 and M1 transitions in the ground state.

The script sh_mchf generate (gen-state) and calculate mchf wave functions (wf-mchf-state).
The folders are gs/MCHF for ground state and es/MCHF for excited state.

Filenames are written as follows:

    2Po-n-Z.* and 2Se-n-Z.*, where n is the active space and Z the isotope.


The script sh_bp generate (gen-bp-state) and calculate BP wave functions (wf-bp-state).
The folders are gs/BP for ground state and es/BP for excited state.

Filenames are written as follows:

    2Po-n-Z-bp.* and 2Se-n-Z-bp.*

The script hfs-cons performs hfs calculations with mchf and Breit-Pauli
wave functions and put the .h files in the folder /hyp.

The scripts for E1, E2 and M1 calculations are in the folder /trans.
