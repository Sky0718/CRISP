# DawnRank-derived routines

DawnRank was developed by Jack Hou. Its source package declares the GNU General
Public License. The TARGET-SL distribution containing the corresponding routines
includes GNU GPL version 3.

The following files contain DawnRank-derived transition, damping, iterative
ranking or output routines and are distributed under GNU GPL version 3, whose
complete text is provided in [GPL-3.0.txt](GPL-3.0.txt):

- `src/crisp/core/r/dawn.R`
- `src/crisp/precision/solver.R`
- `src/crisp/analyses/Allocation_and_Graph_Controls/Graph_Controls/runtime/native/original_adapter.R`
- `src/crisp/analyses/Allocation_and_Graph_Controls/Graph_Controls/runtime/native/functions.R`

The CRISP adaptations add sparse operators, explicit input contracts, convergence
diagnostics, error bounds and allocation interfaces. The Apache licence for
original project code does not replace the licence of these derived routines.

Sources:

- [DawnRank source, commit e34fefd2e0d52ff6fa8b3384cbf6d24eac4766fd](https://github.com/MartinFXP/DawnRank/tree/e34fefd2e0d52ff6fa8b3384cbf6d24eac4766fd)
- [TARGET-SL source, commit ec057f299b3b7017e743660b3a96bcaab9e7fb4f](https://github.com/RhysGillman/TARGET-SL/tree/ec057f299b3b7017e743660b3a96bcaab9e7fb4f)
