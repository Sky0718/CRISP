# Third-party components

Original CRISP project code is distributed under the Apache License 2.0. The bundled DawnRank-derived R routines are an exception, covered by GNU GPL v3; see [DawnRank attribution](../licenses/DawnRank.md) for the affected files and source details. The code incorporates rank-biased overlap routines with the attribution in [RBO.txt](../licenses/RBO.txt).

Native gene-ranking calculations obtain the following implementations separately:

- [TARGET-SL](https://github.com/RhysGillman/TARGET-SL/tree/ec057f299b3b7017e743660b3a96bcaab9e7fb4f), GPL-3.0. PRODIGY uses `data/PRODIGY_pathways.rds` and the `scripts/PRODIGY` sources. DawnRank uses `scripts/DawnRank/DawnNormalize.R`.
- [DawnRank](https://github.com/MartinFXP/DawnRank/tree/e34fefd2e0d52ff6fa8b3384cbf6d24eac4766fd), GPL, for the corresponding implementation lineage.
- [PCSF](https://github.com/IOR-Bioinformatics/PCSF/tree/4e5f2707fdd096b760e5bb4ca60938a3c9836d78), including its solver and licence, for the PRODIGY calculation.

Keep upstream notices with externally obtained software. Network, molecular, screening and annotation data retain their original provider terms. This code repository does not redistribute those provider datasets.
