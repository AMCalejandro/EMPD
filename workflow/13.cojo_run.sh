#!usr/bin/env/ bash

# Add cojo conditional and stepwise model selection in GJA5 locus
gcta64 cojo--bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly
--chr 1 \
--maf 0.01 \
--cojo-file rsloped_v2_forcojo.ma \
--cojo-slct \
--out COJO_CHR1

gcta64 --bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly \
--chr 1 \
--maf 0.01 \
--cojo-file rsloped_v2_forcojo.ma \
--cojo-cond cond.snplist \
--out COJO_CHR1


# Add cojo conditional and stepwise model selection in MAD1L1 locus
gcta64 cojo--bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly
--chr 7 \
--maf 0.01 \
--cojo-file intercept_chr17_forcojo.ma \
--cojo-slct \
--out COJO_CHR7

gcta64 --bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly \
--chr 7 \
--maf 0.01 \
--cojo-file intercept_chr17_forcojo.ma \
--cojo-cond cond.snplist \
--out COJO_CHR7


# Add cojo conditional and stepwise model selection in LINC00511 locus
gcta64 cojo--bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly
--chr 17 \
--maf 0.01 \
--cojo-file intercept_chr17_forcojo.ma \
--cojo-slct \
--out COJO_CHR17

gcta64 --bfile /mnt/rreal/RDS/DATA/AMPPD/AMPPD_w_rsids_rsidsOnly \
--chr 17 \
--maf 0.01 \
--cojo-file intercept_chr17_forcojo.ma \
--cojo-cond cond.snplist \
--out COJO_CHR17
