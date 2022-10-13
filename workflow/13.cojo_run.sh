#!usr/bin/env/ bash

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