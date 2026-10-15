#!/bin/bash
module load Bioinformatics
module load cellranger/9.0.1

# refdata-gex-GRCh38-2024-A: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
# OCM samples in config.csv (time = hours after BMP4; Ctrl = parental TBXT-mNeonGreen H9):
#   sample1 = OB1 = 42h_Ctrl
#   sample2 = OB2 = 42h_shAPLNR
#   sample3 = OB3 = 72h_Ctrl
#   sample4 = OB4 = 72h_shAPLNR
cellranger multi --id=scrnaseq-meso-chicago-2-15-2026 \
--csv=config.csv
