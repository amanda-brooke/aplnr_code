## Cell Ranger scripts for demultiplexing the fastq files using the ocm barcode
OCM barcodes to conditions mapping is noted in the scripts `run_cellranger_chicago.sh` and `run_cellranger_michigan.sh`. Reference genome can be downloaded from the Cell Ranger official site <https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz>.

### Example source data folder
- 01.RawData/
    - Chicago/
        - Chicago-SI_TT_A3_237W2VLT3_S2_L001_I1_001.fastq.gz
        - Chicago-SI_TT_A3_237W2VLT3_S2_L001_I2_001.fastq.gz
        - Chicago-SI_TT_A3_237W2VLT3_S2_L001_R1_001.fastq.gz
        - Chicago-SI_TT_A3_237W2VLT3_S2_L001_R2_001.fastq.gz
    - Michigan/
        - Michigan-SI_TT_A2_237W2VLT3_S1_L001_I1_001.fastq.gz
        - Michigan-SI_TT_A2_237W2VLT3_S1_L001_I2_001.fastq.gz
        - Michigan-SI_TT_A2_237W2VLT3_S1_L001_R1_001.fastq.gz
        - Michigan-SI_TT_A2_237W2VLT3_S1_L001_R2_001.fastq.gz
