#!/bin/sh
#SBATCH -n 1 -c 8 --mem=32g
#SBATCH -e ./homer.err
#SBATCH -D /share/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527

source ~/miniforge3/etc/profile.d/conda.sh
conda activate atac

findMotifsGenome.pl \
    /home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/peaks/CK_ATAC-seq_24h-1_peaks.narrowPeak \
    /home/xiejiaxiao/03.other/00.1/tmp/DM_T2T-Hap1.fa \
    ./homer_motif_results/CK_ATAC-seq_24h-1_motifs \
    -size 200 -p 8

findMotifsGenome.pl \
    /home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/peaks/CK_ATAC-seq_24h-2_peaks.narrowPeak \
    /home/xiejiaxiao/03.other/00.1/tmp/DM_T2T-Hap1.fa \
    ./homer_motif_results/CK_ATAC-seq_24h-2_motifs \
    -size 200 -p 8

findMotifsGenome.pl \
    /home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/peaks/TR_ATAC-seq_24h-1_peaks.narrowPeak \
    /home/xiejiaxiao/03.other/00.1/tmp/DM_T2T-Hap1.fa \
    ./homer_motif_results/TR_ATAC-seq_24h-1_motifs \
    -size 200 -p 8

findMotifsGenome.pl \
    /home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/peaks/TR_ATAC-seq_24h-2_peaks.narrowPeak \
    /home/xiejiaxiao/03.other/00.1/tmp/DM_T2T-Hap1.fa \
    ./homer_motif_results/TR_ATAC-seq_24h-2_motifs \
    -size 200 -p 8

echo "Motif analysis finished."