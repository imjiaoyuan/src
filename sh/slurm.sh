#!/bin/sh
#SBATCH -n 1 -c 4 --mem=8g
#SBATCH -e ./samtools.err
#SBATCH -D /home/xiejiaxiao/05.hifihicont

source ~/miniforge3/etc/profile.d/conda.sh
conda activate hifi

samtools fastq -T '*' I153-1个水稻叶片HiFi/25ZF10738/r84164_250530_001_1_C01/ZHN-pa.bc2026.bc2026.HiFi.bam > hifi_reads.fq