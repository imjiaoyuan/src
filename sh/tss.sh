#!/bin/sh
#SBATCH -n 4 -c 8 --mem=8g
#SBATCH -e ./tssdis.err
#SBATCH -D /share/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/bwfiles

source ~/miniforge3/etc/profile.d/conda.sh
conda activate atac

BED_FILE="/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/DM_hap1.bed"

cd /share/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/bwfiles
# computeMatrix reference-point \
#     --referencePoint TSS \
#     --missingDataAsZero \
#     -b 3000 -a 3000 \
#     -R "$BED_FILE" \
#     -S CK_ATAC-seq_24h-1.bw CK_ATAC-seq_24h-2.bw TR_ATAC-seq_24h-1.bw TR_ATAC-seq_24h-2.bw \
#     -o "TSS_combined_matrix.gz"

computeMatrix scale-regions \
  -S CK_ATAC-seq_24h-1.bw CK_ATAC-seq_24h-2.bw TR_ATAC-seq_24h-1.bw TR_ATAC-seq_24h-2.bw \
  -R /home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/DM_hap1.bed \
  -a 3000 -b 3000 \
  --regionBodyLength 5000 \
  --skipZeros \
  -o matrix_gene_body.gz


#!/bin/sh
#SBATCH -n 4 -c 8 --mem=8g
#SBATCH -e ./tssdis.err
#SBATCH -D /share/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/bwfiles
source ~/miniforge3/etc/profile.d/conda.sh
conda activate atac
cd /share/home/xiejiaxiao/03.other/00.1/tmp/atac-rawdata/ZY20250418-KJ5332-XM01_ATAC_rawdata_20250527/bwfiles

# plotHeatmap \
#   -m TSS_combined_matrix.gz \
#   -out TSS_heatmap.png \
#   --colorMap viridis \
#   --zMin 0 --zMax 10 \
#   --regionsLabel "Peaks" \
#   --samplesLabel CK1 CK2 TR1 TR2

# plotProfile \
#   -m TSS_combined_matrix.gz \
#   -out TSS_profile.png \
#   --colors blue green red purple \
#   --plotHeight 10 --plotWidth 12

plotHeatmap \
  -m matrix_gene_body.gz \
  -out gene_body_heatmap.png \
  --colorMap RdYlBu \
  --zMin -2 --zMax 2 \
  --samplesLabel CK1 CK2 TR1 TR2

# plotProfile \
#   -m matrix_gene_body.gz \
#   -out gene_body_profile.png \
#   --startLabel "TSS" \
#   --endLabel "TES" \
#   --plotType std