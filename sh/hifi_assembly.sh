samtools fastq -T '*' I153-1个水稻叶片HiFi/25ZF10738/r84164_250530_001_1_C01/ZHN-pa.bc2026.bc2026.HiFi.bam > /I153-1个水稻叶片HiFi/25ZF10738/r84164_250530_001_1_C01/hifi_reads.fq

################################################################################################
#!/bin/sh
#SBATCH -n 1 -c 48 --mem=240g
#SBATCH -e ./hifiasm.err
#SBATCH -o ./hifiasm.out
#SBATCH -D /home/xiejiaxiao/05.hifihicont

source ~/miniforge3/etc/profile.d/conda.sh
conda activate hifi

hifiasm -o rice.asm -t 48 \
  --h1 I160-1个水稻叶片HiC/01.rawFq/00.mergeRawFq/ZHN-hic_I119A19/ZHN-hic_I119A19_raw_1.fq.gz \
  --h2 I160-1个水稻叶片HiC/01.rawFq/00.mergeRawFq/ZHN-hic_I119A19/ZHN-hic_I119A19_raw_2.fq.gz \
  I153-1个水稻叶片HiFi/25ZF10738/r84164_250530_001_1_C01/hifi_reads.fq

awk '/^S/{print ">"$2"\n"$3}' rice.asm.hic.p_ctg.gfa > rice.hic.p_ctg.fasta

################################################################################################

/home/xiejiaxiao/05.hifihicont
mkdir -p juicer_work/fastq

ln -sf rice.asm.hic/rice.hic.p_ctg.fasta ./reference.fasta
ln -sf I160-1个水稻叶片HiC/01.rawFq/00.mergeRawFq/ZHN-hic_I119A19/ZHN-hic_I119A19_raw_1.fq.gz ./juicer_work/fastq/
ln -sf I160-1个水稻叶片HiC/01.rawFq/00.mergeRawFq/ZHN-hic_I119A19/ZHN-hic_I119A19_raw_2.fq.gz ./juicer_work/fastq/

################################################################################################

# 生成chrom.sizes文件
awk '/^>/ {if (seqname) print seqname "\t" length(seq); seqname=substr($1,2); seq=""} !/^>/ {seq = seq $0} END {print seqname "\t" length(seq)}' ./reference.fasta > ./rice.chrom.sizes

################################################################################################

#!/bin/sh
#SBATCH -n 1 -c 48 --mem=240g
#SBATCH -e ./juicer.err
#SBATCH -o ./juicer.out
#SBATCH -D /home/xiejiaxiao/05.hifihicont

source ~/miniforge3/etc/profile.d/conda.sh
conda activate hifi
bash /home/xiejiaxiao/juicer/CPU/juicer.sh -g rice_genome -s MobI -z ./reference.fasta -d ./juicer_work -D /home/xiejiaxiao/juicer -t 48 -p ./rice.chrom.sizes

################################################################################################


run-asm-pipeline.sh -r 0 ./reference.fasta ./juicer_work/aligned/merged_nodups.txt

mv reference.final.fasta rice.final_assembly.fasta

cd ..


### --- 步骤 3: 重复序列注释 (请替换为你的封装命令) ---
cd 03_repeat_annotation

ln -sf 02_juicer_3ddna_scaffolding/rice.final_assembly.fasta ./genome.fasta

# 示例命令，请替换:
trf genome.fasta 2 7 7 80 10 50 2000 -d -h
BuildDatabase -name rice_db -engine ncbi genome.fasta
RepeatModeler -database rice_db -engine ncbi -pa 48
RepeatMasker -pa 48 -lib rice_db-families.fa -gff -nolow genome.fasta

cd ..


### --- 步骤 4: 基因注释 (请替换为你的封装命令) ---
cd 04_gene_annotation

# 假设上一步产出了 genome.fasta.masked
ln -sf 03_repeat_annotation/genome.fasta.masked ./genome.masked.fasta

# 示例命令，请替换为你自己的转录组、同源、从头预测及整合流程:
# 例如: my_gene_annotation_pipeline.sh --input genome.masked.fasta --rna /path/to/rna ...

cd ..