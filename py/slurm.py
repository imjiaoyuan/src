#!/usr/bin/env python3
import os, glob, sys
from subprocess import call

INPUT_DIR = "/home/xiejiaxiao/busco/fna"
FILE_SUFFIX = ".fna"
OUTPUT_DIR = "/home/xiejiaxiao/busco/fna/result"
COMMAND_TEMPLATE = "busco -i ${INPUT_DIR}/${FILEID}.fna -l eukaryota_odb12 -o ${OUTPUT_DIR}/${FILEID}.fna -m genome "

JOB_DIR = "./busco_jobs"
LOG_DIR = "./busco_logs"
SLURM_CPUS = "2"
SLURM_MEM = "32"
CONDA_ENV = "busco"

def gen_commands():
    with open("commands.txt", "w") as f:
        for f_path in glob.glob(os.path.join(INPUT_DIR, f"*{FILE_SUFFIX}")):
            fileid = os.path.basename(f_path).replace(FILE_SUFFIX, "")
            cmd = COMMAND_TEMPLATE.replace("${INPUT_DIR}", INPUT_DIR) \
                                .replace("${FILEID}", fileid) \
                                .replace("${OUTPUT_DIR}", OUTPUT_DIR)
            f.write(cmd + "\n")

def build_jobs():
    os.makedirs(JOB_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    for f_path in glob.glob(os.path.join(INPUT_DIR, f"*{FILE_SUFFIX}")):
        fileid = os.path.basename(f_path).replace(FILE_SUFFIX, "")
        job_file = os.path.join(JOB_DIR, f"{fileid}.sl")
        log_prefix = os.path.join(LOG_DIR, f"{fileid}")
        with open(job_file, "w") as f:
            f.write("#!/bin/sh\n")
            f.write(f"#SBATCH -n 1 -c {SLURM_CPUS} --mem={SLURM_MEM}g\n")
            f.write(f"#SBATCH -e {log_prefix}.err\n")
            f.write(f"#SBATCH -D {os.getcwd()}\n\n")
            f.write("source ~/miniforge3/etc/profile.d/conda.sh\n")
            f.write(f"conda activate {CONDA_ENV}\n\n")
            cmd = COMMAND_TEMPLATE.replace("${INPUT_DIR}", INPUT_DIR) \
                                .replace("${FILEID}", fileid) \
                                .replace("${OUTPUT_DIR}", OUTPUT_DIR)
            f.write(cmd + "\n")

def submit_jobs():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for job_file in glob.glob(os.path.join(JOB_DIR, "*.sl")):
        call(["sbatch", job_file])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 script.py [gen_commands|build_jobs|submit_jobs]")
    cmd = {"gen_commands": gen_commands, "build_jobs": build_jobs, "submit_jobs": submit_jobs}
    cmd.get(sys.argv[1], lambda: sys.exit("Invalid command"))()