#!/usr/bin/env python3
import os
import subprocess
import shutil

DOWNLOAD_PROTEIN = False
DOWNLOAD_GFF = False

GENOME_LIST = [
    "GCF_901000725.2",  
]

def main():
    output_dir = "downloaded_genomes"
    os.makedirs(output_dir, exist_ok=True)

    if not shutil.which("datasets"):
        print("Error: 'datasets' command not found.")
        return

    include_opts = ["genome"]
    if DOWNLOAD_PROTEIN: include_opts.append("protein")
    if DOWNLOAD_GFF: include_opts.append("gff3")
    include_str = ",".join(include_opts)

    print("Starting downloads...")

    for accession_id in GENOME_LIST:
        print(f"Downloading {accession_id}...")
        zip_file = os.path.join(output_dir, f"{accession_id}.zip")
        cmd = [
            "datasets", "download", "genome", "accession", accession_id,
            "--include", include_str, "--filename", zip_file
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            if shutil.which("unzip"):
                unzip_dir = os.path.join(output_dir, accession_id)
                subprocess.run(["unzip", "-oq", zip_file, "-d", unzip_dir], check=True)
                os.remove(zip_file)
            print(f"Success: {accession_id}")
        except subprocess.CalledProcessError as e:
            print(f"Failed: {accession_id} | Error: {e.stderr.decode().strip()}")

    print("All tasks finished.")

if __name__ == "__main__":
    main()