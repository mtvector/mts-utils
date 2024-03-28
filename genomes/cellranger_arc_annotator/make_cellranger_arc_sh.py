"""
Author: Matthew Schmitz, Allen Institute, 2024 (Represent)

Writes cellranger arc make scripts for all genomes in a table
usage:
Fill in variables and just run python make_cellranger_arc_sh.py
"""

import pandas as pd
import os
import re
import sys
sys.path.append('/home/matthew.schmitz/mts-utils/genomes/cellranger-arc')
import ChromosomeSplitter


# Directory where the Slurm scripts will be saved
sh_scripts_dir = '/home/matthew.schmitz/subscripts/cellranger-arc/'
# Directory where the debug_gtf.py script is located
script_path = '/home/matthew.schmitz/Matthew/code/SpeciesReferenceProcessing/annotations/debug_gtf.py'
out_dir = '/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/cellranger_arc'
cellranger_bin = '/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-arc-2.0.2/bin'
# Ensure the scripts directory exists
os.makedirs(sh_scripts_dir, exist_ok=True)

# Load the CSV file which must have Species, Common Name, TaxID, FASTA Path and GTF Path columns
df = pd.read_csv('./AllNCBIGenomesData.csv',sep=',',index_col=0)

# Iterate over the rows of the DataFrame
for index, row in df.iterrows():
    species = row['Species'].replace(' ', '_')
    common_name = row['Common Name'].replace(' ', '_')
    tax_spec=str(row['TaxID'])+"-"+row['Species']
    genome_assembly = tax_spec+"__"+re.sub('\.','-',re.sub('_genomic.+','',row['Genome Assembly Name']))
    print(genome_assembly)
    fasta_path = row['FASTA Path']
    gtf_path = row['GTF Path']
    split_fasta_path=re.sub('\.fasta|\.fna|\.fa|\.gz','',fasta_path)+'_splitchr.fa'
    split_gtf_path=re.sub('\.gtf|\.gff|\.gz','',gtf_path)+'_splitchr.gtf'
    ChromosomeSplitter(fasta_path,gtf_path,split_fasta_path,split_gtf_path,length_threshold=5e8)
    gtf_debugged_filename = os.path.basename(split_gtf_path).replace('.gtf', '_debugged.gtf')
    gtf_debugged_path = os.path.join(os.path.dirname(gtf_path), gtf_debugged_filename)
    config_file = f"{out_dir}/refseq_{species}.config"

    # Create the Slurm script content
    script_content = f"""#!/bin/bash
#SBATCH --job-name=mkref_{species}
#SBATCH --output=/home/matthew.schmitz/log/mkref_{species}.out
#SBATCH --error=/home/matthew.schmitz/log/mkref_{species}.err
#SBATCH --time=48:00:00
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

# Load Python environment or module if necessary
source ~/.bashrc
conda activate genomes

# Step 1: Process the GTF file with the Python script
python {script_path} "{gtf_path}" "{gtf_debugged_path}"

cd {out_dir}

# Step 2: Create the config file for cellranger-arc
echo "{{
    organism: \\"{common_name}\\"
    genome: [\\"{genome_assembly}\\"]
    input_fasta: [\\"{split_fasta_path}\\"]
    input_gtf: [\\"{gtf_debugged_path}\\"]
}}" > {config_file}


# Step 3: Run cellranger-arc mkref
export PATH=$PATH:{cellranger_bin}
cellranger-arc mkref --config={config_file}
cp {config_file} {genome_assembly}

{cellranger_bin}/gtf_to_gene_index {out_dir}/{genome_assembly} test.json 
"""

    # Save the script to a file
    script_filename = os.path.join(sh_scripts_dir, f"slurm_{species}.sh")
    with open(script_filename, 'w') as script_file:
        script_file.write(script_content)

    print(f"Slurm script for {species} saved to {script_filename}")

print("All Slurm scripts have been generated.")
print('submit with: for script in *.sh; do sbatch "$script"; done ')