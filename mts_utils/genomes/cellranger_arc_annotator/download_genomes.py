"""
Author: Matthew Schmitz, Allen Institute, 2024 (Represent)

Downloads NCBI reference genomes and annotations for a list of species
usage:
Fill in variables and just run python download_genomes.py Could also replace species list with a dataframe, or species dict as a dictionary of {taxid:scientific_name}
"""
import sys
sys.path.append('/home/matthew.schmitz/mts-utils/genomes/cellranger-arc')
import genome_obtainment
from genome_obtainment import *


genome_directory='/home/matthew.schmitz/Matthew/genome/ncbi2'
genome_csv='./AllNCBIGenomesData.csv'

species_list = [
    "delphinus delphis", "tursiops truncatus", "monodon monoceros"
]

# Fetch TaxIDs
species_dict = {get_taxid(species):species  for species in species_list}


species_list = species_dict.keys()
for taxid in species_list:
    download_genome_data(str(taxid),genome_directory)


ncbi_data = list(process_directory(genome_directory,gtf=True))

ncbi_data = [x for x in ncbi_data if x is not None]

# Create DataFrames
ncbi_df = pd.DataFrame(ncbi_data)

# Combine DataFrames
combined_df = pd.concat([ ncbi_df, other_df])
print(combined_df)
combined_df.to_csv(genome_csv)