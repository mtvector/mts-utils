def max_column_match(v,df):
    max_count=0
    max_ind=''
    for c in df.columns:
        matched_count=df[c].isin(v).sum()
        if matched_count>max_count:
            max_count=matched_count
            max_ind=c
    return(max_ind)

def generate_chromsizes(fasta_file, output_file):
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            sequence_name = record.id
            sequence_length = len(record.seq)
            outfile.write(f'{sequence_name}\t{sequence_length}\n')

def get_taxid(species_name):
    """Get the taxid for a species name using NCBI's Entrez Utilities."""
    Entrez.email = "your-email@example.com"  # Provide your email address
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    
    # Return the first TaxID found
    if record["IdList"]:
        return record["IdList"][0]
    else:
        return None


##############################NCBI 2023 Processing################################################

from Bio import Entrez
import os
import requests
from ftplib import FTP
import re
import ftplib
import requests
import zipfile
import os
import shutil
import gzip
Entrez.email = "matthew.schmitz@alleninstitute.org"

def gzip_file(file_path):
    """Compresses a file using gzip."""
    with open(file_path, 'rb') as f_in:
        with gzip.open(file_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)

def download_file(url, filename):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(filename, 'wb') as f:
            f.write(response.content)
    else:
        print(f"Failed to download {url}")

def download_and_unzip(url, target_path, string_to_check):
    # Download the zip file
    response = requests.get(url)
    zip_file_path = os.path.join(target_path, 'downloaded_file.zip')
    with open(zip_file_path, 'wb') as f:
        f.write(response.content)

    # Unzip the file
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(target_path)

    # Move all files from subdirectories to the target directory and gzip if needed
    for root, dirs, files in os.walk(target_path):
        for file in files:
            file_path = os.path.join(root, file)
            if file_path != zip_file_path:  # Avoid moving the zip file itself
                # Check if the file needs to be gzipped
                shutil.move(file_path, os.path.join(target_path, file))
                if (string_to_check) in file and not file.endswith('.gz'):
                    gzip_file(os.path.join(target_path, file))

    # Remove the unzipped directory tree and the zip file
    for root, dirs, files in os.walk(target_path, topdown=False):
        for name in dirs:
            shutil.rmtree(os.path.join(root, name))
    os.remove(zip_file_path)

def find_files_in_directory(directory, file_extension):
    for file in os.listdir(directory):
        if re.search(file_extension, file):
            return file
    return None



def check_ftp_for_gtf(ftp_path):
    ftp_url = "ftp.ncbi.nlm.nih.gov"
    # Extract the relative path from the full FTP path
    relative_ftp_path = '/'.join(ftp_path.split('/')[3:])  # Remove the 'ftp://ftp.ncbi.nlm.nih.gov/' part
    
    try:
        with ftplib.FTP(ftp_url) as ftp:
            ftp.login()
            ftp.cwd(relative_ftp_path)
            files = ftp.nlst()
            for file in files:
                if file.endswith('.gtf.gz'):
                    return True
            return False
    except ftplib.all_errors as e:
        print(f"FTP error: {e}")
        return False

def sanitize_filename(filename):
    # Replace spaces with underscores
    filename = filename.replace(" ", "_")

    # Replace special characters with dashes
    # Define a regex pattern for unwanted characters
    pattern = r'[^\w\-_.]'
    filename = re.sub(pattern, '-', filename)

    return filename

def download_ftp_file(ftp_url, destination_filename):
    # Parse the FTP URL
    ftp_host, ftp_path = ftp_url.replace("ftp://", "").split("/", 1)

    # Connect to the FTP server
    with FTP(ftp_host) as ftp:
        ftp.login()
        
        # Ensure the directory for the file exists
        os.makedirs(os.path.dirname(destination_filename), exist_ok=True)

        # Download the file
        try:
            with open(destination_filename, 'wb') as local_file:
                ftp.retrbinary(f'RETR {ftp_path}', local_file.write)
        except Exception as e:
            print(f"Error downloading {ftp_url}: {e}")


def get_best_assembly_and_download(tax_id,overall_path='',redownload=False):
    handle = Entrez.esummary(db="taxonomy", id=tax_id)
    record = Entrez.read(handle)
    scientific_name = record[0]['ScientificName']

    handle = Entrez.esearch(db="assembly", term=f"{tax_id}[Taxonomy ID]  AND (latest[filter] AND all[filter] NOT anomalous[filter]) AND (\"chromosome\"[Assembly Level] OR \"complete genome\"[Assembly Level] OR )", retmax=50)
    record = Entrez.read(handle)
    best_assembly = None
    best_assembly_gtf = None
    max_contig_n50 = 0
    max_contig_n50_ref = 0
    for assembly_id in record['IdList']:
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(summary_handle,validate=False)
        doc_summary = summary['DocumentSummarySet']['DocumentSummary'][0]
        contig_n50 = int(doc_summary['ContigN50'])
        if contig_n50 > max_contig_n50:
            best_assembly = doc_summary
            max_contig_n50 = contig_n50
        if check_ftp_for_gtf(doc_summary['FtpPath_GenBank']) and (contig_n50 > max_contig_n50_ref):
            best_assembly_gtf=doc_summary
            max_contig_n50_ref = contig_n50

    
    if not best_assembly:
        print(f"No suitable assembly found for {tax_id}")
        return

    make_two=False
    if best_assembly_gtf is not None:
        if best_assembly['FtpPath_Assembly_rpt'] != best_assembly_gtf['FtpPath_Assembly_rpt']:
            print('No GTF')
            make_two=True
    
    print(f"Best assembly for {tax_id}: {best_assembly['AssemblyAccession']}")


    if make_two:
        assemblies={'best':best_assembly,'reference':best_assembly_gtf}
    elif best_assembly_gtf is None:
        assemblies={'best':best_assembly}
    else:
        assemblies={'reference':best_assembly_gtf}
        
    for d in assemblies.keys():
        print(d)
        dirname=os.path.join(str(tax_id)+"-"+sanitize_filename(scientific_name),d)
        cur_assembly=assemblies[d]
           
        # Downloading the FASTA file
        print(cur_assembly['FtpPath_Assembly_rpt'])
        stat_ftp_url = cur_assembly['FtpPath_Assembly_rpt']
        stat_name=sanitize_filename(os.path.basename(stat_ftp_url))#sanitize_filename(f"{tax_id}_assembly_report.txt")
        stat_destination = os.path.join(overall_path,dirname,stat_name)
        download_ftp_file(stat_ftp_url, stat_destination)

        fasta_ftp_url = re.sub('_assembly_report\.txt','_genomic.fna.gz',cur_assembly['FtpPath_Assembly_rpt'])
        fasta_name=sanitize_filename(os.path.basename(fasta_ftp_url))#sanitize_filename(f"{tax_id}_genomic.fna.gz")
        fasta_destination = os.path.join(overall_path,dirname,fasta_name)
        if not redownload and os.path.exists(fasta_destination):
            continue
        download_ftp_file(fasta_ftp_url, fasta_destination)
        if d != 'best':
            gtf_ftp_url = re.sub('_assembly_report\.txt','_genomic.gtf.gz',cur_assembly['FtpPath_Assembly_rpt'])
            gtf_name=sanitize_filename(os.path.basename(gtf_ftp_url))#sanitize_filename(f"{tax_id}_genomic.gtf.gz")
            gtf_destination = os.path.join(overall_path,dirname,gtf_name)
            #gtf_destination = f"{tax_id}/{tax_id}_genomic.gtf.gz"
            download_ftp_file(gtf_ftp_url, gtf_destination)
        else:
            try:
                url=f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{cur_assembly['AssemblyAccession']}/download?include_annotation_type=GENOME_GTF,SEQUENCE_REPORT&filename={cur_assembly['AssemblyAccession']}.zip"
                print(url)
                print(os.path.join(overall_path,dirname))
                download_and_unzip(url, os.path.join(overall_path,dirname),'.gtf')
            except:
                pass
                


##############################ENSEMBL 2023 Processing################################################

import requests
import os
import re
from ftplib import FTP

def sanitize_filename(filename):
    # Replace spaces with underscores and special characters with dashes
    filename = filename.replace(" ", "_")
    pattern = r'[^\w\-_.]'
    return re.sub(pattern, '-', filename)

def download_ftp_file(ftp_url, destination_filename):
    # Parse and download from FTP URL
    ftp_host, ftp_path = ftp_url.replace("ftp://", "").split("/", 1)
    with FTP(ftp_host) as ftp:
        ftp.login()
        os.makedirs(os.path.dirname(destination_filename), exist_ok=True)
        try:
            with open(destination_filename, 'wb') as local_file:
                ftp.retrbinary(f'RETR {ftp_path}', local_file.write)
        except Exception as e:
            print(f"Error downloading {ftp_url}: {e}")

def get_species_name(tax_id):
    """ Fetch species name using taxonomic ID from Ensembl """
    url = f"https://rest.ensembl.org/taxonomy/id/{tax_id}?content-type=application/json"
    response = requests.get(url)
    return response.json()['scientific_name']

def find_matching_files(ftp, path, patterns):
    """ Find files matching given patterns """
    ftp.cwd(path)
    files = ftp.nlst()
    matched_files = [file for file in files if any(pattern in file for pattern in patterns)]
    return matched_files

def download_ensembl_files(tax_id, path,redownload=False):
    species_name = get_species_name(tax_id)
    species_path = species_name.lower().replace(' ', '_')

    # Fetch genome information from Ensembl
    url = f"https://rest.ensembl.org/info/genomes/{species_path}?content-type=application/json"
    print(url)
    response = requests.get(url)
    data = response.json()
    assembly_name = data['assembly_name']

    # Prepare directory name
    dir_name = sanitize_filename(f"{tax_id}-{species_name}")
    dir_path = os.path.join(path, dir_name)
    os.makedirs(dir_path, exist_ok=True)

    # Prepare FTP paths
    ensembl_ftp_base = "ftp.ensembl.org"
    fasta_ftp_dir = f"/pub/release-110/fasta/{species_path}/dna/"
    gtf_ftp_dir = f"/pub/release-110/gtf/{species_path}/"

    with FTP(ensembl_ftp_base) as ftp:
        ftp.login()

        # Find matching FASTA and GTF files
        fa_patterns = [".dna_sm.toplevel.fa.gz"]
        gtf_patterns = [".gtf.gz"]
        unwanted_patterns = ["chr", "hapl", "abinitio"]

        fasta_files = find_matching_files(ftp, fasta_ftp_dir, fa_patterns)
        gtf_files = find_matching_files(ftp, gtf_ftp_dir, gtf_patterns)
        fasta_files = [f for f in fasta_files if not any(p in f for p in unwanted_patterns)]
        gtf_files = [f for f in gtf_files if not any(p in f for p in unwanted_patterns)]
        if not redownload and os.path.exists(os.path.join(dir_path, fasta_files[0])):
            return None
        # Download files
        if fasta_files:
            download_ftp_file(f"ftp://{ensembl_ftp_base}{fasta_ftp_dir}{fasta_files[0]}", os.path.join(dir_path, fasta_files[0]))
        if gtf_files:
            download_ftp_file(f"ftp://{ensembl_ftp_base}{gtf_ftp_dir}{gtf_files[0]}", os.path.join(dir_path, gtf_files[0]))


############################ASSESSING REFERENCE QUALITY##################################

import os
import pandas as pd
from Bio import SeqIO
import gtfparse
import re
import gzip

def is_gzipped(file_path):
    """ Check if the file is gzipped """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def calculate_N50(fasta_path):
    file_open_func = gzip.open if is_gzipped(fasta_path) else open
    with file_open_func(fasta_path, "rt") as handle:
        lengths = [len(rec) for rec in SeqIO.parse(handle, "fasta")]
    total_len = sum(lengths)
    lengths.sort(reverse=True)
    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum >= total_len / 2:
            return length
    return 0

"""
def count_genes_transcripts(gtf_path):
    db = gffutils.create_db(gtf_path, dbfn=":memory:", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    num_genes = len(set(db.features_of_type('gene_id')))
    num_transcripts = len(set(db.features_of_type('transcript_id')))
    return num_genes, num_transcripts
"""

def count_genes_transcripts(gtf_path):
    try:
        gtf_df = gtfparse.read_gtf(gtf_path)
        num_genes = len(gtf_df['gene_id'].unique())
        try:
            num_transcripts = len(gtf_df['transcript_id'].unique())
        except:
            num_transcripts = num_genes
        return num_genes, num_transcripts
    except:
        return 0,0

def find_files_in_directory(directory, file_extension):
    for file in os.listdir(directory):
        if re.search(file_extension, file):
            return file
    return None

def process_directory(dir_path,gtf=False):
    for root, dirs, files in os.walk(dir_path):
        for dir_name in dirs:
            #try:
            dir_full_path = os.path.join(root, dir_name)
            #taxid, species_name = dir_name.split('-', 1)
            try:
                taxid, species_name =re.search('([0-9]+-[A-Za-z_]+)',dir_full_path).group(0).split('-',1)
                common_name = species_name.replace("_", " ")
            except:
                continue
            print(species_name)
            # Find the FASTA and GTF files
            fasta_file = find_files_in_directory(dir_full_path, r'\.f(n)?a(\.gz)?$')
            print(fasta_file)
            if not fasta_file:
                continue

            fasta_path = os.path.join(dir_full_path, fasta_file)

            # Read assembly name for NCBI
            assembly_name = fasta_file
            assembly_report_path = os.path.join(dir_full_path, f"{taxid}_assembly_report.txt")
            if os.path.exists(assembly_report_path):
                with open(assembly_report_path) as f:
                    for line in f:
                        if 'Assembly name' in line:
                            assembly_name = line.split(':')[1].strip()
                            break
                      
            n50 = calculate_N50(fasta_path)
            gtf_file = find_files_in_directory(dir_full_path, r'\.gtf')
            if gtf and (gtf_file is not None):
                gtf_path = os.path.join(dir_full_path, gtf_file)
                num_genes, num_transcripts = count_genes_transcripts(gtf_path)
            else:
                gtf_file = gtf_file
                gtf_path = os.path.join(dir_full_path, gtf_file if gtf_file is not None else '')
                num_genes = None
                num_transcripts = None
            
            stat_dict = {
                'Species': species_name,
                'TaxID': taxid,
                'Common Name': common_name,
                'Genome Assembly Name': assembly_name,
                'N50': n50,
                'Number of Genes': num_genes,
                'Number of Transcripts': num_transcripts,
                'FASTA File': fasta_file,
                'FASTA Path': fasta_path,
                'GTF File': gtf_file,
                'GTF Path': gtf_path
            }
            yield stat_dict
            #except:
            #    yield None


import pandas as pd
from ete3 import Tree
import copy

def get_genera_from_species_list(species_list):
    """Extract genera from a list of species' scientific names."""
    return set(species.split()[0] for species in species_list)

def prune_tree(newick_tree_path, genera):
    """Prune the tree to retain only species from the specified genera."""
    tree = Tree(newick_tree_path,format=1)
    original_tree = copy.deepcopy(tree)
    for leaf in tree.iter_leaves():
        genus = leaf.name.split()[0]
        if genus not in genera:
            leaf.delete()
    return tree,original_tree

def get_species_not_in_tree(species_list, tree):
    # Extract species names from the tree
    species_in_tree = set(leaf.name for leaf in tree.iter_leaves())

    # Find species not in the tree
    species_not_in_tree = [species for species in species_list if species not in species_in_tree]

    return species_not_in_tree

def add_node_to_common_ancestor(tree, new_node_name, existing_node_names):
    all_nodes = [node.name for node in tree.traverse()]
    if new_node_name in all_nodes:
        print(new_node_name,'is in the tree already!')
        return None
    # Find the common ancestor of the existing nodes
    common_ancestor = tree.get_common_ancestor(existing_node_names)

    # Calculate the mean branch length to the common ancestor
    branch_lengths = [tree.get_distance(name,common_ancestor) for name in existing_node_names]
    mean_branch_length = sum(branch_lengths) / len(branch_lengths)

    # Create the new node and add it to the common ancestor
    new_node = common_ancestor.add_child(name=new_node_name, dist=mean_branch_length)
    
    return new_node

def filter_genome_dataframe(df):
    # Sort the dataframe first by 'Species', then by 'Number of Transcripts' (descending), and then by 'N50' (descending)
    df_sorted = df.sort_values(by=['Species', 'Number of Transcripts', 'N50'], ascending=[True, False, False])

    # Drop duplicates, keeping the first occurrence (which is the one with the highest 'Number of Transcripts' or 'N50' in case of a tie)
    df_filtered = df_sorted.drop_duplicates(subset='Species', keep='first')

    return df_filtered

def write_tree_and_dataframe_to_file(tree, df, filename, col1, col2):
    """
    Writes a Newick tree, a space, and two DataFrame columns to a file.

    :param tree: Newick tree string.
    :param df: Pandas DataFrame.
    :param filename: Output file name.
    :param col1: First column name from DataFrame.
    :param col2: Second column name from DataFrame.
    """
    with open(filename, 'w') as f:
        # Write the Newick tree
        f.write(tree + '\n\n')  # Newick tree followed by a newline and a space

        # Write the specified DataFrame columns
        df[[col1, col2]].to_csv(f, sep='\t', index=False,header=None)


###################################Ortho Tables#################



def download_biomart_table(dataset, filename):
    #all_attributes = [x for x in dataset.attributes.keys() if 'hsapiens' in x]
    all_attributes = [
        'ensembl_gene_id',
        'external_gene_name',
        'hsapiens_homolog_associated_gene_name',
        'hsapiens_homolog_ensembl_gene',
        'hsapiens_homolog_orthology_confidence'
    ]

    response = dataset.search({
        'attributes': all_attributes
    }, header=1)  # header=1 will include the column names

    with open(filename, 'wb') as f:
        f.write('\t'.join(all_attributes).encode('ascii')+b'\n')
        for line in response.iter_lines():
            f.write(line + b'\n')

def find_matching_rows(strings, df,colnames=None):
    matches = {}
    if colnames is None:
        colnames=df.columns
    use_df=df.loc[:,colnames]
    for s in tqdm.tqdm(strings):
        mask = use_df.isin([s]).any(axis=1)
        matches[s] = df[mask]
    return matches

def get_human_orthologs(species_name, identifiers, cache_path):
    #Returns a dictionary of {original_key: dataframe of matched rows}
    server = BiomartServer("http://www.ensembl.org/biomart")
    dataset_name = f"{species_name.lower()}_gene_ensembl"
    
    # Construct cache_path based on dataset_name
    cache_name = os.path.join(cache_path,dataset_name + "_table.txt")
    
    dataset = server.datasets[dataset_name]

    # Check cache and download if necessary
    if not os.path.exists(cache_name):
        print("Cache not found. Downloading BioMart table...")
        download_biomart_table(dataset, cache_name)
        print(f"Downloaded and saved to {cache_name}")

    # Read the table into a DataFrame
    df = pd.read_csv(cache_name, sep='\t', dtype=str)
    
    # Filter using Ensembl IDs and gene symbols and add 'original_identifier' column
    #filtered_df = df[df['ensembl_gene_id'].isin(identifiers) | df['external_gene_name'].isin(identifiers)].copy()
    #filtered_df['original_identifier'] = filtered_df.apply(lambda row: row['ensembl_gene_id'] if row['ensembl_gene_id'] in identifiers else row['external_gene_name'], axis=1)
    #orthologs = {row['original_identifier']: row.drop('original_identifier').to_dict() for _, row in filtered_df.iterrows()}

    orthologs=find_matching_rows(identifiers, df,df.columns[df.columns.str.contains('gene')])        
    return orthologs



