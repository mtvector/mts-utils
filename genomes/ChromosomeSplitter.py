from Bio import SeqIO
import csv
import re
from collections import defaultdict

def zcumsum(iterable):
    cumulative_sum=[0]
    total=0
    for item in iterable:
        total += item
        cumulative_sum.append(total)
    cumulative_sum=cumulative_sum+[float("inf")]
    return cumulative_sum


class ChromosomeSplitter:
    """
    Takes a dictionary of chromosomes and split locations and splits into multiple chromosomes, and corrects associated GTF file

    # Example usage
    import os
    top_path='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/Opossum/raw_files/'
    fasta_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic.fna')
    new_fasta_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic_split.fna')
    gtf_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic.gtf')
    new_gtf_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic_split.gtf')
    
    split_locations={
        'NC_077227.1': [495286654],
        'NC_077228.1': [502619815]
    }
    
    new_sequences=ChromosomeSplitter.read_fasta_and_split(fasta_file, split_locations)
    ChromosomeSplitter.write_new_fasta(new_sequences, new_fasta_file)
    ChromosomeSplitter.gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)
    
    ##or
    ChromosomeSplitter(split_locations,fasta_file,gtf_file,new_fasta_file,new_gtf_file)
    """
    def __init__(self,split_locations,fasta_file,gtf_file,new_fasta_file,new_gtf_file):
        new_sequences=self.read_fasta_and_split(fasta_file, split_locations)
        self.write_new_fasta(new_sequences, new_fasta_file)
        self.gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)
        
    @staticmethod
    def read_fasta_and_split(file_path, split_locations):
        """
        Reads in a fasta and split chromosomes at points specicified indicated in split_locations dict
        """
        new_sequences={}
        for record in SeqIO.parse(file_path, "fasta"):
            chrom=record.id
            if chrom in split_locations:
                starts=[0] + split_locations[chrom]
                ends=split_locations[chrom] + [len(record.seq)]
                for i, (start, end) in enumerate(zip(starts, ends)):
                    new_id=f"{chrom}__part{i+1}"
                    new_sequences[new_id]=record.seq[start:end]
            else:
                new_sequences[chrom]=record.seq
        return new_sequences

    @staticmethod
    def write_new_fasta(new_sequences, output_file):
        """
        Writes a seqio genome object to fasta
        """
        with open(output_file, "w") as output_handle:
            for chrom, seq in new_sequences.items():
                SeqIO.write(SeqIO.SeqRecord(seq, id=chrom, description=""), output_handle, "fasta")
    
    @staticmethod
    def gtf_chrom_name_change(gtf_file, split_locations, output_file):
        """
        Replaces split chromosome names in associated gtf and corrects their start and end point
        """
        with open(gtf_file) as gtf_handle, open(output_file, "w") as out_handle:
            gtf_reader=csv.reader(gtf_handle, delimiter="\t")
            gtf_writer=csv.writer(out_handle, delimiter="\t")
            split_starts={x:zcumsum(split_locations[x]) for x in split_locations.keys()}
            
            for row in gtf_reader:
                if row[0].startswith('#'):
                    gtf_writer.writerow(row)
                    continue
                chrom=row[0]
                start=int(row[3])
                end=int(row[4])
    
                if chrom in split_locations.keys():
                    for idx in range(len(split_starts[chrom])-1):
                        if start < split_starts[chrom][idx+1] and start > split_starts[chrom][idx]:
                            new_chrom_name=f'{chrom}__part{idx+1}'
                            row[0]=new_chrom_name
                            row[3]=str(start - split_starts[chrom][idx])
                            row[4]=str(end - split_starts[chrom][idx])
                            break
    
                gtf_writer.writerow(row)


