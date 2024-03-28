import pandas as pd
import gtfparse
import polars as pl
import os
import csv

def gtfparse_gtf_file(file_path):
    '''Read the GTF file into a Polars DataFrame'''
    return gtfparse.read_gtf(file_path)

def polars_to_pandas(df_polars):
    '''Convert Polars DataFrame to pandas DataFrame'''
    return df_polars.to_pandas()

def gtf_add_missing_transcripts(df):
    '''creates transcript annotations for each exon annotated to transcript, encompassing all assigned exons. for gtfs missing transcript entries'''
    # Group by gene_id and transcript_id
    grouped = df.groupby(['gene_id', 'transcript_id'])

    new_rows = []
    for (gene_id, transcript_id), group in grouped:
        if transcript_id and 'transcript' not in group['feature'].values:
            # Calculate the range for the transcript
            start = group['start'].min()
            end = group['end'].max()

            # Get the first row of the group as a template for new row
            template_row = group.iloc[0].copy()

            # Modify the template row for the transcript
            template_row['feature'] = 'transcript'
            template_row['start'] = start
            template_row['end'] = end

            new_rows.append(template_row)

    # Append new rows to the original DataFrame
    new_df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    # Sort the DataFrame based on start position
    new_df.sort_values(by=['seqname', 'start'], inplace=True)

    return new_df

def gtf_attribute_string_generator(columns_to_condense):
    '''combines specified columns into gtf/gff2 attributes format'''
    def format_string(row):
        return '; '.join([f'{col} "{row[col]}"' for col in columns_to_condense if row[col] != ''])
    return format_string

def gtf_df_sort(df):
    '''sort a gtf df in the format required by cellranger, gene-transcript, exon nested and sorted'''
    def feature_sort(row):
        if row['feature'] == 'gene':
            return 1
        elif row['feature'] == 'transcript':
            return 2
        elif row['feature'] == 'exon':
            return 3
        else:
            return 4
    
    # Apply the custom sorting function
    df['feature_sort'] = df.apply(feature_sort, axis=1)
    
    # Sort by gene_id, feature_sort, and then start position
    df_sorted = df.sort_values(by=['gene_id','transcript_id','feature_sort','start','end'])
    
    # Drop the temporary sorting column
    df_sorted = df_sorted.drop('feature_sort', axis=1)
    return df_sorted

def write_gtf_df(df, output_file_path):
    '''write gtf df columns to gtf file format'''
    if len(df.columns)>8:
        if 'attribute' in df.columns:
            df.drop('attribute',axis=1,inplace=True)
        format_string=gtf_attribute_string_generator(df.columns[8:])
        df.loc[:,'attribute']=df.apply(format_string, axis=1)
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    #####FOR SOME REASON GTF MUST BE IN ASCII ENCODING FOR CELLRANGER!!!!!!!############
    df.loc[:,columns].to_csv(output_file_path, sep='\t', index=False,
                             header=False,quoting=csv.QUOTE_NONE,encoding='ascii')

#Example usage of above functions
'''
top_path='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/African_green_monkey/fixed_genome/ncbi_dataset/data/GCF_015252025.1'
file_path = os.path.join(top_path,'genomic.gtf')
output_file_path = os.path.join(top_path,'genomic_transcripts.gtf')

polars_df = gtfparse_gtf_file(file_path)
pandas_df = polars_to_pandas(polars_df)
df_with_transcripts = gtf_add_missing_transcripts(pandas_df)
df_with_transcripts=df_with_transcripts.loc[~df_with_transcripts['transcript_id'].str.contains('unknown'),:]
df_with_transcripts = gtf_df_sort(df_with_transcripts)
df_with_transcripts['score']='.'
write_gtf_df(df_with_transcripts, output_file_path)
'''



