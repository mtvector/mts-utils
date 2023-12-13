import anndata
import pandas as pd
import scanpy as sc


def rename_merged_table(df,suffix_x='_x',suffix_y='_y'):
    '''
    renames a merged dataframe with suffixes x and y to tidy duplicate columns.
    !Drops the second if they're duplicated!
    '''
    # Get columns as a list
    columns = df.columns.tolist()
    
    # Dictionary to hold the columns to be renamed
    rename_dict = {}
    
    # Set to hold names of columns to drop
    drop_columns = set()
    
    for column in columns:
        # If column ends with '_x', it takes precedence, so we rename it by stripping '_x' and drop '_y'
        if column.endswith(suffix_x):
            new_column_name = column.rstrip(suffix_x)
            rename_dict[column] = new_column_name
            y_column = new_column_name + suffix_y
            if y_column in columns:
                drop_columns.add(y_column)
        # If column ends with '_y' but '_x' is not present, we simply rename '_y' to the base name
        elif column.endswith(suffix_y) and (column.rstrip(suffix_y) + suffix_x) not in columns:
            rename_dict[column] = column.rstrip(suffix_y)
    
    # Rename columns based on rename_dict
    df.rename(columns=rename_dict, inplace=True)
    
    # Drop the '_y' columns
    df.drop(columns=list(drop_columns), inplace=True)
    
    return df

def get_mouse_human_genemapping(hom_table='/home/matthew.schmitz/Matthew/genome/archaic/mm10/HOM_AllOrganism.rpt'):
    '''
    Get the dictionary of human to mouse orthologs from the MGI homology table
    '''
    orthos=pd.read_csv(hom_table,sep='\t')
    orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
    classcounts=orthos['DB Class Key'].value_counts()
    one2one=classcounts.index[list(classcounts==2)]
    orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]
    
    htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
    mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
    genemapping=dict(zip(mtab['Symbol'].str.upper(),htab['Symbol'].str.upper()))
    return(genemapping)
