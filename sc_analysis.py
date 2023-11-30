#General scRNAseq analysis functions
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn
import anndata
import scipy
import tqdm
import matplotlib
import matplotlib.pyplot as plt

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


def select_k_cells(adata, obs_column, k):
    """
    Select k cells from an AnnData object with as equal representation from each category in the obs column as possible.

    Parameters:
    - adata: AnnData object
    - obs_column: The column in adata.obs to use for categories
    - k: Total number of cells to select

    Returns:
    - An AnnData object containing the selected cells
    """
    
    categories = adata.obs[obs_column].value_counts().index
    n_categories = len(categories)
    base_num = k // n_categories

    cells_to_take = {}
    selected_indices = []
    leftovers = 0

    # First pass: Allocate cells according to min(base_num, available_cells_in_category)
    for cat in categories:
        cat_cells = adata.obs.index[adata.obs[obs_column] == cat].tolist()
        if len(cat_cells) <= base_num:
            selected_from_cat = cat_cells
            leftovers += base_num - len(cat_cells)
        else:
            selected_from_cat = np.random.choice(cat_cells, size=base_num, replace=False)
        selected_indices.extend(selected_from_cat)

    # Second pass: Distribute the leftovers
    for cat in categories:
        cat_cells = list(set(adata.obs.index[adata.obs[obs_column] == cat].tolist()) - set(selected_indices))
        if leftovers > 0 and len(cat_cells) > 0:
            additional = min(leftovers, len(cat_cells))
            selected_additional = np.random.choice(cat_cells, size=additional, replace=False)
            selected_indices.extend(selected_additional)
            leftovers -= additional

    return adata[selected_indices].copy()


def sanitize_anndata(adata):
    '''
    Fix for saving anndata objects with nans
    '''
    def convert_columns(df):
        for column in df.columns:
            df[column].astype('string')
        df=df.convert_dtypes()
        #df[df.select_dtypes(['object']).columns] = df[df.select_dtypes(['object']).columns].astype('string')
        for column in df.columns:
            df[column] = df[column].fillna(np.nan)
            # Check if column is of type 'object' or 'category' with object categories
            if (df[column].dtype == 'object') or (df[column].dtype == 'string') or (df[column].dtype.name == 'category' and issubclass(df[column].cat.categories.dtype.type, object)):
                # Use to_numeric for conversion, set errors='coerce' to handle non-convertible values
                converted_column = pd.to_numeric(df[column], errors='coerce')
    
                # Check if the conversion was successful (i.e., not all values are NaN)
                if converted_column.notna().any():
                    df[column] = converted_column
                else:
                    # Convert to string if conversion to float fails
                    df[column] = df[column].astype(str).fillna('nan')
            if 'float' in df[column].dtype.name.lower():
                df[column]=df[column].astype('string').fillna('nan').astype('float')
        return df

    adata.obs = convert_columns(adata.obs)
    adata.var = convert_columns(adata.var)
    adata.obs.index=adata.obs.index.astype(str)
    adata.var.index=adata.var.index.astype(str)
    return adata


def harmonize_and_concatenate(adatas,label=None):
    """
    Concatenate a list of AnnData objects after harmonizing their obs columns.

    Parameters:
    - adatas: List of AnnData objects

    Returns:
    - An AnnData object resulting from the concatenation
    """
    
    # Identify all unique columns across the AnnData obs
    all_columns = set()
    for adata in adatas:
        all_columns.update(adata.obs.columns)
    
    # Ensure each AnnData obs has every column, filling with 'nan' if not
    for adata in adatas:
        for col in all_columns:
            if col not in adata.obs.columns:
                adata.obs[col] = np.nan
                   
    adatas=[sanitize_anndata(x) for x in adatas]
    # Concatenate the AnnData objects
    if label is not None:
        combined = anndata.concat(adatas, join='outer',label=label, merge="first")
    else:
        combined = anndata.concat(adatas, join='outer', merge="first")
    return combined


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

#Keeps the rows in the same order, reorders columns
def maximize_diagonality(matrix, n_clusters=None):
    '''
    sorts the columns of the matrix to maximize the diagonality of the matrix.
    '''
    reordered_rows=list(range(matrix.shape[0]))
    reordered_cols=np.argsort(np.argmax(matrix,0))
    optimized_matrix=matrix[:,reordered_cols]   
    return optimized_matrix, reordered_rows, reordered_cols

def get_cross_category_jaccard(df,obs1,obs2,normalize_jaccard=True,diagonalize=True):
    '''
    Get the jaccard matrix from adataframe with two categorical columns obs1 and obs2
    Normalize jaccard is false removes denominator of jaccard index calculation.
    Diagonalize sorts the columns of the matrix to maximize the diagonality of the matrix.
    '''
    df[obs1]=df[obs1].astype('category')
    df[obs2]=df[obs2].astype('category')
    pca_oh=one_hot(df[obs1].cat.codes)
    pca_categories=df[obs1].cat.categories
    scvi_oh=one_hot(df[obs2].cat.codes)
    scvi_categories=df[obs2].cat.categories
    if normalize_jaccard:
        jac_matrix=jaccard_index_matrix(pca_oh,scvi_oh)
    else:
        jac_matrix=pca_oh.T@scvi_oh
    if diagonalize:
        jac_diag=maximize_diagonality(jac_matrix)
        jac_out=pd.DataFrame(jac_diag[0])
        jac_out.index=pca_categories[jac_diag[1]]
        jac_out.columns=scvi_categories[jac_diag[2]]
    else:
        jac_out=pd.DataFrame(jac_matrix,index=pca_categories,columns=scvi_categories)
    return(jac_out)

def create_heatmap(df, **kwargs):
    '''
    Resize text to make heatmap size
    Determine the number of columns and rows in the DataFrame
    '''
    num_cols = len(df.columns)
    num_rows = len(df.index)
    
    # Estimate font size based on the number of columns and rows
    font_size = min(20, max(8, int(400 / (num_cols + num_rows))))
    
    # Estimate figure size based on the number of columns and rows
    fig_size = (min(max(6, num_cols * 0.6),48), min(max(6, num_cols * 0.6),48))
    
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Create the heatmap
    seaborn.heatmap(df, ax=ax, annot_kws={"size": font_size}, **kwargs)# annot=True, fmt=".1f"
    
    # Resize the labels
    ax.set_xticklabels(ax.get_xticklabels(), size=font_size)
    ax.set_yticklabels(ax.get_yticklabels(), size=font_size, rotation=0)
    
    # Show the plot
    #plt.show()

    
def one_hot(a,size=None):
    '''
    converts numpy vector of integers to onehot matrix
    '''
    size=a.max() if size is None else size
    a = np.array(a)
    b = np.zeros((a.size,size + 1))
    b[np.arange(a.size), a] = 1
    return b


def jaccard_index_matrix(matrix1, matrix2):
    '''
    Calculate jaccard index from two onehot matrices
    '''
    # Matrix multiplication to get the intersection
    intersection = matrix1.T@matrix2
    
    # Addition to get the union, followed by clipping at 1
    union = np.expand_dims(np.sum(matrix1, axis=0), axis=1) + np.expand_dims(np.sum(matrix2, axis=0), axis=0)
    
    # Avoid division by zero
    union = np.where(union == 0, 1, union)
    
    return intersection / union


def diagonality(matrix,ex=2):
    '''
    Diagonality score for the weight of matrix values that are on the diagonal
    '''
    m, n = matrix.shape
    if m == 0 or n == 0:
        return 0

    # Distance from a point (i, j) to the line y = (m/n)x
    def distance_to_diagonal(i, j,ex):
        return abs(m * j - n * i) / ((m**ex + n**ex)**(1/ex))

    # Sum distances for all nonzero elements
    total_distance = 0
    for i in range(m):
        for j in range(n):
            if matrix[i, j] != 0:
                total_distance += distance_to_diagonal(i, j,ex=ex)

    # Normalize the sum
    # Maximum distance is corner elements to the diagonal
    max_distance = max(distance_to_diagonal(0, n - 1,ex=ex), distance_to_diagonal(m - 1, 0,ex=ex))
    max_total_distance = max_distance * np.count_nonzero(matrix)

    return 1-(total_distance / max_total_distance if max_total_distance != 0 else 0)




def matrix_multiplication_and_rename(A, B):
    # Get the common columns between A and B
    common_cols = set(A.columns) & set(B.columns)
    
    # Restructure A and B to only have the common columns
    A_common = A[list(common_cols)]
    B_common = B[list(common_cols)]
    
    # Perform matrix multiplication
    result = np.dot(A_common.values, B_common.values.T)
    
    # Convert the result to a DataFrame
    result_df = pd.DataFrame(result, index=A.index, columns=B.index)
    
    # Rename rows and columns
    #result_df.index.name = 'Rows_from_A'
    #result_df.columns.name = 'Rows_from_B'
    
    return result_df


def crazyshuffle(arr):
    '''
    Shuffles rows relative for columns e.g. for empirical correlation significance testing
    '''
    x, y = arr.shape
    rows = np.indices((x,y))[0]
    cols = [np.random.permutation(y) for _ in range(x)]
    return arr[rows, cols]


def corr2_coeff(A, B):
    '''
    Calculate pearson correlation of row vectors of A vs row vectors of B
    '''
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)
    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))


def get_top_n_indices_and_values(df, columns, n):
    '''
    For pandas dataframes, self explanatory.
    '''
    top_indices_values = {}
    for column in columns:
        # Get the top N indices
        top_n_indices = df[column].nlargest(n).index.tolist()
        # Get the values for these indices
        top_n_values = df[column][top_n_indices].tolist()
        # Combine indices and values into a list of tuples
        top_indices_values[column] = list(zip(top_n_indices, top_n_values))
    return top_indices_values

def get_matching_indices(a, b):
    """Return the indices to reorder b to match the sequence in a."""
    # Create an index map of b
    index_map = {value: index for index, value in enumerate(b)}
    
    # Use the index map to get indices of b that match the sequence in a
    indices = [index_map[val] for val in a if val in index_map]
    
    return indices


def sum_duplicate_var(adata):
    '''
    Sums the duplicate var names of an anndata instead of making them unique.
    Very slow.
    '''
    import tqdm
    #Instead of removing duplicate vars, sum them
    #Make the first duplicate the sum of all duplicate
    for g in tqdm.tqdm(adata.var.index[adata.var.index.duplicated(keep='first')]):
        inds=np.where(adata.var.index.isin([g]))[0]
        adata.X[:,inds[0]]=adata.X[:,inds].sum(1)
    return(adata[:,~adata.var.index.duplicated(keep='first')])

"""
def sum_duplicate_var(adata):
    #adata.var.sort_index(inplace=True)
    adata.X = adata.X[:, adata.var.index.argsort()]

    # Identify the duplicated indices
    duplicated = adata.var.index.duplicated(keep='first')
    unique_duplicated = adata.var.index[duplicated].unique()

    # Loop through each unique duplicated index and sum the columns
    for gene in tqdm.tqdm(unique_duplicated):
        # Find the columns that need to be summed
        cols = np.where(adata.var.index == gene)[0]
        summed_col = adata.X[:, cols].sum(axis=1)

        # Replace the first occurrence with the sum
        adata.X[:, cols[0]] = summed_col

    # Drop all but the first of each duplicated column
    adata = adata[:, ~duplicated]
    return(adata)
"""

def most_frequent(List): 
    return max(set(List), key = List.count)


def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

def cell_cycle_score(adata,plot=False):
    '''
    Calculates the cell cycle score using cell cycle stage markers
    '''
    s_genes=['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes,g2m_genes=g2m_genes)
    if plot:
        sc.pl.violin(adata, ['G2M_score','S_score'], groupby='leiden')
        sc.pl.umap(adata, color=['G2M_score','S_score','phase','leiden'],save='_cc')
    return(adata)


def marker_analysis(adata,variables=['leiden','region'],markerpath='https://docs.google.com/spreadsheets/d/e/2PACX-1vTz5a6QncpOOO-f3FHW2Edomn7YM5mOJu4z_y07OE3Q4TzcRr14iZuVyXWHv8rQuejzhhPlEBBH1y0V/pub?gid=1154528422&single=true&output=tsv',subclass=True,plotall=True,save=False,prefix=''):
    '''
    Old. Gets marker list scores from personal marker list using scanpy score_genes
    '''
    sc.set_figure_params(color_map="Purples")
    import random
    markerpath=os.path.expanduser(markerpath)
    markers=pd.read_csv(markerpath,sep="\t")
    markers[markers.keys()[0]]=[str(x) for x in markers[markers.keys()[0]]]
    markers[markers.keys()[2]]=[str(x).split(',') for x in markers[markers.keys()[2]]]
    markers[markers.keys()[3]]=[str(x).split(';') for x in markers[markers.keys()[3]]]
    markers[markers.keys()[3]]=[[str(x).split(',') for x in y] for y in markers[markers.keys()[3]]]
    uniqueClasses=set([y for x in markers[markers.keys()[2]] for y in x if y!='nan'])
    uniqueSubClasses=set([z for x in markers[markers.keys()[3]] for y in x for z in y if z!='nan'])
    comboClasses=[]
    #print(markers)
    markers[markers.keys()[2]]=[[y.lstrip() for y in x ]for x in markers[markers.keys()[2]]]
    markers[markers.keys()[3]]=[[[z.lstrip() for z in  y] for y in x] for x in markers[markers.keys()[3]]]
    if subclass:
        for i in range(markers.shape[0]):
            rowlist=[]
            for j in range(len(markers[markers.keys()[2]][i])):
                for k in markers[markers.keys()[3]][i][j]:
                    rowlist.append(re.sub('^ ','',' '.join(filter(lambda x: x != 'nan',[k,markers[markers.keys()[2]][i][j]]))))
            comboClasses.append(rowlist)
    else:
        for i in range(markers.shape[0]):
            rowlist=[]
            for j in range(len(markers[markers.keys()[2]][i])):
                rowlist.append(markers[markers.keys()[2]][i][j])
            comboClasses.append(rowlist)

    markers['fullclass']=comboClasses
    markers.set_index(markers.keys()[0],inplace=True,drop=False)
    markers=markers.loc[ [x for x in markers[markers.keys()[0]] if x in adata.var_names],:]
    uniqueFullClasses=set([y for x in markers['fullclass'] for y in x if y!='nan'])
    from collections import defaultdict
    markerDict = defaultdict(list)

    for x in uniqueFullClasses:
        for y in markers[markers.keys()[0]]:
            if x in markers.loc[y,'fullclass']:
                markerDict[x].append(y)
    markerDictClass = defaultdict(list)
    for x in uniqueClasses:
        for y in markers[markers.keys()[0]]:
            if x in markers.loc[y,'fullclass']:
                markerDictClass[x].append(y)

    markerPlotGroups=[]
    for k in markerDict.keys():
        if len(markerDict[k])>1:
            sc.tl.score_genes(adata,gene_list=markerDict[k],score_name='mkrscore'+k,gene_pool= markerDict[k]+random.sample(adata.var.index.tolist(),min(1000,adata.var.index.shape[0])),use_raw=False)
            markerPlotGroups.append('mkrscore'+k)
    adata.uns['marker_groups']=list(markerDict.keys())
    for tag in variables:
        pd.DataFrame(adata.obs.groupby(tag).describe()).to_csv(os.path.join(sc.settings.figdir, tag+prefix+"MarkerSumStats.csv"))

    sc.pl.umap(adata, color=markerPlotGroups,save=prefix+"_Marker_Group")

    #sc.pl.violin(adata, markerPlotGroups, groupby='leiden',save=prefix+"_Marker_Group_violins")
    if plotall:
        for i in markerDictClass:
            genes=sorted(markerDictClass[i])
            genes=[g for g in genes if g in adata.var.index]
            sc.pl.umap(adata, color=genes,save=prefix+"_"+str(i)+"_Marker",use_raw=False)
    #adata.uns['markers']=markers
    adata.obs['subclassname']=[re.sub('mkrscore','',x) for x in adata.obs.loc[:,['mkrscore' in x for x in adata.obs.columns]].astype('float').idxmax(axis=1)]
    def most_frequent(List): 
        return max(set(List), key = List.count) 
    classlist=[]
    for c in adata.obs['subclassname']:
        fullbool=[c in x for x in markers['fullclass']]
        flatclass=[item for sublist in markers.loc[fullbool,'type'] for item in sublist]
        classlist.append(most_frequent(flatclass))
    adata.obs['classname']=classlist
    sc.pl.umap(adata,color=['classname'],save='classname')
    sc.pl.umap(adata,color=['subclassname'],save='subclassname')
    #adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    return(adata)


import numba
@numba.njit()
def abscorrelation(x, y):
    '''
    Calculate abs correlation distance (used for X.T nn clustering)
    '''
    mu_x = 0.0
    mu_y = 0.0
    norm_x = 0.0
    norm_y = 0.0
    dot_product = 0.0

    for i in range(x.shape[0]):
        mu_x += x[i]
        mu_y += y[i]

    mu_x /= x.shape[0]
    mu_y /= x.shape[0]

    for i in range(x.shape[0]):
        shifted_x = x[i] - mu_x
        shifted_y = y[i] - mu_y
        norm_x += shifted_x ** 2
        norm_y += shifted_y ** 2
        dot_product += shifted_x * shifted_y

    if norm_x == 0.0 and norm_y == 0.0:
        return 0.0
    elif dot_product == 0.0:
        return 1.0
    else:
        return 1.0 - np.absolute(dot_product / np.sqrt(norm_x * norm_y))


def op_by_obs(adata, obs_col, op=lambda x: x.mean(axis=0)):
    '''
    Groups anndata by an obs column, then gets the mean value of adata.X by column
    modified from Valentine_Svensson
    '''
    groupby_object = adata.obs.groupby([obs_col], observed = True)
    
    X = adata.X
    
    N_obs = groupby_object.ngroups
    N_var = X.shape[1]
    X_operated = scipy.sparse.lil_matrix((N_obs, N_var))
    
    group_names = []
    index_names = []
    row = 0
    for group_columns, idx_ in tqdm.tqdm(groupby_object.indices.items()):
        X_operated[row] = op(X[idx_])
        row += 1
        group_names.append(group_columns)
        index_names.append(group_columns)
        #index_names.append('_'.join(map(str, group_columns)))
    
    if scipy.sparse.isspmatrix_csr(X):
        X_operated = X_operated.tocsr()
    else:
        X_operated = np.array(X_operated.todense())
    
    obs = pd.DataFrame(group_names, columns=[obs_col], index=index_names)
    
    new_adata = anndata.AnnData(X=X_operated, obs=obs, var=adata.var)
    return new_adata

def standard_scale(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata,n_top_genes=10000,subset=False)
    sc.pp.scale(adata, max_value=10)

def gini_coefficient(counts):
    """Calculate the Gini coefficient for a vector of counts."""
    total = sum(counts)
    if total == 0:
        return 0
    squared_sum = sum((count / total) ** 2 for count in counts)
    return squared_sum

def calculate_function_per_cluster(df, cluster_col, individual_col,fun=gini_coefficient):
    """Calculate a function (originally gini coefficient) for each cluster and return as a DataFrame."""
    gini_scores = []
    grouped = df.groupby(cluster_col)
    for cluster, group in grouped:
        counts = group[individual_col].value_counts()
        gini_score = fun(counts)
        gini_scores.append({'Cluster': cluster, 'Value': gini_score, 'Type': cluster_col})
    return pd.DataFrame(gini_scores)

