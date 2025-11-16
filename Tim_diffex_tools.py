"""
A collection of functions I used to perform my analysis. 
These were called within jupyter notebooks, but defined here to make the notebooks more concise. 
"""

# Standard libraries
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
import seaborn as sns
import numpy as np 
import scanpy as sc
import re

# diffex libraries 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
from gseapy.plot import gseaplot
from sanbomics.tools import id_map
from sanbomics.plots import volcano

# PCA libraries 
from sklearn.decomposition import PCA
from sklearn import preprocessing

from scipy.stats import spearmanr


# ==============
# String formatting functions
# ==============
def cleanname(colname):
    """
    A function which cleans the column names in our matrix
    This is specific to the RNA-seq analysis.... dont use in others...
    """
    match = re.search(r"_\d\b", colname)    
    if match:
        clean = colname.replace(match.group(),"").replace("n","-").replace("p","+").replace("_","") 
        if clean == "DOX-DPI-":
            clean = "Control"
        return clean
    else:
        return colname
    

def peptide_namecleaner(peptide_name, pri=False):
    """
    Cleans peptide names before they are displayed on plots. 
    Prints Changes to quality check
    """    
    istailmod = re.search("_[A-Z]{5,}",peptide_name) # <-- determines if this is a tail modification
    if istailmod:           # if it is, then we go and replace that part. 
        peptide_name = peptide_name.replace(istailmod.group(),"_TailVar")
    
    
    removal_part = re.search(pattern="_[0-9]{,3}_[0-9]{,3}_",string=peptide_name)
    clean_pepname = peptide_name.replace(removal_part.group(),"_")
    if pri:
        print(f"{peptide_name} --> {clean_pepname}")
    return clean_pepname




# ==============
# RNA-seq differential expression analysis
# ==============
def diffex_pipeline(counts,metadata,design_contrast,treatment,control, species="human"):
    """ 
    A function which automatically runs Deseq2 for our analysis, making the jupyter notebook cleaner. 

    Args:
        - Counts table (samples = rows, genes = cols)
        - Metadata table (Samples = rows, metadata = cols)
        - design_contrast - the header from the metadata table we want to focus on. (will be input to design contrast and other. )

    Returns: 
        - 
    """
    # 0. safety checks
    if design_contrast not in metadata.columns:
        print("Invalid design_contrast selected")
        return
    
    elif treatment not in list(metadata[design_contrast]) or control not in list(metadata[design_contrast]):
        print("Invalid treatment/control selected")
        return  
    
    # 1. Filter out 0s from the tranposed table 
    nonzero=(counts.sum(axis=0) >0)
    counts=counts[nonzero[nonzero].index]

    # 2. Create dds object
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors=design_contrast
    )

    # 3. Run Deseq2
    # this adds a bunch of stats to our dds object. 
    dds.deseq2()

    # 4. We now create a Deseq_Stats object  
    stat_res = DeseqStats(
        dds,
        contrast=[design_contrast,treatment,control]
    )
    # We then run the calculation
    stat_res.summary()

    # we then call the results dataframe from stat_res
    res = stat_res.results_df

    
    mapper = id_map(
        species=species
    )
    res["Symbol"] = res.index.map(mapper.mapper)

    return dds, res



def asterisker(x):
    """ 
    Places significance *** at the correct positions in the heatmap
    """
    asterisk = ''
    if x <=0.001:
        asterisk = "***"
    elif x <=0.01:
        asterisk = "**"
    elif x <=0.05:
        asterisk = "*"
        
    return asterisk 


def asterisk_fcheatmap(pvals,fcs,vmax=2,vmin=-2):
    """
    Args:
        - a dataframe of p-values
        - a dataframe of fold-changes
    """
    sigmap=pvals.applymap(asterisker)
    plt.figure(figsize=(5,len(pvals)/2))
    sns.heatmap(
        data=fcs,
        cmap="icefire",
        center=0,vmin=vmin,vmax=vmax,
        annot=sigmap,fmt="s"
    )


def plot_genesets(df,geneset_paths,symbolcol="Symbol",fc_col="log2FoldChange", sigcol="padj", siglevel=0.05):
    """
    Reads in gene set .tsv files from the GSEA website and plots them for our data
    """
    for path in geneset_paths:
        # Retrieve the genes
        gene_file = pd.read_csv(fr"{path}",delimiter="\t")
        geneset_name= gene_file.columns[1]
        geneset=gene_file.set_index("STANDARD_NAME").loc["GENE_SYMBOLS"].values
        geneset = geneset[0].split(",")

        data=df[df[symbolcol].isin(geneset)].set_index(symbolcol)


        asterisk_fcheatmap(
            pvals=data[["padj"]],
            fcs=data[["log2FoldChange"]]
                   )

        plt.title(geneset_name)
        plt.show()
    return




# ==============
# PCA functions
# ==============
def pc_autom(t_dataset, metadata=None,mergetype="left"):
    """ 
    Automatically scales and fits a PCA model on the dataset
    returns that PCA model, and the dataset containing mappings for the principal components

    Optional** 
    Can also take a metadata table. 
    The index of this table must match up to the index of the dataset
    """
    scaled = preprocessing.scale(t_dataset)
    pca = PCA()
    pca.fit(scaled)
    pca_data = pca.transform(scaled)

    # reformat the dataframe
    data = pd.DataFrame(pca_data)
    data.index = t_dataset.index
    data.columns = ["PC"+str(i) for i in list(range(1,len(data.columns)+1))]

    if metadata is not None:
        # By default we merge left, but can pass a different arg in if you want.
        data = data.join(metadata, how=mergetype)

    return pca, data


def screeplot(pca):
    plt.figure()
    explained_var = pca.explained_variance_ratio_

    var_ex = np.round(
        explained_var*100,
        decimals=1
    )
    sns.barplot(y=var_ex,
                x=list(range(1,len(var_ex)+1)))
    plt.ylabel("% variation")
    plt.xlabel("principal component")
    plt.title("Scree plot")
    plt.show()


def comp_contributions(pca,t_data,components=[1],num_shown=False):
    """
    Requires as input, PCA, T_data and the number of components you want (default = only PC1)
    """
    # -1 from all the components to properly index PCA
    components = [i-1 for i in components]

    contribution_list = []
    for comp in components:
        pc_contributions=pd.Series(
            pca.components_[comp],
            index=t_data.columns
        )
        contribution_list.append(pc_contributions)

    # Concatenate all of these into a Dataframe 
    conc = pd.DataFrame(pd.concat(contribution_list,axis=1))
    

    conc=conc
    summed_contributions = conc.apply(abs).sum(axis=1).sort_values(ascending=False)

    
    if num_shown:
        show=num_shown
    
    else:
        show=len(summed_contributions)

    fig = plt.figure(figsize=(8,7))
    sns.barplot(
        x=summed_contributions.values[:show],
        y=map(peptide_namecleaner,summed_contributions.index[:show]) 
    )

    return summed_contributions, conc, fig