import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import re
import numpy as np



# statistical packages 
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_ind

from sklearn.decomposition import PCA
from sklearn import preprocessing



def load_rename(EpiProfile_output,SamplesToAbbrev):
    """
    Automates the process of loading and formatting the dataset
    Inputs:
    - The EpiProfile Output file
    - The Samples to Abbreviation file 
    """
    xic_areas=pd.read_csv(EpiProfile_output,index_col=0)

    # Select ONLY the columns which contain Area 
    mask= xic_areas.iloc[0]=="Area"
    area_cols = xic_areas.columns[mask]
    xic_areas =xic_areas[area_cols]

    # Remove the first row (this contains only the word area and we have already selected using it)
    xic_areas = xic_areas.iloc[1:]
    # Drop all rows that contain only Nan values - this will de-select all the rows that specify 
    # the peptide type (e.g., TKQTAR(H3_3_8)) and also peptides that were not idetified in any sample
    xic_areas = xic_areas.dropna(axis=0)

    # remove any potential trailing characters pandas may have added on....
    xic_areas.columns = [re.sub(r"\.\d$", "", col) for col in xic_areas.columns] 

    # Rename the samples using the excel file that maps the names.
    samples = pd.read_excel(SamplesToAbbrev)
    samples = {sample:abbrev for sample, abbrev in zip(samples["Sample_Names"],samples["Abbreviated_names"])}
    xic_areas.columns = [samples[col] for col in xic_areas.columns]

    return xic_areas






def combine_similar(histone_data):
    # I am personally skeptical over whether this is a necessary step
    # as in my experience EpiProfile was reasonably good at differentiating 
    # these mofifications (data unpublished). Also Benjamin Garcia's lab often 
    # doesnt do this step
    # however you may choose to do so.
    combin_mods =    {
    ('H3_9_17_K9ma', 'H3_9_17_K14ma'): 'H3_9_17_K9/K14ma',
    ('H3_9_17_K9ac', 'H3_9_17_K14ac'): 'H3_9_17_K9/K14ac',
    ('H3_9_17_K9pr', 'H3_9_17_K14pr'): 'H3_9_17_K9/K14pr',
    ('H3_9_17_K9gl', 'H3_9_17_K14gl'): 'H3_9_17_K9/K14gl',
    ('H3_9_17_K9cr', 'H3_9_17_K14cr'): 'H3_9_17_K9/K14cr',
    ('H3_9_17_K9su', 'H3_9_17_K14su'): 'H3_9_17_K9/K14su',
    ('H3_9_17_K9bu', 'H3_9_17_K14bu'): 'H3_9_17_K9/K14bu',
    ('H3_18_26_K18ac', 'H3_18_26_K23ac'): 'H3_18_26_K18/K23ac',
    ('H3_18_26_K18pr', 'H3_18_26_K23pr'): 'H3_18_26_K18/K23pr',
    ('H3_18_26_K18cr', 'H3_18_26_K23cr'): 'H3_18_26_K18/K23cr',
    ('H3_18_26_K18gl', 'H3_18_26_K23gl'): 'H3_18_26_K18/K23gl',
    ('H3_18_26_K18hi', 'H3_18_26_K23hi'): 'H3_18_26_K18/K23hi',
    ('H3_18_26_K18ma', 'H3_18_26_K23ma'): 'H3_18_26_K18/K23ma',
    ('H3_18_26_K18su', 'H3_18_26_K23su'): 'H3_18_26_K18/K23su',
    ('H3_18_26_K18bu', 'H3_18_26_K23bu'): 'H3_18_26_K18/K23bu',
    ('H3_27_40_K27hi', 'H3_27_40_K36hi'): 'H3_27_40_K27/K36hi',
    ('H3_27_40_K27ma', 'H3_27_40_K36ma'): 'H3_27_40_K27/K36ma',
    ('H3_27_40_K27cr', 'H3_27_40_K36cr'): 'H3_27_40_K27/K36cr',
    ('H3_27_40_K27bu', 'H3_27_40_K36bu'): 'H3_27_40_K27/K36bu',
    ('H3_27_40_K27pr', 'H3_27_40_K36pr'): 'H3_27_40_K27/K36pr',
    ('H2A1_4_11_K5ac', 'H2A1_4_11_K9ac'): 'H2A1_4_11_K5/K9ac',
    ('H2AX_4_11_K5ac', 'H2AX_4_11_K9ac'): 'H2AX_4_11_K5/K9ac',
    ('H2AJ_4_11_K5ac', 'H2AJ_4_11_K9ac'): 'H2AJ_4_11_K5/K9ac'
    }
    # combine these mods 
    for mod_pair in combin_mods.keys():
        if mod_pair[0] in histone_data.index and mod_pair[1] in histone_data.index:
            combined_mod = histone_data.loc[list(mod_pair)].sum()
            combined_mod.name = combin_mods[mod_pair]

            histone_data=histone_data.drop(list(mod_pair))
            histone_data=pd.concat([histone_data,combined_mod.to_frame().T])

    # now average the H4 peptides 
    h4 = histone_data[histone_data.index.str.contains(r"H4.*4.*17.*ac",regex=True)].copy()
    h4["Acetyl_count"] = h4.index.str.count("ac")
    h4 = h4.groupby("Acetyl_count").sum().reset_index().sort_values("Acetyl_count").set_index("Acetyl_count")
    h4.index= ['H4_4_17_1ac', 'H4_4_17_2ac', 'H4_4_17_3ac', 'H4_4_17_4ac']
    histone_data = pd.concat(
        [
        h4,
        histone_data[~histone_data.index.str.contains(r"H4.*4.*17.*ac", regex=True)]
        ],
    )

    return histone_data



def peptide_namecleaner(pep_name, rem_number=False):
    """
    Reformats the peptide strings presented by EpiProfile 
    in a manner that is better for presenting

    - For analysis you may want to leave rem_number as False
    - but for presenting the data you may want to check it as 
      True to remove the numbering which will look nicer
    """
    # Replace all " " OR "." WITH "_" 
    clean_name =re.sub(pattern=r"(\s)|(\.)", 
           repl="_", 
           string=pep_name.strip())
    
    # The tail variant names are very long so we reformat them. 
    # .+ means one or more of anything
    tailvar=re.search(pattern=r"H.+_\d+_\d+_(.+)_",
                 string=clean_name)
    if tailvar:
        clean_name = clean_name.replace(
            clean_name[clean_name.find(f"_{tailvar[1]}_"):],
            f"_Var-{tailvar[1]}"
            )

    if rem_number and not re.search(pattern="unmod",string=pep_name):   
        removal_part = re.search(pattern="_[0-9]{,3}_[0-9]{,3}_",string=clean_name)
        clean_name = clean_name.replace(removal_part.group(),"_")

    return clean_name



def signal_plot(peptide_totals):
    """
    Plots total signal per sample and also the relative abundances of 
    individual peptides allowing the user to identify samples that had 
    low overall signal or samples that had unusual abundance of a particular 
    peptide. 
    (note some peptide families will ionize with greater efficiency than others so the 
    relative abundances of peptides is not biologically relevant, however this can tell 
    you if there was a technical error)
    
    Args:
    - a dataframe with 
        index = peptide_id
        columns = sample
        values = total signal per peptide_id in that sample
    """
    fig, axs= plt.subplots(
    nrows=2,
    figsize=(len(peptide_totals.columns)/2.5,24/2.5)   # <- height will be constant, but width proportional to number of samples
    )

    # Plot the raw intensities
    peptide_totals.T.plot(
        kind="bar",stacked=True,
        edgecolor="black",
        color=sns.color_palette("colorblind",len(peptide_totals)),
        ax=axs[0]
    )
    axs[0].set_xticks([])
    axs[0].legend().remove()
    axs[0].set_ylabel("Total Area")


    # Normalise the peptide totals.
    norm = peptide_totals.T.div(
        peptide_totals.T.sum(axis=1),
        axis=0
    )
    norm.plot(
        kind="bar",stacked=True,
        edgecolor="black",
        color=sns.color_palette("colorblind",len(peptide_totals)),
        ax=axs[1]
    )
    axs[1].legend(facecolor="black")
    handles, labels  = axs[1].get_legend_handles_labels()
    axs[1].legend(facecolor="black").remove()
    axs[1].set_ylabel("Normalised Area")








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
    if len(t_dataset) > len(t_dataset.columns):
        print("Warning: please ensure rows=genes and columns=samples")


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



# ========= 
# Fold-changes 
# ========= 
def histone_foldchanges(normalised,control,metadata):
    normalised_groups = normalised.T
    normalised_groups["Group"] = [metadata.loc[sample].values[0] for sample in normalised_groups.index]

    group_means = normalised_groups.groupby("Group").mean()

    ctrl_means=group_means.loc[control]
    fc=pd.DataFrame()
    fc.index=ctrl_means.index
    for sample in [sample for sample in group_means.index if sample !=control]:
        fc[sample]=group_means.loc[sample]/ctrl_means


    return fc



# ========= 
# Anova-Tukey tests 
# ========= 
def histone_AnovaTukey(normalised,meta,control,fdr_correct=True):
    normalised_groups = normalised.T
    normalised_groups["Group"] = [meta.loc[sample].values[0] for sample in normalised_groups.index]
    
    pval_ser = pd.Series(
        index=normalised.index,
    )

    for peptide in normalised.index:
        # IMPORTANT! this method does not check for r^2 = 1 as it is not using general linearm odels
        # so always check the raw scores to confirm they arent mostly imputed

        # select that peptides column and the groups
        vals=normalised_groups[[peptide,"Group"]]
        
        # extract the values for each group
        groups = [group[peptide].values for _, group in vals.groupby("Group")]

        # Run the ANOVA test 
        fstat, pval= stats.f_oneway(*groups)
        
        pval_ser.loc[peptide] = pval


    sig ,fdr_adj = fdrcorrection(
        pvals=pval_ser,
        alpha=0.05,   
    )

    stats_df = pd.DataFrame(
        {
        "p_value":pval_ser,
        "fdr_adjusted":fdr_adj
        },
        index=pval_ser.index
        )
    
    # Decide whether to run tukey tests on the raw or corrected p-values
    fdr_correct =True
    if fdr_correct:
        p="fdr_adjusted"
    else:
        p="p_value"

    cols = (normalised_groups.loc[normalised_groups["Group"] != control, "Group"].unique().tolist())
    tukey_df = pd.DataFrame(columns=cols,index=normalised.index)

    for peptide in stats_df[stats_df[p] <0.05].index:
        vals = normalised_groups[[peptide,"Group"]]
        tukey = pairwise_tukeyhsd(endog=vals[peptide],
                            groups=vals["Group"],
                            alpha=0.05
                            )
        
        # Process the results of the tukey test 
        res = tukey.summary()

        res=pd.DataFrame(
                    res[1:],
                    columns=res.data[0]  
                    )

        res=pd.DataFrame(
            data=tukey._results_table.data[1:],
            columns=tukey._results_table.data[0]
        )

        res["p-adj"] = tukey.pvalues

        # filter for the results that involve the control sample and another
        res=res[
            (res["group1"].astype(str)==str(control))|
            (res["group2"].astype(str)==str(control)) 
            ]

        # set the index to the non-control groups
        res["Group_compared"]=res[["group1","group2"]].apply(lambda row: row["group1"] if str(row["group1"])!=str(control) else row["group2"],axis=1)
        res=res.set_index("Group_compared")

        # select the p-vals from this and append them to the dataframe 
        tukey_pvals = res["p-adj"]
        
        tukey_df.loc[peptide,tukey_pvals.index] = tukey_pvals

    return stats_df, tukey_df




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



def heat_sigmap(fc,TukeyPval_df):
    # Extract the data we need
    heat_fcs = fc.applymap(np.log2)
    sigmap = TukeyPval_df[heat_fcs.columns].applymap(asterisker)

    # clean the peptide names 
    heat_fcs.index = [peptide_namecleaner(pep,rem_number=True) for pep in heat_fcs.index]
    sigmap.index = [peptide_namecleaner(pep,rem_number=True) for pep in sigmap.index]

    return heat_fcs, sigmap

