import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import re
import numpy as np
import statistics

from scipy.stats import f_oneway

imr_lfq = pd.read_csv(r"c:\Users\Tim\OneDrive\__ThesisProject\_DryLab\Datasets\IMR5-75_totalprot\matrix.tsv", delimiter="\t")
sample_data = pd.read_excel("C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Datasets\\IMR5-75_totalprot\\Sample_descriptions.xlsx")
sample_groups = dict(
    zip(
        sample_data["Sample ID"],
        sample_data["Group name"].replace("DOX-","Control").values
        )
    )
# delete the sample metadata df by re-assigning the variable
sample_data = imr_lfq.set_index("Genes")[sample_groups.keys()]

# Create our color pallete for our samples (we keep this constant for all figures)
pal =["#3b3b3b"] + [mcolors.to_hex(c) for c in sns.color_palette("Blues",3)] + ["Red","purple"]
pal_r = ["#3b3b3b"] + ["#2070b4"] + ["Red","purple"]

def get_iqr_houskeep(gene,housekeep_data):
    
    # Calculate Interquartile Range (IQR) using numpy
    q1 = np.percentile(housekeep_data[gene], 25)
    q3 = np.percentile(housekeep_data[gene], 75)
    iqr = q3 - q1

    # Calculate Quartile Deviation
    quartile_deviation = (q3 - q1) / 2

    print(gene)
    print("Interquartile Range (IQR):", iqr)
    print("Quartile Deviation:", quartile_deviation)
    return iqr




# ===============================
# ===============================

# This code was used to generate the figure validating mycn knockdown
# Create the treatment column
mycn_data = pd.DataFrame(sample_data.loc["MYCN"])
mycn_data["Treatment"] = [sample_groups[sample] for sample in mycn_data.index]
# As we are interested in the % knockdown, we take the average of the control
mean_mycn = mycn_data["MYCN"][mycn_data["Treatment"]=="Control"].mean()
mycn_data["% of Control"] = (mycn_data["MYCN"]/mean_mycn)*100

plt.figure(figsize=(10,6))

sns.barplot(mycn_data,
            x="Treatment",
            y="% of Control",
            palette=pal,  
            capsize=0.05,
        err_kws={
            "color": "black",    
            "linewidth": 1.2, 
        },
        edgecolor="black"
            )

plt.show()










# ===============================
# ===============================

sample_data_restricted = sample_data[[col for col in sample_data.columns if sample_groups[col] in ["Control","DOX+ 72h","DOX-DPI+","DOX+DPI+"]]]

sns.set_style("white")
housekeepers = ["ACTB","GAPDH","HSPD1","VCL"]
fig, axs =plt.subplots(ncols=2,nrows=2,figsize=(18,12))
lim_dif=80000

for gene,ax in zip(housekeepers,axs.flatten()):
    lfq_vals = sample_data_restricted.loc[gene]
    if isinstance(lfq_vals,pd.DataFrame):
        # GAPDH has two rows...possibly different isoforms?
        # so we sum them as presumably antibodies bind the same epitope on both forms
        lfq_vals = lfq_vals.sum()
        lfq_vals.name = gene

    housekeep_data = pd.DataFrame(lfq_vals) 
    housekeep_data["Treatment"] = [sample_groups[sample] for sample in housekeep_data.index]

    # Run an ANOVA test
    group_vals = []
    group_cvs = []
    print(gene,"cv values")
    for treatment in set(housekeep_data["Treatment"]):
        group_data = list(housekeep_data[housekeep_data["Treatment"]==treatment][gene])
        group_cvs.append(statistics.stdev(group_data)/statistics.mean(group_data))
        print(statistics.stdev(group_data)/statistics.mean(group_data))

        group_vals.append(group_data)
    fstat,pval= f_oneway(*group_vals)
    mean_cv = statistics.mean(group_cvs)


    plt.figure(figsize=(10,6))

    sns.boxplot(housekeep_data,
                x="Treatment",
                y=gene,
                palette=pal_r,  
                ax = ax
                    )
    
    sns.stripplot(housekeep_data,
                x="Treatment",
                y=gene,
                ax = ax,
                color="black"
                    )
                    
    iqr = round(get_iqr_houskeep(gene,housekeep_data),1)
    ax.legend(
            title=f"ANOVA_pval: {round(pval,4)} \nIQR: {iqr} \nmean_CV: {round(mean_cv,3)}",
            facecolor="lightgrey",
            loc="lower left"

    )


    ax.set_ylim(housekeep_data[gene].min()-lim_dif, housekeep_data[gene].max()+lim_dif)
    ax.set_title(gene,fontweight="bold")
    ax.set_xlabel(None)
    ax.set_ylabel("LFQ")


fig.savefig(r"C:\Users\Tim\OneDrive\__ThesisProject\Figures\Appendix2_loadingc\LoadingPlot_boxed_newest.svg",dpi=600,transparent=False)


