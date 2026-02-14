# Diphenyleneiodonium Induces Global changes to Histone regulation in MYCN-amplified Neuroblastoma
If you are reading this, thank you for visiting my GitHub repo :). 

These are the scripts I used to analyse data for my undergrad thesis.

I analysed histone mass spectrometry, total proteomics and RNA-seq data. I also queried the ChEMBL database


### Histone MS
This analysis pipeline consisted of pre-processing, post-processing and analysis steps. 
- Preprocessing was carried out using Epiprofile 2.2 to convert raw mass spec files to Area (Yuan et al, 2018)
- Post-processing largely followed the recomendations laid out by Thomas et al. 2020. 
- Analysis was performed mostly in Python.


### RNA-seq & Total Proteomics
- Analysis was performed in python with PyDeseq2 for differential expression, the differentially expressed genes were analysed using reactome pathway analysis.

### ChEMBL analysis 
- I also queried the ChEMBL database using MySQL to extract a table of pChEMBL values
- This is included in Appendix 1. 




<img width="6000" height="4200" alt="Copy of Thesis_Figure" src="https://github.com/user-attachments/assets/e80106e5-f6c8-466c-b63c-d9e599b3642c" />



### **References**
Yuan, Z.F., Sidoli, S., Marchione, D.M., Simithy, J., Janssen, K.A., Szurgot, M.R. and Garcia, B.A., 2018. EpiProfile 2.0: a computational platform for processing epi-proteomics mass spectrometry data. Journal of proteome research, 17(7), pp.2533-2541.

Thomas, S.P., Haws, S.A., Borth, L.E. and Denu, J.M., 2020. A practical guide for analysis of histone post-translational modifications by mass spectrometry: best practices and pitfalls. Methods, 184, pp.53-60.

Muzellec, B., Teleńczuk, M., Cabeli, V. & Andreux, M. 2023. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39, btad547.
