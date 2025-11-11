# Effects of Diphenyleneiodonium Chloride on the Metabolic Epigenetic Axis in MYCN-Amplified Neuroblastoma
If you are reading this, thank you for visiting my GitHub repo :)

Found here are the scripts I used to analyse data for my thesis.

I analysed histone mass spectrometry, total proteomics and RNA-seq data. 

### Histone MS
This analysis pipeline consisted of pre-processing, post-processing and analysis steps. 

- Preprocessing was carried out using Epiprofile 2.2 to convert raw mass spec files to Area (Yuan et al, 2018)
- Post-processing largely followed the recomendations laid out by Thomas et al. 2020 and was performed in R. 
- Analysis was performed mostly in Python.


### RNA-seq
- Analysis was performed in python with PyDeseq2 for differential expression.

### **References**
Yuan, Z.F., Sidoli, S., Marchione, D.M., Simithy, J., Janssen, K.A., Szurgot, M.R. and Garcia, B.A., 2018. EpiProfile 2.0: a computational platform for processing epi-proteomics mass spectrometry data. Journal of proteome research, 17(7), pp.2533-2541.

Thomas, S.P., Haws, S.A., Borth, L.E. and Denu, J.M., 2020. A practical guide for analysis of histone post-translational modifications by mass spectrometry: best practices and pitfalls. Methods, 184, pp.53-60.

Muzellec, B., Teleńczuk, M., Cabeli, V. & Andreux, M. 2023. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39, btad547.
