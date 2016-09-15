#**DHS-Mining**#
dhs-mining gives epigenetic information for SNPs-eQTL pairs from [CoDeS3D](https://git.com/alcamerone/codes3d).

The dhs.db (created by init_dhsDB.py) includes the following tables from the 
[Regulatory Elements DB](http://dnase.genome.duke.edu)(Sheffield et al, 2013):
'''  
+ Gene Correlations (of known genes with DHSs, p values <= 0.05)  
+ Overlap (of DHS, CpG sites and promoters)  
+ Concordance (cell lines, types and tissues of the ENCODE samples in the database)  
+ Open Cell Types (DHS clusters with the samples, cell types and tissues in which they are found)   
+ DHS Predictors (coordinates of DHSs, their clusters and tissues)  
+ Open Samples (DHS clusters; counts of DHS, open samples; samples, tissues with max counts)  
+ Motif (assigned to transcription factors based on the JASPAR database)  
+ DHS (coordinates of DHSs and their raw signals)  
+ Gene Expression (data of the 112 samples used in the study)   
'''