#**DHS-Mining**#
dhs-mining gives epigenetic information for SNPs-eQTL pairs from [CoDeS3D](https://git.com/alcamerone/codes3d).

The dhs.db (created by init_dhsDB.py) includes the following tables from the 
[Regulatory Elements DB](http://dnase.genome.duke.edu)(Sheffield et al, 2013):
1. Gene Correlations (of known genes with DHSs, p values <= 0.05)
2. Overlap (of DHS, CpG sites and promoters)
3. Concordance (cell lines, types and tissues of the ENCODE samples in the database)
4. Open Cell Types (DHS clusters with the samples, cell types and tissues in which they are found) 
5. DHS Predictors (coordinates of DHSs, their clusters and tissues)
6. Open Samples (DHS clusters; counts of DHS, open samples; samples, tissues with max counts)
7. Motif (assigned to transcription factors based on the JASPAR database)
8. DHS (coordinates of DHSs and their raw signals)
9. Gene Expression (data of the 112 samples used in the study) 