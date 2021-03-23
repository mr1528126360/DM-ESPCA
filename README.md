# DM-ESPCA

#Abstract

Motivation: The latest research shows that each type of cancer can be divided into multiple sub-types, which is one of the key reasons that make cancer difficult to cure. 
Therefore, finding a new target gene of cancer subtypes has great significance on developing new tumor drugs and person-alized treatment plan. Because the genetic data of cancer 
are usually high-dimensional and have multiple potential subtype classifications, the sparse principal component analysis (PCA) model has become a natural candidate method for 
identifying cancer subtype targets and subtype cluster-ing. Researchers has been proposed some cancer subtype gene probe screening models based on sparse PCA. However, the 
performance of existing models is not good enough. We have no-ticed that none of the existing models using the existing classification information of cancer sub-types as prior 
information and the existing models also cannot solve the sample quality problem.

Results: In order to solve the problems above, we propose Dynamic Metadata Edge-group Sparse PCA (DM-ESPCA) model. DM-ESPCA model first uses the known classification information 
of cancer subtypes as prior knowledge. Next, this model filters the metadata set to build a dynamic gene network for each cancer subtype. Finally, the dynamic network is used to 
establish a sparse PCA model to screen the potential targets corresponding to each subtype. The experiments in this paper show that the DM-ESPCA model is better than the 
existing sparse PCA model. Compared with the existing models, the re-clustering results of the DM-ESPCA model can be improved by up to 23% of accuracy. The results of machine 
learning classification models can be improved by up to 22% of accuracy. We believe that the DM-ESPCA model is an important supplement to the can-cer subtype target screening 
model. The DM-ESPCA model can help researchers better conduct cancer subtype research.

#process

1.filter metadata for each subtype in the cancer dataset. 

2.based on metadata, use known subtype classification information as prior information, calculate the correlation degree of each gene probe corresponding to each subtype.

3.using the quantita-tive value of correlation as a parameter to generate a unique biolog-ical network for each subtype. 

4.build DM-ESPCA model using a dynamic gene network to select biomarkers for each sub-type.

#function

DM_ESPCA = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=20, err=0.01, Num.init=5, w_l)

#parameter

X: gene expression data, n*p, n is sample and p is gene.

k: the amount of PCs and PC loadings.

overlap.group: data of gene pathway. 

k.group: the amount of edges reserved.

we:

t:

niter:

err:

Num.init:

w_l: t-value data of each gene corresponding to subtype.

#output

return (list(U, D, V))

#example

There is a gastric cancer data set (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-35809/) as an example in the 'data' folder, which can be used directly after 
decompression.

And the result of the example is saved in the 'result' folder.
