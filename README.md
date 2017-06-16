# LatNet
LatNet: Latent network-based representations for large-scale gene expression data analysis.

LatNet is a generic framework for deriving latent representations of input signals based on appropriate measures and a network structure that defines the relations between the features.

In the following is an example of such latent signals. We transform an input dataset of gene expressions in different conditions into 1) regulator activity signals or 2) gene perturbations based on a regulatory network structure that defines the relations between the genes of the input dataset.

## Quick user guide

1- Load LatNet functions
```````{r}
yourPath <- '/home/myHomeDirectory/LatNet/'
source(paste(yourPath,"LatNet.R", sep = ""))
```````

2- Load your network structure (a gene regulatory network in this example), your reference data that define the reference state and your target data on which you want to measure gene perturbations
``````
GRN = your_regulator_network
refData = your_reference_data
targetData = your_target_data
``````

3- Generating regulator activity signals from the input expression data "refData" and based on the network structure defined in the regulatory network "GRN". The output is a dataset with same experimental conditions but with only the regulators (as defined in GRN) at the level of genes. Values represent the activity of each regulator at each condition.  
```
activities = .regulatorActivity(GRN, refData)
```

4- Generating gene perturbation signals from the input target expression data "targetData" and based on the reference states of regulations defined by the regulatory network structure "GRN" and the reference data "refData". The output is a dataset with same experimental conditions but with only the perturbation signals for the target genes (as defined in GRN). Values represent the perturbation of each target gene at each condition with respect to the expected expression.

a- Estimate perturbations for one gene "MyGene"
```
perturbations_for_MyGene = .oneGenePerturbations(MyGene, GRN, refData, targetData)
```
b- Perturbations for multiple genes "AllGenes" could be estimated easily through parallel computation
```
library(parallel)
allPerturbations <- t(simplify2array(mclapply(AllGenes, .oneGenePerturbations, GRN, refData, targetData)))
rownames(allPerturbations) <- AllGenes
```
