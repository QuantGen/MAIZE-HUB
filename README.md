# MAIZE-HUB

This repository has the pipelies used for curating [Genomes to Fields](https://www.genomes2fields.org/) data and to generate Environmental Covariates (EC) for each of the year-locations in the data set. We also provide pipelines we used to analyze this data using Bayesian models.

If you have questions about the content of this repository, please contact Marco Lopez-Cruz (lopezcru@msu.edu) and Gustavo de los Campos (gustavoc@msu.edu).
   

## Data source

The phenotypic and genotypic data can be downloaded from the G2F repoisotry
  - **[Phenotypic data](https://www.genomes2fields.org/resources/)**
  - **[Genotypic data](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023)**

The pipelines use the files in the above links as inputs.


### Final curated data set

The final curated data set can be donwloaded form this Figshare [DOI](https://figshare.com/s/5d730ac680a4f6926a4a).

Once the files are downloaded and uncompressed, you can read them into an R environment using the following code.

```r
 PHENO=read.csv('PHENO.csv') 
 ECOV=read.csv('ECOV.csv', row.names=1)
 GENO=data.table::fread('GENO.csv',sep=',',data.table=FALSE) 
 rownames(GENO)=GENO[,1]
 GENO=as.matrix(GENO[,-1])
```

###  Data curation and environmental covariates derivation

The folloiwng [wiki page](https://github.com/QuantGen/MAIZE-HUB/wiki/Pipeline-data-curatio) shows the workflow we use to curate geneotypes and phenotypes and to genrate the enviromental covariates. 

## Benmchmarks

The following pipeline provides tools to perform multivariate analyses to characterize the data set’s genetic and environmental structure, study the association of key environmental factors with grain yield and flowering traits, and provide benchmarks using state-of-the-art genomic prediction models. 

