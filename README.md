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

Once the files are downloaded and uncompressed, it will generate the following files


```
data
├── ECOV.csv
├── ECOV_KL.csv
├── ECOV_layered.csv
├── ECOV_period.csv
├── GENO.csv
├── MAP.csv
└── PHENO.csv
```

You can read them into an R environment using the following code.

```r
 PHENO=read.csv('PHENO.csv') 
 ECOV=read.csv('ECOV.csv', row.names=1)
 GENO=data.table::fread('GENO.csv',sep=',',data.table=FALSE) 
 rownames(GENO)=GENO[,1]
 GENO=as.matrix(GENO[,-1])
```

###  Data curation and environmental covariates derivation

The pipeline used to generate the curated data set can be downloaded from this (link](https://github.com/QuantGen/MAIZE-HUB/blob/main/pipeline_analysis.zip).

Once the file is uncompressed, it will generate a collection of folders, each containing a module of the workflow. The number that preceeds the folder name indicate the order in which the modules need to be run. For instance the folder `2_`

The overall structure of the workflow is as follows

```


```

## Benmchmarks

The following pipeline provides tools to perform multivariate analyses to characterize the data set’s genetic and environmental structure, study the association of key environmental factors with grain yield and flowering traits, and provide benchmarks using state-of-the-art genomic prediction models. 

