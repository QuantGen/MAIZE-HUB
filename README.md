# MAIZE-HUB

## Data collection

**Phenotypic data:** Is released to public along with metadata, agronomic information every year since 2014. These data sets have their own year-specific data identifier (DOI) and are available in the G2F repository at [https://www.genomes2fields.org/resources](https://www.genomes2fields.org/resources/).

**Genotypic data:** Markers genotypes common for years 2014-2022 can be downloaded from [https://doi.org/10.25739/tq5e-ak26](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023)


## Data curation and environmental covariates derivation
A modular code workflow as below was used to clean phenotypic and genotypic data, and to develop environmental covariates (EC) derived from crop modeling.
This workflow can be found in this [wiki](https://github.com/QuantGen/MAIZE-HUB/wiki/Pipeline-data-curation).

## Curated phenotypic, genotypic, and EC data
The curated phenotype data set includes 78,686 records of 4,372 maize hybrids, tested over 8 years (from 2014 to 2021) and 38 locations. The final set of marker genotypes includes 4,372 hybrids and 98,026 SNPs. The EC file includes 189 ECs evaluated in 136 unique year-location combinations. The following snippet shows how to read the data into an R-session.

```
PHENO=read.csv('PHENO.csv') 
ECOV=read.csv('ECOV.csv', row.names=1)
GENO=data.table::fread('GENO.csv',sep=',',data.table=FALSE) 
rownames(GENO)=GENO[,1]
GENO=as.matrix(GENO[,-1])
```

These data sets can be found in the Figshare repository at [https://doi.org/10.6084/m9.figshare.22776806](https://figshare.com/s/5d730ac680a4f6926a4a)

## Data analysis
The following pipeline provides tools to perform multivariate analyses to characterize the data set’s genetic and environmental structure, study the association of key environmental factors with grain yield and flowering traits, and provide benchmarks using state-of-the-art genomic prediction models. 

