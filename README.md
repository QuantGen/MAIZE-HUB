# MAIZE-HUB

This repository has the pipelines used in [Lopez-Cruz *et al*. (2023)](https://doi.org/10.1038/s41467-023-42687-4) for curating the [Genomes to Fields](https://www.genomes2fields.org/) data and to generate Environmental Covariates (EC) for each of the year-locations in the data set. We also provide pipelines we used to analyze this data using Bayesian models.

Funding: [NSF PGRP-Tech Grant 2035472](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2035472&HistoricalAwards=false).

If you have questions about the content of this repository, please contact Marco Lopez-Cruz (lopezcru@msu.edu) and Gustavo de los Campos (gustavoc@msu.edu).

## Menu
 1. [Requeriments and software](#soft_req)
 2. [Data source](#data_source)
 3. [Data curation and environmental covariates workflow](#data_curation)
 4. [Data analysis and benchmarks modules](#data_analysis)
   
## Final curated data set

The final curated data set can be downloaded form this Figshare repository https://doi.org/10.6084/m9.figshare.22776806.

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

<a name="soft_req"></a>
### 1) Software requeriments

Environmental Covariates derivation:
  - [APSIM Next Generation](https://apsimnextgeneration.netlify.app/) crop modeling tool (>=2021.11.3.6921)
  - [apsimx](https://cran.r-project.org/web/packages/apsimx/vignettes/apsimx.html) R-package (>=2.3.1)
  - [tidygeocoder](https://jessecambon.github.io/tidygeocoder/) R-package (>=1.0.5)
  - [SFSI](https://github.com/MarcooLopez/SFSI) R-package (>=1.3.1)

Data analysis:
  - [BGLR](https://github.com/gdlc/BGLR-R) R-package (>=1.1.0)
  - [tensorEVD](https://github.com/MarcooLopez/tensorEVD) R-package (>=0.1.1)

<a name="data_source"></a>
### 2) Data source

The phenotypic and genotypic data can be downloaded from the G2F repository
  - Phenotypic data: https://www.genomes2fields.org/resources
  - Genotypic data: https://doi.org/10.25739/tq5e-ak26

These data must be downloaded and saved in the `source` folder of the following pipeline.

<a name="data_curation"></a>
###  3) Data curation and environmental covariates workflow

The pipeline used to generate the curated data set can be downloaded from this [link](https://github.com/QuantGen/MAIZE-HUB/blob/main/data_curation_and_ecov.zip).


Once the file is uncompressed, it will generate a collection of folders, each containing a module of the workflow (see folder tree below). The number that preceeds the folder name indicate the order in which the modules need to be run. Because some modules use outputs from other modules as inputs, as a rule of thumb we recommend running the modules sequentially.

Before running any module you need to populate the `source` folder with the data downloaded from G2F (see [Data source](#data_source)). 


The overall structure of the workflow is as follows

```
├── 1_phenotypes
│   ├── code
│   │   └── 1_get_phenotype.R
│   └── output
├── 2_genotypes
│   ├── code
│   │   ├── 1_SNP_filter_codify.sub
│   │   └── 2_LD_prune_genotypes.R
│   └── output
├── 3_year_loc_summary
│   ├── code
│   │   └── 1_year_location_summary.R
│   └── output
├── 4_weather
│   ├── code
│   │   └── 1_get_met_apsim.R
│   └── output
├── 5_APSIM
│   ├── code
│   │   ├── 1_apsim_sim.R
│   │   └── 1_apsim_sim.sub
│   └── output
├── 6_ecov
│   ├── code
│   │   └── 1_get_env_cov.R
│   └── output
├── 7_GDD
│   └── code
│       └── 1_add_GDD_to_pheno.R
├── README
├── source
│   ├── Agronomic_information
│   │   | # save here the agronomic information files downloaded from G2F (see Source above), one file per year.
│   ├── Genotype
│   │    # save here the compressed vcf file with genotypes, downloaded from G2F
│   ├── Metadata
│   │   | # save here the metadata downloaded from G2F (see Source above), one file per year.
│   └── Phenotype
│       | # save here the phenotype files downloaded from G2F (see Source above), one file per year.
└── tools
    ├── APSIM_functions.R
    ├── ecov_utils.R
    ├── Functions.R
    ├── LD_prune.R
    ├── read_metadata.R
    └── read_phenotype.R
```

<a name="data_analysis"></a>
### 4) Data analysis and benchmarks modules

We used the curated data set to run vairous analysis including, principal components analysis, environmental-covariate-phenotype association analysis, variance components estimation, and assessment of prediction accuracy using Bayesian models.

The workflow used for analysis can be downloaded from the following [link](https://github.com/QuantGen/MAIZE-HUB/blob/main/analysis.zip).

The uncompressed file will generate the following folders. Unlike the data curation, the analysis modules can be run in any order, since there are no dependencies between modules of the analysis workflow.


```
pipeline_analysis
├── 1_PC_genotypes
│   ├── code
│   │   ├── 1_get_G_matrix_PC.R
│   │   └── 2_make_Figure_1.R
│   └── output
├── 2_PC_ecov
│   ├── code
│   │   ├── 1_get_WWt_matrix_PC.R
│   │   └── 2_make_Figure_2.R
│   └── output
├── 3_ecov_varcomps
│   ├── code
│   │   ├── 1_ecov_variance_components.R
│   │   └── 2_make_Figure_3.R
│   └── output
├── 4_ecov_WAS
│   ├── code
│   │   ├── 1_ecov_WAS_SMR.R
│   │   └── 2_make_Figures_4_&_S6.R
│   └── output
├── 5_drought_heat_stress
│   ├── code
│   │   ├── 1_make_Figure_5.R
│   │   ├── 2_make_Figure_6.R
│   │   ├── 3_HI30_SDR_mixed_model.R
│   │   └── 4_make_Figures_7_&_S7.R
│   └── output
├── 6_genomic_prediction
│   ├── code
│   │   ├── 1_get_model_components.R
│   │   ├── 2_ANOVA_mixed_model.R
│   │   ├── 3_make_Tables_3_&_S4.R
│   │   ├── 4_10F_CV_model.R
│   │   ├── 5_LYO_CV_model.R
│   │   └── 6_make_Figures_8_&_S8.R
│   └── output
├── README
└── tools
    ├── c_utils.c
    ├── c_utils.o
    ├── c_utils.so
    ├── ecov_utils.R
    ├── Functions.R
    ├── get_folds.R
    ├── get_variance.R
    └── TENSOR_EVD.R

```

