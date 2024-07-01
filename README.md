# Analysis code for "A compendium of Amplification-Related Gain Of Sensitivity (ARGOS) genes in human cancer"

> Chromosomal gains are among the most frequent somatic genetic alterations
> occurring in cancer. While the effect of sustained oncogene expression has
> been characterized, the impact of copy-number gains affecting
> collaterally-amplified “bystander” genes on cellular fitness remains less
> understood. To investigate this, we built a comprehensive map of dosage
> compensations across human cancers by integrating expression and copy number
> profiles from over 8,000 TCGA tumors and CCLE cell lines. Further, we
> analyzed the effect of gene overexpression across 17 human cancer ORF screens
> to provide an overview of genes that prove toxic to cancer cells when
> overexpressed. Combining these two independent approaches we propose a class
> of ‘Amplification-Related Gain Of Sensitivity’ (ARGOS) genes. These genes are
> located in commonly amplified regions of the genome, have lower expression
> levels than expected by their copy-number status, and are toxic to cancer
> cells when overexpressed. We experimentally validated CDKN1A and RBM14 as
> high-confidence pan-cancer ARGOS genes in lung and breast cancer cell line
> models. We additionally suggest that RBM14’s mechanism of toxicity involves
> altered DNA damage response and innate immune signaling processes following
> gene overexpression. Finally, we provide a comprehensive catalog of
> compensated, toxic, and ARGOS genes as a community resource.

## Analysis directories

The main Analysis steps are performed in the following directories:

* `data` - Creating the analytical data sets
* `model_compensation` - The compensation analysis
* `model_orf` - The ORF toxicity analysis
* `report` - The scripts to create all figures in the manuscript

Those can be run via the `snakemake` workflow manager once all dependencies are
met.

## Dependencies

* All required R packages and their versions are listed in `renv.lock`.
* The `ebits` toolkit and setup according to the project README
* The ORF and TCGA data sets from their respective publications
