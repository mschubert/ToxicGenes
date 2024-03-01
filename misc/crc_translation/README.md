Read Me
=======

Files:

Supplementary_Table_01.xlsx: Clinical data file including survival information (OS and RFS)

CRC-SW.FACETS.facets-suite.1063_DNBSEQ.20210706.gene_level.gencode_v35.ProteinCoding_IG_TR.txt: copy number variation file by gene level. File columns include:

* sample (same as DNA_Tumor_Sample_Barcode in the clinical file)
* gene (gene name)
* chrom (chromossome)
* gene_start
* gene_end
* strand (gene strand)
* ensg_id
* hgnc_id
* gene_type
* source (of gene information)
* tsg (tumour supressor gene TRUE or FALSE)
* oncogene (TRUE or FALSE)
* seg (CNV segment ID)
* median_cnlr_seg (Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal)
* segclust (Segment cluster to which the segment belongs)
* seg_start
* seg_end
* cf.em (Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid)
* tcn.em (Total copy number. 2 for normal diploid)
* lcn.em (Lesser (minor) copy number. 1 for normal diploid)
* cf (Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid)
* tcn (Total copy number. 2 for normal diploid)
* lcn (Lesser (minor) copy number. 1 for normal diploid)
* seg_length
* mcn (not sure, I think is minor copy number same as lcn)
* genes_on_seg (Number of genes in the segment)
* gene_snps (Number of SNPs in the segment)
* gene_het_snps (Number of SNPs that are deemed heterozygous)
* spans_segs (If gene spans more than one segment TRUE or FALSE)
* cn_state (Final copy number state according to tcn and lcn: deletion (tcn=0), LOH (tcn=1), cnLOH (tcn=2 and lcn=0), normal (tcn=2 and lcn=1), gainLOH (tcn=3 or 4 and lcn=0), ampLOH (tcn>=5 and lcn=0), gain (tcn=3 or 4 and lcn>=1) and amplification (tcn>=5 and lcn>=1)).
* filter

There is also data uploaded to EVA
(https://www.ebi.ac.uk/eva/?eva-study=PRJEB61514) for CNV, SV and small
mutations (coding and non coding) per patient sample. In that CNV files at EVA
they are not gene annotated so that’s why I uploaded this one here instead
since it’s easier to use.
