# Genomic Passport
In this repository one can find the code that was written for the genomic passport project carried out at the [Fellay Lab](http://fellay-lab.epfl.ch/) at [EPFL](www.epfl.ch).

## List of files in the repository
### Analysis_Walkthrough.ipynb
This notebook contains a walkthrough of the analysis pipeline that was implemented in the main.py files, to explain the steps that were taken to obtain the appropriate results.
### Data
This files contains the pharmacogenomic recommendation data that was used for the project.
1. **total_var_ann.csv** contains drug recommendations over all variants, takent from the [PharmGKB](https://www.pharmgkb.org/) database.
2. **HLA_var_ann.csv** contains drug recommendations from the same database but only for HLA alleles.

### drugSNP.txt
This files contains a list of all relevant SNP's that are in the drug recommendation database, as to be able to extract them using [Plink](http://pngu.mgh.harvard.edu/~purcell/plink/) software, and to speed up the subsequent analysis with Python.
### European-HLA4-hg19.RData & HLA_from_patient.R
These two files are necessary for the R script needed to perform the [HIBAG](https://www.ncbi.nlm.nih.gov/pubmed/23712092)  method for HLA imputation. It relies on a classifier designed by a group at the [University of Washington](http://www.biostat.washington.edu/~bsweir/HIBAG/). The code also relies on a [Bioconductor HIBAG](http://www.biostat.washington.edu/~bsweir/HIBAG/) library. 
### main
Automatized analysis pipeline. It was codod with [Python 3.5.2](https://www.python.org/downloads/release/python-352/)
##### Dependancies
 - [Pandas](http://pandas.pydata.org/)
 - R and Bioconductor as cited above
 - [Plink](http://pngu.mgh.harvard.edu/~purcell/plink/) version 1.9
 
Input : patientName.vcf , 23&Me genotype data

How to run the program (From the Terminal)
`python main.py -f patientName.vcf`
 
 
### Results
Folder where the results of the pipeline are stored. One folder per day with one folder per patient.