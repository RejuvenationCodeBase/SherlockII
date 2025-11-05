# Authors:
Chris Fuller, Jialiang Gu, Jiashun Zheng @ UCSF

# Sherlock-II: calculating gene-phenotype p-values using GWAS and eQTL data

# Python environment setup

Sherlock_II runs on Python3 using a few very common Python modules including: Numpy, Scipy, pytabix, and Matplotlib. 

## python 3 environment with pip

```
python3 -m venv sherlock_env
source sherlock_env/bin/activate
pip install -r venv_python3_requirements.txt
```

# Usage example:
**IMPORTANT! DOWNLOAD the sample data first (see below)**

```python sherlock_II.py Schizophrenia_PGC_SCZ2 Merged_041918 sherlock_II.cfg```

1st argument: GWAS name (Schizophrenia_PGC_SCZ2)

2nd argument: eQTL name (Merged_041918)

3rd argument: a configuration file containing data_folder_location, output_folder_location and 3 numeric parameters 

**IMPORTANT! before running the script, make sure you edit the the configuration (sherlock_II.cfg in the example) to point to the correct directory**

## Data folder structure and sample data download
Data folder should contain a folder name "eQTL" and a folder name "GWAS", please use the link below to download the data and put it in the appropriate location.

### Downloading eQTL data
```wget http://genome.ucsf.edu/~jiashun/genome_browser/SherlockII/eQTL.tar.gz ```

Size: 2GB, md5: 0427fbf717076030db678f8469ac0f84, unzip the file in the data folder. This is the eQTL data we used in the manuscript.  

### Downloading GWAS data
```wget http://genome.ucsf.edu/~jiashun/genome_browser/SherlockII/GWAS.tar.gz```

Size: 781MB, md5: 470f713e2d552137049be762cc14d4e1, unzip the file in the data folder.

For reviewing purpose, we are providing the following GWAS in the GWAS.tar.gz, which includes the summary statistics for the three pairs: CD and RA,  Schizophrenia and Height, and Alzheimer’s disease and Breast cancer (mentioned in the manuscript). We also provided a GWAS dataset on Autism. Reviewers can also run the GWAS data they choose by generating the data in the same form and placing them in the same directory  (see detailed instruction on the format and preprocessing below).


We will provide download for the full set of summary statistics data we used upon publication of the paper.  

### Download the GWAS catalog data for gene annotation:

In the data folder download the file:

```wget http://genome.ucsf.edu/~jiashun/genome_browser/SherlockII/gwas_catalog_v1.0-associations_e98_r2019-12-16.tsv```

## The data directory structure should looks something like this
```
data/
├── eQTL
│   ├── allSNPs.LDDB
│   ├── Merged_041918.gene
│   ├── Merged_041918.info
│   ├── Merged_041918.LD02
│   ├── Merged_041918.LD02.alias
│   ├── Merged_041918.LD02.gz
│   ├── Merged_041918.LD02.gz.tbi
│   ├── Merged_041918.LD85
│   ├── Merged_041918.LD85.alias
│   ├── Merged_041918.LD85.gz
│   ├── Merged_041918.LD85.gz.tbi
│   ├── Merged_041918.meta
│   ├── Merged_041918.pickle
│   └── Merged_041918.txt
├── GWAS
│   ├── Alzheimer_Neuropathologic_Neuritic_np_status_relax.meta
│   ├── Alzheimer_Neuropathologic_Neuritic_np_status_relax.tab.gz
│   ├── Alzheimer_Neuropathologic_Neuritic_np_status_relax.tab.gz.tbi
│   ├── Alzheimer_Neuropathologic_Neuritic_np_status_relax.txt
│   ├── BreastCancer_pha001338.meta
│   ├── BreastCancer_pha001338.tab.gz
│   ├── BreastCancer_pha001338.tab.gz.tbi
│   ├── BreastCancer_pha001338.txt
│   ├── Franke_Cronhs_meta.meta
│   ├── Franke_Cronhs_meta.tab.gz
│   ├── Franke_Cronhs_meta.tab.gz.tbi
│   ├── Franke_Cronhs_meta.txt
│   ├── GIANT_HEIGHT_Wood_et_al_2014.meta
│   ├── GIANT_HEIGHT_Wood_et_al_2014.tab.gz
│   ├── GIANT_HEIGHT_Wood_et_al_2014.tab.gz.tbi
│   ├── GIANT_HEIGHT_Wood_et_al_2014.txt
│   ├── RA_GWASmeta_TransEthnic.meta
│   ├── RA_GWASmeta_TransEthnic.tab.gz
│   ├── RA_GWASmeta_TransEthnic.tab.gz.tbi
│   ├── RA_GWASmeta_TransEthnic.txt
│   ├── Schizophrenia_PGC_SCZ2.meta
│   ├── Schizophrenia_PGC_SCZ2.tab.gz
│   ├── Schizophrenia_PGC_SCZ2.tab.gz.tbi
│   └── Schizophrenia_PGC_SCZ2.txt
└── gwas_catalog_v1.0-associations_e98_r2019-12-16.tsv
```

## Expected output
We expect to see a pickle file and a directory in the output_folder/html/ folder. The html folder contain human readable results and the pickle store the data structure used by the program. You can ignore the pickle file as it exists mainly to store intermediate results when we were developing Sherlock-II. 

## Recommend computation hardware requirement
Sherlock-II take quite amount of memory when performing alignment of SNPs between eQTL and GWAS and tagging (finding independent SNP). A computer with at least 32G of memory is recommend to run the tool smoothly. Typically the run time varies from a few minutes to an hour depends on the size of the GWAS and eQTL data (and the speed of CPU). 

We also provided a lite version of SherlockII which take the pre-alignment and tagged file for computation. Check the ../SherlockII_lite/ folder.

## Brief explaination on GWAS data format

For GWAS data, we expect to see four files: 
1. GWAS_name.meta file store the meta information of the GWAS, we recommend to take one of the example data as an example and change the GWAS name in the first three lines. 
2. GWAS_name.txt, this is converted from the GWAS summary statistics, the columns are:
rsid<tab>pvalue<tab>chromosome<tab>position<tab>maf(minor allele frequency)
The chromosome and position are taken from dbsnp with hg19 annotation.
3. GWAS_name.tab.gz, reorder the GWAS_name.txt to have columns:
chrom<tab>position<tab>rsid<tab>pval<tab>maf
Then sorted by chromosome and position, then use tabix to index it (http://www.htslib.org/doc/tabix.html)
```
bgzip GWAS_name.tab
tabix -b 2 -e 2 GWAS_name.tab.gz
```
4. GWAS_name.tab.gz.tbi, tabix index (generated from tabix -b 2 -e 2 GWAS_name.tab.gz) 

## eQTL data format:
```
├── eQTL
│   ├── allSNPs.LDDB : list of all SNPs we have LD data
│   ├── Merged_041918.gene : All genes in the eQTL
│   ├── Merged_041918.info : gene<tab>rsid<tab>pval<tab>cis/trans(1/2)<tab>chrom<tab>position<tab>maf<tab>source tissue
│   ├── Merged_041918.LD02 : eQTL SNPs with LD>=0.2 to any other eQTL SNPs. 
│   ├── Merged_041918.LD02.alias : alias of SNPs in LD02 (same SNP with different names)
│   ├── Merged_041918.LD02.gz : bgzip result of LD02 file
│   ├── Merged_041918.LD02.gz.tbi: tabix index
│   ├── Merged_041918.LD85 : SNPs with LD>=0.85 to any eQTL SNPs
│   ├── Merged_041918.LD85.alias: alias of SNPs in LD85
│   ├── Merged_041918.LD85.gz : bgzip of the LD85 file
│   ├── Merged_041918.LD85.gz.tbi: tabix index
│   ├── Merged_041918.meta : meta data
│   ├── Merged_041918.pickle : pickle file (preloaded eQTL data into SherlockII format)
│   └── Merged_041918.txt : gene<tab>rsid<tab>pval<tab>cis/trans
```

The purpose of all these files is to speed up the searching through SNPs in the alignment and tagging steps. Please contact us if you need help generating input for your own dataset. 
