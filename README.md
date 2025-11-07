# Enhancement and validation of the antibiotic resistance prediction performance of a cloud-based genetics processing platform for Mycobacteria

This repository contains the data and code necessary to reproduce the tables and figures in the below preprint that is currently under review.

> Westhead J, Baker CS, Brouard M, Colpus M, Constantinides B, Hall A, Knaggs J, Lopes Alves M, Spies R, Thai H, Surrell S, Govender K, Peto TEA, Crook DW, Omar SV, Turner R, Fowler PW
> Characterising the performance of an antibiotic resistance prediction tool, gnomonicus, using a diverse testset of 2,663 Mycobacterium tuberculosis samples.
> bioRxiv preprint [doi:10.1101/2024.11.08.622466](https://doi.org/10.1101/2024.11.08.622466)


## Reproduce AMR predictions from the VCF files using `gnomonicus`

First install `gnomonicus`; the easiest way to do this is using `pip`

```
pip install gnomonicus
```

This should automatically place `gnomonicus` in your `$PATH`. Now you need to get the WHOv2 catalogue and H37Rv version 3 Genbank file via

```
git clone git@github.com:fowler-lab/validate-myco-amr.git
cd ..
git clone git@github.com:oxfordmmm/tuberculosis_amr_catalogues.git
cd validate-myco-amr/
```

Lastly you'll need to have installed [GNU Parallel](https://www.gnu.org/software/parallel/) for the below to work (this nicely uses all the cores on your machine to speed up the processing). On a Mac this is easiest via MacPorts or Brew. The below should take 1-2 hours on a Mac laptop with an M-series CPU.

```
cd dat/outputs/ukmyc/
find dat/ -name '*vcf' | parallel --bar gnomonicus --vcf_file {} --catalogue_file ../../../../tuberculosis_amr_catalogues/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv --json --genome_object ../../../../tuberculosis_amr_catalogues/catalogues/NC_000962.3/NC_000962.3.gbk --min_dp 3
cd ../mgit/
find dat/ -name '*vcf' | parallel --bar gnomonicus --vcf_file {} --catalogue_file ../../../../tuberculosis_amr_catalogues/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv --json --genome_object ../../../../tuberculosis_amr_catalogues/catalogues/NC_000962.3/NC_000962.3.gbk --min_dp 3
```

The above steps will create JSON files in the same directories as the VCF files. These JSON files contain the AMR predictions that will be used in the downstream analysis.

## Parse the `gnomonicus` output JSON files and create the results tables

For simplicity there is a Python script that will detect the output JSON files and recreate `dat/RAW_EFFECTS.csv` and `dat/RAW_PREDICTIONS.csv`. To recreate these tables issue

```
python bin/parse_gnomonicus.py
```

## Reproduce AMR predictions from TB-Profiler

The output JSON files from TB-Profiler are included in the repository in `dat/tbprofiler`. If you wish to reproduce these you will need to download the FASTQ files for all 2,663 samples using these scripts

```
bash bin/UKMYC_1000_samples_download.sh
bash bin/MGIT_1663_samples_download.sh
```

Then you can either use the web portal at https://tbdr.lshtm.ac.uk/ or install TB-Profiler locally yourself, or contact the TB-Profiler team for assistance. Whichever way you choose, you'll need to download and keep all the output JSON files (the ones used in the analysis are stored in `dat/outputs/tbprofiler/`).

## Parse the TB-Profiler output JSON files and recreate the results tables 

The script below will parse the TB-Profiler output JSON files and recreate `dat/tbprofiler_EFFECTS.csv` and `dat/tbprofiler_PREDICTIONS.csv`.

```
python bin/parse_tbprofiler.py
```

## Create the results tables used to drive the analysis

The first two jupyter notebooks (`01-create-phenotypes-table.ipynb` and `02-create-results-tables.ipynb`) are designed to step you through (i) creating the one row per sample per pDST phenotype table (`dat/PHENOTYPES.csv`, includes Table 1) and (ii) reading in the raw `EFFECTS` and `PREDICTIONS` tables from both `gnomonicus` and `TB-Profiler` and creating e.g. the `RESULTS`  table.

The last two notebooks are designed to then step you through re-creating all the figures and tables in the manuscript. The first notebook (`03-main-analysis.ipynb`) creates the main results figures and tables, and the second notebook (`04-discrepancy-analysis.ipynb`) creates the discrepancy analysis tables.

