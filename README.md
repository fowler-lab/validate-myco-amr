# Enhancement and validation of the antibiotic resistance prediction performance of a cloud-based genetics processing platform for Mycobacteria

This repository contains the data and code necessary to reproduce the tables and figures in the below preprint that is currently under review.

> Westhead J, Baker CS, Brouard M, Colpus M, Constantinides B, Hall A, Knaggs J, Lopes Alves M, Spies R, Thai H, Surrell S, Govender K, Peto TEA, Crook DW, Omar SV, Turner R, Fowler PW
> Enhancement and validation of the antibiotic resistance prediction performance of a cloud-based genetics processing platform for Mycobacteria
> bioRxiv preprint [doi:10.1101/2024.11.08.622466](https://doi.org/10.1101/2024.11.08.622466)


## Reproduce AMR predictions from the VCF files

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

