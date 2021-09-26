# VariantsQC

A use friendly tool for single sample variants quality control with `vcf` format. Note: This tool is only for VCF file
generate by [GATK 4.x](https://gatk.broadinstitute.org/hc/en-us).

## request

- please add [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) script to `bin` directory
  include:`annotate_variation.pl`, `coding_change.pl`, `convert2annovar.pl`, `table_annovar.pl`
- python3.x
- python package: pandas，feather-format, numpy ,tenserflow 2.x, tqdm
- samtools (set `samtools` as executable command)

## quickly start

0. __Separating multi-allelic variants__

~~~bash
bcftools norm -m -any -f Homo_sapiens_assembly38.fasta raw.vcf.gz -Oz --threads 6 > raw.norm.vcf.gz
~~~

- There are [several ways](https://genome.sph.umich.edu/wiki/Variant_Normalization) to do this. We
  use [Bcftools]((https://samtools.github.io/bcftools/bcftools.html)) to accomplish these step.
- Please note that if you do not properly separate out multi-allelic variants, `VariantsQC` will automatically remove
  that variant in later steps.

1. __get annotated matrix__

~~~bash
python 1_vcf2matrix.py -i raw.norm.vcf.gz -o annoWithSeq.matrix.tsv -r Homo_sapiens_assembly38.fasta -d /Path_of_humandb/ --thread 6 -j 32 --reserved
~~~

- There are several parameters you can find by `python 1_vcf2matrix.py -h`
    - -i, --vcf: you must have a variants file generate by GATK with vcf format.
    - -o, ---annofile: you can set result name or generated named `result.matrix`.
    - -r, --reference: the chromosome in the reference must start with `chr` or you change the script `getRef.sh`
      in `bin`: remove `chr` in line 8.
    - -d, --humandb: database of annovar, and have: hgxx_refGene, hgxx_rmsk, hgxx_cpgIslandExt, hgxx_genomicSuperDups
    - -t, --tmp: temp directory, [./tmp/]
    - --dp: depth for filter variants. [0]
    - --refTag: reference of variants hg19/hg38. [hg38]
    - --thread: thread of annovar annotation. [1]
    - -j, --job: multiple running to get reference, if set `-j` > 1, please install parallel. [1]
    - -n, --paranum: number of lines that get reference one job, please set `-j` > 1. Recommend `50k-1m`. [100000]
    - --reserved: Whether to reserved temp documents.

2. __filter with multi dnn__

~~~bash

~~~

