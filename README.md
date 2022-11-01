Dockerized DNA-seq-pipeline
==================================

This Docker version snakemake pipeline implements the DNA-seq pipeline for calling SNV/CNV/SV variants.

## Authors

* Jingyu Yang

## Usage


#### Step 1: clone the docker version DNA-seq-pipeline.
Clone the project-repository:

```
git clone [link]
```

#### Step 2: Download reference data for DNA-seq-pipeline:
1. Download following GDC genome reference data into your `[ref file path]` directory.
   
   1. [Reference genome(fastq)](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834).
   2. [BWA Index](https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225).
   3. [GATK Index](https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5)

2. Download following GATK reference data into your `[ref file path]` directory.

   1. know_dbsnp_vcf ([Homo_sapiens_assembly38.dbsnp138.vcf](https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf), [Homo_sapiens_assembly38.dbsnp138.vcf.idx](https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx)) 
   2. interval_list ([wgs_calling_regions.hg38.interval_list](https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list))
   3. gatk_germline_resource ([af-only-gnomad.hg38.vcf.gz](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz), [af-only-gnomad.hg38.vcf.gz.tbi](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi)), gatk_panel_of_normal([1000g_pon.hg38.vcf.gz](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz),[1000g_pon.hg38.vcf.gz.tbi](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi)) 
   4. exac_common_knownsite ([small_exac_common_3.hg38.vcf.gz](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz), [small_exac_common_3.hg38.vcf.gz.tbi](https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi))


#### Step 3: Download control bam data for DNA-seq-pipeline:
Store the corresponding control bam and index file as well as the reference genes interval bed file in your `[control bam path]` directory. You can directly use our default gene intervel bed `genesHG38.bed` been provided.
If there is no corresponding control bam, you can download the public control cram file [NA12878](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram) into your `[control bam path]` directory.
Modify the 'CramToBam.sh' script with your `[control bam path]` and `[ref file path]` and run it to generate the control bam file.
Start terminal and run the command below:
	
	sh CramToBam.sh

#### Step 4: Prepare the fastq raw data
Make sure there only one sample's fastq data located in your `[raw fastq data path]` directory.
    
    Eg. sample.R1.fq, sample.R2.fq


#### Step 5: Build the docker image.
1. Switch to the docker DNA-seq-pipeline project-folder:

```
cd ./dna_seq_pipeline
```

2. Build the Docker Image (named onkopus):

```
docker build -t snakemake_local .

```

#### Step 5: Starting and Usage of the Application:
1. Check if the snakemake available on docker container(set the reference path,fastq input path and pipeline path in advance).

```
docker run -v [ref file path]:/data/ref/ \
-v [raw fastq data path]:/data/InputFastqDir/ \
-v [control bam path]:/data/control_bam/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -v

```

2. Dry run of the snakemake DNA-Seq-Pipeline.(check the processing steps)

```
docker run -v [ref file path]:/data/ref/ \
-v [raw fastq data path]:/data/InputFastqDir/ \
-v [control bam path]:/data/control_bam/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -j all --use-conda -n

```
3. Run pipeline and get the bam and vcf files in local.

```
docker run -v [ref file path]:/data/ref/ \
-v [raw fastq data path]:/data/InputFastqDir/ \
-v [control bam path]:/data/control_bam/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -j all --use-conda

```

For example:

```
 docker run -v /home/kingyang/harddisk/DNA_Seq_pipeline/dna_seq_pipeline_snakemakeconda/Ref_data/:/data/ref/ \
-v /home/kingyang/harddisk/DNA_Seq_pipeline/dna_seq_pipeline_snakemakeconda/control_bam/:/data/control_bam/ \
-v /home/kingyang/harddisk/SSH_server_Data/Fastq_test:/data/InputFastqDir/ \
-v $(pwd):/work snakemake_local snakemake -j all --use-conda

```
4. Modify the Snakemake file to generate the specific results in need. 
