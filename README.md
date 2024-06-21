Dockerized DNA-seq-pipeline
==================================

This Docker version snakemake pipeline implements the DNA-seq pipeline for calling SNV/CNV/SV variants from raw FASTQ  data.

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
   
   1. [Reference genome(FASTQ)](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834).
   2. [BWA Index](https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225).
   3. [GATK Index](https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5)

2. Download following GATK reference data into your `[ref file path]` directory.

   1. know_dbsnp_vcf ([dbsnp_146.hg38.vcf.gz](wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz), [dbsnp_146.hg38.vcf.gz.tbi](wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi))
   2. know_dbsnp_vcf ([Homo_sapiens_assembly38.known_indels.vcf.gz](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz), [Homo_sapiens_assembly38.known_indels.vcf.gz.tbi](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi))
   3. gold_standard_indels ([Mills_and_1000G_gold_standard.indels.hg38.vcf.gz](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz), [Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi))
   4. 1000G_snpshigh_confidence ([1000G_phase1.snps.high_confidence.hg38.vcf.gz](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz), [1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi))
   5. interval_list ([wgs_calling_regions.hg38.interval_list](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list))
   6. af_only_gnomad ([af-only-gnomad.hg38.vcf.gz](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz), [af-only-gnomad.hg38.vcf.gz.tbi](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi))
   7. gatk_panel_of_normal ([1000g_pon.hg38.vcf.gz](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz), [1000g_pon.hg38.vcf.gz.tbi](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi))
   8. exac_common_knownsite ([small_exac_common_3.hg38.vcf.gz](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz), [small_exac_common_3.hg38.vcf.gz.tbi](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi))
   

#### Step 3: Download bed files for CNVkit:
Download the target and antitarget bed files for CNVkit in CNV detection. Store bed files in your `[CNVkit bed path]` directory and configure the file name and path in config.yaml file (hg38gene_bed: bed_files/`target.bed`
hg38access_bed: bed_files/`antitarget.bed`). Alternatively, users can directly use the default bed files we provided.


#### Step 4: Prepare the FASTQ raw data
Make sure there only one sample's FASTQ file located in your `[raw FASTQ data path]` directory. The FASTQ file name must follow illumina naming convention rule depicted in this [website](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)
    
    Eg. SampleName_S1_L001_R1_001.fastq.gz, SampleName_S1_L001_R1_001.fastq.gz


#### Step 5: Build the docker image.
1. Switch to the docker DNA-seq-pipeline project-folder:

```
cd ./Onkopipe-main
```

2. Build the Docker Image (named onkopus):

```
docker build -t snakemake_local .

```

#### Step 5: Starting and Usage of the Application:
1. Check if the snakemake available on docker container(set the reference path,fastq input path and pipeline path in advance).

```
docker run -v [ref file path]:/data/ref/ \
-v [raw FASTQ data path]:/data/InputFastqDir/ \
-v [CNVkit bed path]:/data/bed_files/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -v

```

2. Dry run of the snakemake DNA-Seq-Pipeline.(check the processing steps)

```
docker run -v [ref file path]:/data/ref/ \
-v [raw FASTQ data path]:/data/InputFastqDir/ \
-v [CNVkit bed path]:/data/bed_files/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -j all --use-conda -n

```
3. Run pipeline and get the bam and vcf files.

```
docker run -v [ref file path]:/data/ref/ \
-v [raw fastq data path]:/data/InputFastqDir/ \
-v [CNVkit bed path]:/data/bed_files/ \
-v [local snakemake pipeline path]:/work snakemake_local snakemake -j all --use-conda

```

For example:

```
docker run -v /sybig/scratch/Jingyu/DNA_Seq_pipeline/dna_seq_pipeline_snakemakeconda/Ref_data/:/data/ref/ \
-v /sybig/scratch/Jingyu/DNA_Seq_pipeline/NA12878/Fastq_test/:/data/InputFastqDir/ \
-v /sybig/scratch/Jingyu/DNA_Seq_pipeline/Onkopipe_docker/dna_seq_pipeline-master/bed_files/:/data/bed_files/ \
-v $(pwd):/work snakemake_local snakemake -j all --use-conda

```
4. Modify the Snakemake file to meet the specific requirement in need. 
