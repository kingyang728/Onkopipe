CNVkit target/antitarget bed file generation
==================================

CNVkit command to generate sequence-accessible coordinates file 

## Authors

* Jingyu Yang

## Usage


#### Step 1: Download gene annotation flat file.
1. refFlat.txt for hg19 ([hg19_refFlat.txt.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz))
2. refFlat.txt for hg38  ([hg38_refFlat.txt.gz](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz))

#### Step 2: Download exome target bed file.
1. hg19_Exome_bed ([Twist_Exome_Core_Covered_Targets_hg19_liftover.bed](https://www.twistbioscience.com/sites/default/files/resources/2022-01/Twist_Exome_Core_Covered_Targets_hg19_liftover.bed))
2. hg38_Exome_bed ([Twist_Exome_Core_Covered_Targets_hg38.bed](https://www.twistbioscience.com/sites/default/files/resources/2022-01/Twist_Exome_Core_Covered_Targets_hg38.bed))

#### Step 3: Download excluded bed file.
1. hg19_excluded_bed ([wgEncodeDukeMapabilityRegionsExcludable.bed](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz))
2. hg38_excluded_bed ([GRCh38_unified_blacklist.bed](https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz))

#### Step 4: Generate target/antitarget bed file
1. Generate access.hg19.bed:

```
cnvkit.py access hg19.fa -x Hg19_dukeExcludeRegions.bed -o access-excludes.hg19.bed

```

```
cnvkit.py access /sybig/scratch/Jingyu/DNA_Seq_pipeline/Ref_dataHg19/genome.fa -x Hg19_dukeExcludeRegions.bed -o access-excludes.hg19.bed

```
2. Generate access.hg38.bed:

```
cnvkit.py access hg38.fa -x GRCh38_unified_blacklist.bed -o access-excludes.hg38.bed

```

```
cnvkit.py access /sybig/scratch/Jingyu/DNA_Seq_pipeline/Ref_data/GRCh38.d1.vd1.fa -x GRCh38_unified_blacklist.bed -o access-excludes.hg38.bed

```


#### Step 5: Generate target/antitarget bed file
1. Generate target_hg19.bed:

```
cnvkit.py target Twist_Exome_Core_Covered_Targets_hg19_liftover.bed --annotate hg19_refFlat.txt -o target_hg19.bed --short-names

```

2. Generate target_hg38.bed:

```
cnvkit.py target Twist_Exome_Core_Covered_Targets_hg38.bed --annotate hg38_refFlat.txt -o target_hg38.bed --short-names

```

3. Generate hg19_antitarget.bed:

```
cnvkit.py antitarget target_hg19.bed -g access-excludes.hg19.bed -o hg19_antitarget.bed

```

4. Generate hg38_antitarget.bed:

```
cnvkit.py antitarget target_hg38.bed -g access-excludes.hg38.bed -o hg38_antitarget.bed

```


