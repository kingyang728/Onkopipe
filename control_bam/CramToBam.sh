samtools view -@ 2 -b  -T /[[ref file path]]/GRCh38.d1.vd1.fa -o /[control bam path]/NA12878.exome.bam /[control bam path]/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram
samtools index /[control bam path]/NA12878.exome.bam
