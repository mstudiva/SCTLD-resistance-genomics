#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s029_R1_001.fastq.gz s029_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s030_R1_001.fastq.gz s030_R2_001.fastq.gz --fastqc -o trimmedReads/
