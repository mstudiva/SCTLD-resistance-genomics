#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s047_R1_001.fastq.gz s047_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s048_R1_001.fastq.gz s048_R2_001.fastq.gz --fastqc -o trimmedReads/
