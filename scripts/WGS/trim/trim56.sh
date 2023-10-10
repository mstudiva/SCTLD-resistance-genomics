#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s113_R1_001.fastq.gz s113_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s114_R1_001.fastq.gz s114_R2_001.fastq.gz --fastqc -o trimmedReads/
