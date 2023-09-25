#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s151_R1_001.fastq.gz s151_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s152_R1_001.fastq.gz s152_R2_001.fastq.gz --fastqc -o trimmedReads/
