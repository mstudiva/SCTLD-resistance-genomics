#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s041_R1_001.fastq.gz s041_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s042_R1_001.fastq.gz s042_R2_001.fastq.gz --fastqc -o trimmedReads/
