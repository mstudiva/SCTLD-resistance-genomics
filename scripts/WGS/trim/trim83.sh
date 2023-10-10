#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s167_R1_001.fastq.gz s167_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s168_R1_001.fastq.gz s168_R2_001.fastq.gz --fastqc -o trimmedReads/
