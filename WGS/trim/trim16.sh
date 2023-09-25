#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s033_R1_001.fastq.gz s033_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s034_R1_001.fastq.gz s034_R2_001.fastq.gz --fastqc -o trimmedReads/
