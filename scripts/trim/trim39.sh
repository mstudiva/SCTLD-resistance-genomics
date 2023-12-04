#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s079_R1_001.fastq.gz s079_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s080_R1_001.fastq.gz s080_R2_001.fastq.gz --fastqc -o trimmedReads/
