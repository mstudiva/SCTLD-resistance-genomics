#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s015_R1_001.fastq.gz s015_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s016_R1_001.fastq.gz s016_R2_001.fastq.gz --fastqc -o trimmedReads/
