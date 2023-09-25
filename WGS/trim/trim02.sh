#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s005_R1_001.fastq.gz s005_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s006_R1_001.fastq.gz s006_R2_001.fastq.gz --fastqc -o trimmedReads/
