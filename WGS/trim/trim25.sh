#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s051_R1_001.fastq.gz s051_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s052_R1_001.fastq.gz s052_R2_001.fastq.gz --fastqc -o trimmedReads/
