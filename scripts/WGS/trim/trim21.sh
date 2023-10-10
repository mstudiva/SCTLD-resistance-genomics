#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s043_R1_001.fastq.gz s043_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s044_R1_001.fastq.gz s044_R2_001.fastq.gz --fastqc -o trimmedReads/
