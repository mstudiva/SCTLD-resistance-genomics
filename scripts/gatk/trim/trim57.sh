#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s115_R1_001.fastq.gz s115_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s116_R1_001.fastq.gz s116_R2_001.fastq.gz --fastqc -o trimmedReads/
