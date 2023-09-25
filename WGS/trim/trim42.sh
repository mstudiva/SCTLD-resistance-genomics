#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s085_R1_001.fastq.gz s085_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s086_R1_001.fastq.gz s086_R2_001.fastq.gz --fastqc -o trimmedReads/
