#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s075_R1_001.fastq.gz s075_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s076_R1_001.fastq.gz s076_R2_001.fastq.gz --fastqc -o trimmedReads/
