#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s175_R1_001.fastq.gz s175_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s176_R1_001.fastq.gz s176_R2_001.fastq.gz --fastqc -o trimmedReads/
