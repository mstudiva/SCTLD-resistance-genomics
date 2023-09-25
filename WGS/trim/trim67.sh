#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s135_R1_001.fastq.gz s135_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s136_R1_001.fastq.gz s136_R2_001.fastq.gz --fastqc -o trimmedReads/
