#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s073_R1_001.fastq.gz s073_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s074_R1_001.fastq.gz s074_R2_001.fastq.gz --fastqc -o trimmedReads/
