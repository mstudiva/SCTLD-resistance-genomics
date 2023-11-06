#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s123_R1_001.fastq.gz s123_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s124_R1_001.fastq.gz s124_R2_001.fastq.gz --fastqc -o trimmedReads/
