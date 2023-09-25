#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s009_R1_001.fastq.gz s009_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s010_R1_001.fastq.gz s010_R2_001.fastq.gz --fastqc -o trimmedReads/
