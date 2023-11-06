#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s155_R1_001.fastq.gz s155_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s156_R1_001.fastq.gz s156_R2_001.fastq.gz --fastqc -o trimmedReads/
