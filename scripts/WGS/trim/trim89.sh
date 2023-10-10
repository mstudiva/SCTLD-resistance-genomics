#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s179_R1_001.fastq.gz s179_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s180_R1_001.fastq.gz s180_R2_001.fastq.gz --fastqc -o trimmedReads/
