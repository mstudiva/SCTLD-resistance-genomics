#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s137_R1_001.fastq.gz s137_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s138_R1_001.fastq.gz s138_R2_001.fastq.gz --fastqc -o trimmedReads/
