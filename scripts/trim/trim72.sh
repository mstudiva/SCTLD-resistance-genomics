#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s145_R1_001.fastq.gz s145_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s146_R1_001.fastq.gz s146_R2_001.fastq.gz --fastqc -o trimmedReads/
