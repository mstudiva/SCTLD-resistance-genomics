#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s143_R1_001.fastq.gz s143_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s144_R1_001.fastq.gz s144_R2_001.fastq.gz --fastqc -o trimmedReads/
