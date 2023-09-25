#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s093_R1_001.fastq.gz s093_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s094_R1_001.fastq.gz s094_R2_001.fastq.gz --fastqc -o trimmedReads/
