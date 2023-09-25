#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s007_R1_001.fastq.gz s007_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s008_R1_001.fastq.gz s008_R2_001.fastq.gz --fastqc -o trimmedReads/
