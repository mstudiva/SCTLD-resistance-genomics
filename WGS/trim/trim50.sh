#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s101_R1_001.fastq.gz s101_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s102_R1_001.fastq.gz s102_R2_001.fastq.gz --fastqc -o trimmedReads/
