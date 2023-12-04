#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s037_R1_001.fastq.gz s037_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s038_R1_001.fastq.gz s038_R2_001.fastq.gz --fastqc -o trimmedReads/
