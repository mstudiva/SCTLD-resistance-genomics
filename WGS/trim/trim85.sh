#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s171_R1_001.fastq.gz s171_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s172_R1_001.fastq.gz s172_R2_001.fastq.gz --fastqc -o trimmedReads/
