#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s107_R1_001.fastq.gz s107_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s108_R1_001.fastq.gz s108_R2_001.fastq.gz --fastqc -o trimmedReads/
