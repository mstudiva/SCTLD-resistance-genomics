#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s147_R1_001.fastq.gz s147_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s148_R1_001.fastq.gz s148_R2_001.fastq.gz --fastqc -o trimmedReads/
