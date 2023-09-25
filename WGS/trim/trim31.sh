#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s063_R1_001.fastq.gz s063_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s064_R1_001.fastq.gz s064_R2_001.fastq.gz --fastqc -o trimmedReads/
