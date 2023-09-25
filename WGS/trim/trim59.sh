#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s119_R1_001.fastq.gz s119_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s120_R1_001.fastq.gz s120_R2_001.fastq.gz --fastqc -o trimmedReads/
