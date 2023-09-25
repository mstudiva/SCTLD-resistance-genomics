#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s133_R1_001.fastq.gz s133_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s134_R1_001.fastq.gz s134_R2_001.fastq.gz --fastqc -o trimmedReads/
