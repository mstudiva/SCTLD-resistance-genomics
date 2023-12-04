#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s127_R1_001.fastq.gz s127_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s128_R1_001.fastq.gz s128_R2_001.fastq.gz --fastqc -o trimmedReads/
