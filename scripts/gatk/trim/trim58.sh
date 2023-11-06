#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s117_R1_001.fastq.gz s117_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s118_R1_001.fastq.gz s118_R2_001.fastq.gz --fastqc -o trimmedReads/
