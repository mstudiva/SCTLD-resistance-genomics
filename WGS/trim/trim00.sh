#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s001_R1_001.fastq.gz s001_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s002_R1_001.fastq.gz s002_R2_001.fastq.gz --fastqc -o trimmedReads/
