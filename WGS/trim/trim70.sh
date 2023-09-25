#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s141_R1_001.fastq.gz s141_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s142_R1_001.fastq.gz s142_R2_001.fastq.gz --fastqc -o trimmedReads/
