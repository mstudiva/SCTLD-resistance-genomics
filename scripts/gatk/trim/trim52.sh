#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s105_R1_001.fastq.gz s105_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s106_R1_001.fastq.gz s106_R2_001.fastq.gz --fastqc -o trimmedReads/
