#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s131_R1_001.fastq.gz s131_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s132_R1_001.fastq.gz s132_R2_001.fastq.gz --fastqc -o trimmedReads/
