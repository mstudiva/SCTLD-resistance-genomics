#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s049_R1_001.fastq.gz s049_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s050_R1_001.fastq.gz s050_R2_001.fastq.gz --fastqc -o trimmedReads/
