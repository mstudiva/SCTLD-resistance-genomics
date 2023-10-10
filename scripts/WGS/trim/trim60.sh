#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s121_R1_001.fastq.gz s121_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s122_R1_001.fastq.gz s122_R2_001.fastq.gz --fastqc -o trimmedReads/
