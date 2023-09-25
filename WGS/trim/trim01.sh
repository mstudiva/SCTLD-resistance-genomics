#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s003_R1_001.fastq.gz s003_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s004_R1_001.fastq.gz s004_R2_001.fastq.gz --fastqc -o trimmedReads/
