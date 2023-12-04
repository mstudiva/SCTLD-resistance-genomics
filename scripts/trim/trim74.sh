#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s149_R1_001.fastq.gz s149_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s150_R1_001.fastq.gz s150_R2_001.fastq.gz --fastqc -o trimmedReads/
