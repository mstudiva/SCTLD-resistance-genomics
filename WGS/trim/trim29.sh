#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s059_R1_001.fastq.gz s059_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s060_R1_001.fastq.gz s060_R2_001.fastq.gz --fastqc -o trimmedReads/
