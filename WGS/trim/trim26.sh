#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s053_R1_001.fastq.gz s053_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s054_R1_001.fastq.gz s054_R2_001.fastq.gz --fastqc -o trimmedReads/
