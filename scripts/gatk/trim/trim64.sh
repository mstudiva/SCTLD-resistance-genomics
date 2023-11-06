#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s129_R1_001.fastq.gz s129_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s130_R1_001.fastq.gz s130_R2_001.fastq.gz --fastqc -o trimmedReads/
