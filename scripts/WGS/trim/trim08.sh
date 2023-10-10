#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s017_R1_001.fastq.gz s017_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s018_R1_001.fastq.gz s018_R2_001.fastq.gz --fastqc -o trimmedReads/
