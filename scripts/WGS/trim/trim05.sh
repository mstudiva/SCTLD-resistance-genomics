#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s011_R1_001.fastq.gz s011_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s012_R1_001.fastq.gz s012_R2_001.fastq.gz --fastqc -o trimmedReads/
