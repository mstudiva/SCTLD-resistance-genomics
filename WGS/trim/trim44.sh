#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s089_R1_001.fastq.gz s089_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s090_R1_001.fastq.gz s090_R2_001.fastq.gz --fastqc -o trimmedReads/
