#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s173_R1_001.fastq.gz s173_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s174_R1_001.fastq.gz s174_R2_001.fastq.gz --fastqc -o trimmedReads/
