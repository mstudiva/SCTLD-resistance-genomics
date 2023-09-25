#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s165_R1_001.fastq.gz s165_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s166_R1_001.fastq.gz s166_R2_001.fastq.gz --fastqc -o trimmedReads/
