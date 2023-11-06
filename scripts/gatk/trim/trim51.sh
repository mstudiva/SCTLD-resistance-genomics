#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s103_R1_001.fastq.gz s103_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s104_R1_001.fastq.gz s104_R2_001.fastq.gz --fastqc -o trimmedReads/
