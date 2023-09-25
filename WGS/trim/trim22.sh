#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s045_R1_001.fastq.gz s045_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s046_R1_001.fastq.gz s046_R2_001.fastq.gz --fastqc -o trimmedReads/
