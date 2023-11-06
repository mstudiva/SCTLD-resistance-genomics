#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s071_R1_001.fastq.gz s071_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s072_R1_001.fastq.gz s072_R2_001.fastq.gz --fastqc -o trimmedReads/
