#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s081_R1_001.fastq.gz s081_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s082_R1_001.fastq.gz s082_R2_001.fastq.gz --fastqc -o trimmedReads/
