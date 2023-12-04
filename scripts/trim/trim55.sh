#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s111_R1_001.fastq.gz s111_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s112_R1_001.fastq.gz s112_R2_001.fastq.gz --fastqc -o trimmedReads/
