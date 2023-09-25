#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s091_R1_001.fastq.gz s091_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s092_R1_001.fastq.gz s092_R2_001.fastq.gz --fastqc -o trimmedReads/
