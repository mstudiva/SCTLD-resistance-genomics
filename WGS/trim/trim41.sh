#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s083_R1_001.fastq.gz s083_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s084_R1_001.fastq.gz s084_R2_001.fastq.gz --fastqc -o trimmedReads/
