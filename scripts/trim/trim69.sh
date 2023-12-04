#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s139_R1_001.fastq.gz s139_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s140_R1_001.fastq.gz s140_R2_001.fastq.gz --fastqc -o trimmedReads/
