#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s031_R1_001.fastq.gz s031_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s032_R1_001.fastq.gz s032_R2_001.fastq.gz --fastqc -o trimmedReads/
