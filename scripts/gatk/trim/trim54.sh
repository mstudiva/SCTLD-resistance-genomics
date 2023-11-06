#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s109_R1_001.fastq.gz s109_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s110_R1_001.fastq.gz s110_R2_001.fastq.gz --fastqc -o trimmedReads/
