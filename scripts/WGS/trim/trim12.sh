#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s025_R1_001.fastq.gz s025_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s026_R1_001.fastq.gz s026_R2_001.fastq.gz --fastqc -o trimmedReads/
