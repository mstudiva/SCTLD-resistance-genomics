#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s177_R1_001.fastq.gz s177_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s178_R1_001.fastq.gz s178_R2_001.fastq.gz --fastqc -o trimmedReads/
