#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s153_R1_001.fastq.gz s153_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s154_R1_001.fastq.gz s154_R2_001.fastq.gz --fastqc -o trimmedReads/
