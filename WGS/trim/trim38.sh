#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s077_R1_001.fastq.gz s077_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s078_R1_001.fastq.gz s078_R2_001.fastq.gz --fastqc -o trimmedReads/
