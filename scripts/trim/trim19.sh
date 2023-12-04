#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s039_R1_001.fastq.gz s039_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s040_R1_001.fastq.gz s040_R2_001.fastq.gz --fastqc -o trimmedReads/
