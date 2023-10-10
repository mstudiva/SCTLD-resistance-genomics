#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s169_R1_001.fastq.gz s169_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s170_R1_001.fastq.gz s170_R2_001.fastq.gz --fastqc -o trimmedReads/
