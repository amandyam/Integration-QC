#!/usr/bin/bash

#To use absolute path for fastq files, if the folder was more than 1 sample, create a list as follow:
JOBID="drgcellranger"
SAMPLE_IDS="SRR21958173"
#, SRR21958171"
TRANSCRIPTOME="/storage/chentemp/atulyadm/cellranger-7.1.0/refdata-gex-mm10-2020-A"
FASTQS="/storage/chentemp/atulyadm/cellrangertutorial/fastq"

cellranger count --id=$JOBID \
                 --transcriptome=$TRANSCRIPTOME \
                 --fastqs=$FASTQS \
                 --sample=$SAMPLE_IDS \
                 --localcores=16 \
                 --localmem=64 \
                 --expect-cells=10000


echo "All processes for all samples were done !!"
