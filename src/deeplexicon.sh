#!/bin/bash
#BSUB -R "rusage[mem=128] span[hosts=1]"
#BSUB -J deeplexicon[1]
#BSUB -o deep_o_%J.log
#BSUB -e deep_e_%J.log
#BSUB -n 4
#BSUB -q cranio 

# to demultiplex dRNA-seq libraries made with Deeplexicon barcodes per PMID: 32907883

# Deeplexicon docker image gave us a lot of trouble, so this is the run script for
# Deeplexicon singularity run on an LSF job server
# running this does not make the whole cluster unstable, yay

cd # insert path to directory where you want to run deeplexicon here 

. /usr/share/Modules/init/bash
module load modules modules-init
module load singularity

singularity exec ~/deeplex_sing.sif /deeplexicon/deeplexicon_multi.py dmux --threads 4 \
    -p path/to/nanopore_run/fast5_all/ \
    -m /deeplexicon/models/resnet20-final.h5 > dmuxed.tsv


