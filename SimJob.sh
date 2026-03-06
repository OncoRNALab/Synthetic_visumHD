#!/bin/bash
#PBS -N 100ksim
#PBS -l nodes=1:ppn=8
#PBS -l walltime=06:00:00
#PBS -l mem=64gb

source activate /user/gent/446/vsc44685/DataVO_dir/miniconda2/envs/stdata_sim
cd /user/gent/446/vsc44685/ScratchVO_dir/STdata_sim
ml Nextflow/25.04.8

export APPTAINER_TMPDIR="${VSC_SCRATCH_VO_USER}/.apptainer/tmp"
export APPTAINER_CACHEDIR="${VSC_SCRATCH_VO_USER}/.apptainer/cache"
export SINGULARITY_CACHEDIR="${VSC_SCRATCH_VO_USER}/singularity"
export NXF_SINGULARITY_CACHEDIR="${VSC_SCRATCH_VO_USER}/.apptainer/cache"

nextflow run main.nf \
    -profile conda_named \
    --feature_h5 /user/gent/446/vsc44685/ScratchVO_dir/spaceranger/Visium_HD_Mouse_Brain_Fresh_Frozen_feature_slice.h5 \
    --sim_shards 8 


