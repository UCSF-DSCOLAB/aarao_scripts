#PBS -e /data/shared/krummellab/arrao/logs/run_scRNASeq_stg3.err
#PBS -o /data/shared/krummellab/arrao/logs/run_scRNASeq_stg3.out
#PBS -N run_scRNASeq
#PBS -l nodes=1:ppn=1
#PBS -l vmem=150gb

/home/arrao/venvs/scRNASeq/bin/python /home/arrao/src/scRNASeq/process_snRNAseq.py \
    --input_list /data/shared/krummellab/arrao/projects/immunox/sample_default_3.0.0.tsv \
    --reference_folder /data/shared/krummellab/arrao/projects/immunox/references \
    --important_celltypes_file /data/shared/krummellab/arrao/projects/immunox/important_cell_types.tsv \
    --stage 3 \
    --results_dir /data/shared/krummellab/arrao/projects/immunox/results/ \
    --overwrite \
    --snRNAseq
