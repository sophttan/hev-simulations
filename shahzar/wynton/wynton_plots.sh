#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -pe smp 30         # 30 slots requested for multiprocessing
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=4G     # job requires up to x GiB of RAM per slot
#$ -l scratch=8G      # job requires up to x GiB of local /scratch space
#$ -l h_rt=100:00:00   # job requires up to x hours of runtime
##$ -t 1-10        # array (serial) job with x tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted

## If you array jobs (option -t), this script will run T times, once per task.
## For each run, $SGE_TASK_ID is set to the corresponding task index (here 1-10).
## To configure different parameters for each task index, one can use a Bash 
## array to map from the task index to a parameter string.

## All possible parameters
# params=(1bac 2xyz 3ijk 4abc 5def 6ghi 7jkl 8mno 9pqr 10stu)

## Select the parameter for the current task index
## Arrays are indexed from 0, so we subtract one from the task index
# param="${params[$((SGE_TASK_ID - 1))]}"

date
hostname

conda deactivate

module load CBI
module load r

for coverage in 75 100
do for static_force in 100 0 50
do for prev in 15 30 50
do for drug in single novel double
do Rscript Run_Schisto_Simulation_Table.R $prev $drug $coverage $static_force TRUE "${NSLOTS-1}"
done
done
done
done

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
