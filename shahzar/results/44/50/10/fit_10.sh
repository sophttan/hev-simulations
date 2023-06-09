#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -pe smp 64
#$ -l mem_free=2G     # job requires up to 2 GiB of RAM per slot
#$ -l scratch=4G      # job requires up to 4 GiB of local /scratch space
#$ -l h_rt=360:00:00   # job requires up to 100 hours of runtime
##$ -t 1-3           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted


date
hostname

module load CBI
module load r

Rscript fit_10.R

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"