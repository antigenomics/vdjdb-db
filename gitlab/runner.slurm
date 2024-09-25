#!/bin/bash -l

#SBATCH --job-name=gitlab-ci-runner
#SBATCH --output=gitlab-ci-runner.log
#SBATCH --partition=short
#SBATCH --constraint=compute

# use time expected to run you jobs
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

cd ~/gitlab
~/bin/gitlab-runner run
pwd > ~/my_path.txt

exit 0
