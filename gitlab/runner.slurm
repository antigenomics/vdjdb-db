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

pwd > ~/dluppov/my_path.txt
~/bin/gitlab-runner run
pwd > ~/dluppov/my_path.txt
mkdir ~/sinc_test

exit 0
