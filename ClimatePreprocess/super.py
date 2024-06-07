import subprocess

# Number of files to process in each job
files_per_job = 2

# Total number of files
total_files = 876

# Calculate the number of jobs
# num_jobs = (total_files + files_per_job - 1) // files_per_job
num_jobs = 3

# Loop over the jobs
for i in range(num_jobs):
    # Calculate the start and end indices for this job
    start_index = i * files_per_job
    end_index = start_index + files_per_job#  - 1

    # Create a job script for this job
    job_script = f"""#!/bin/bash
#SBATCH --qos=priority
#SBATCH --job-name=mT_{i}
#SBATCH --account=impactee
#SBATCH --output=mT_{i}-%j.out
#SBATCH --error=mT_{i}-%j.err
#SBATCH --mail-user=claussar@pik-potsdam.de
#SBATCH --mail-type=fail
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
srun -n 1 python -u CalcClimateVars/calc_monthlyT_ERA5_chunks.py {start_index} {end_index}
"""

    # Write the job script to a file
    with open(f'job_{i}.sh', 'w') as f:
        f.write(job_script)

    # Submit the job
    subprocess.run(['sbatch', f'job_{i}.sh'])