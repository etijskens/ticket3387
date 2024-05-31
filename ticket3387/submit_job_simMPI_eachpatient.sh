#!/bin/bash
#SBATCH --job-name=batchrun_job       # Job name
#SBATCH --account=ap_test             # Project account
#SBATCH --partition=skylake           # Use the 'skylake' partition
#SBATCH --output=output_%j.txt        # Standard output and error log
#SBATCH --error=error_%j.txt          # Error log
#SBATCH --ntasks=1                    # Number of MPI tasks
#SBATCH --cpus-per-task=28            # Number of CPU cores per task
#SBATCH --time=7-00:00:00             # Time limit days-hrs:min:sec

# Load necessary modules
module load calcua/2023a
module load Python/3.11.3-GCCcore-12.3.0
module load SciPy-bundle/2023.07-gfbf-2023a
#module load matplotlib/3.7.2-gfbf-2023a
module load openpyxl/3.1.2-GCCcore-12.3.0

export PYTHONPATH="${VSC_DATA}/mesa_lib/lib/python3.11/site-packages/:${PYTHONPATH}"

# Run the Python script
python batch_run_simMPI_eachpatient.py

