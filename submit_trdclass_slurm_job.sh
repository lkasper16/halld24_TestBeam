#!/bin/bash

if [ -z "$1" ]; then
    echo "Error: Please provide the run number."
    echo "Usage: $0 <run_number>"
    exit 1
fi

RUNNUMBER=${1}

# Define executable path
EXECUTABLE=/work/halld2/home/nseptian/halld24_TestBeam/run_trdclass_halld24.sh

INPUT_DIR="/work/halld2/home/nseptian/halld24_TestBeam/ROOT/"

echo "INPUT_DIR: $INPUT_DIR"

# Create a list of input files
cd $INPUT_DIR
INPUT_FILES=($(ls Run_00${RUNNUMBER}_*.root))
# INPUT_FILES=($(ls Run_00${RUNNUMBER}_001.root))

# Create a cancel script
CANCEL_SCRIPT="/work/halld2/home/nseptian/halld24_TestBeam/cancel_trdclass_jobs_${RUNNUMBER}.sh"
echo "#!/bin/bash" > $CANCEL_SCRIPT

# Loop through each input file and create a separate job script
for INPUT_FILE in "${INPUT_FILES[@]}"; do
    INPUT_PATH="$INPUT_DIR/$INPUT_FILE"
    echo "Making job script for $INPUT_FILE"
    JOB_SCRIPT=$(mktemp /tmp/submit_trdclass_job_XXXXXX.sh)
    
    cat <<EOT > $JOB_SCRIPT
#!/bin/bash
#SBATCH --job-name=TRDClass_${INPUT_FILE}  # Job name
#SBATCH --output=$INPUT_DIR/output_trdclass_${INPUT_FILE}.txt  # Output file name
#SBATCH --error=$INPUT_DIR/error_trdclass_${INPUT_FILE}.txt    # Error file name
#SBATCH --time=1:00:00                 # Time limit (hh:mm:ss)
#SBATCH --partition=production               # Partition name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Total number of tasks
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=1G                        # Memory per task

# Run the executable
cd /work/halld2/home/nseptian/halld24_TestBeam/
$EXECUTABLE $INPUT_DIR $INPUT_FILE
EOT

    # Submit the job script and capture the job ID
    echo "Submitting job for $INPUT_FILE"
    echo "job command = "$EXECUTABLE $INPUT_DIR $INPUT_FILE
    JOB_ID=$(sbatch $JOB_SCRIPT | awk '{print $4}')

    # Add the job ID to the cancel script
    echo "scancel $JOB_ID" >> $CANCEL_SCRIPT

    # Remove the job script
    rm $JOB_SCRIPT
done

chmod +x $CANCEL_SCRIPT