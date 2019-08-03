#! /bin/bash
###
# Sample Torque script for Abaqus
###

# ======= PBS OPTIONS ======= (user input required)
#
### Specify queue to run in
#PBS -q copperhead 
### Set the job name
#PBS -N SubmitINPfile.py
### Specify the # of cpus for your job.
#PBS -l procs=10
#PBS -l mem=30GB
### Adjust walltime below (default walltime = 7 days, or 168 hours)
### if you require > 7 days, INCREASE to estimated # hours needed
### if you DON'T require 7 days DECREASE to estimated # hours needed
### (hint: jobs with smaller walltime value tend to run sooner)
#PBS -l walltime=300:00:00
#
# ===== END PBS OPTIONS =====

# ======= APP OPTIONS ======= (user input required)
#
### (REQUIRED) define inputfile here
#job="SubmitINPfile.py"
job="SubmitINPfile.py"
### uncomment and define user subroutine here (optional)
#user="your-subroutine-name-here"
### specify additional options 
opts=""
#
# ===== END APP OPTIONS =====

# ===== No changes required below here =====

### Get the short $PBS_JOBID
SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`

### Use this to redirect STDOUT and STDERR to working dir
exec 1>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out  2>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.err

### Extract the number of cores requested from the PBS nodefile
cores=$(awk 'END {print NR}' $PBS_NODEFILE)
### Print some node info to the output file
echo -e "\nPBS_NODEFILE contains:\n"
cat $PBS_NODEFILE
echo -e "\nI ran on $HOSTNAME\n"

### Choose TCP or Infiniband by
### Commenting/Uncommenting one of these
export MPI_IC_ORDER=tcp
#export MPI_IC_ORDER=psm

opts="cae nogui=$job" #  $opts interactive ask_delete=off"
if [ ! -z "$user" ]; then
        opts="user=$user $opts"
fi
module load abaqus/6.13-4

# =========== Main Program ===========
# Run Abaqus with  options

cd $PBS_O_WORKDIR
echo "/apps/pkg/abaqus/Commands/abq6134 $opts"
/apps/pkg/abaqus/Commands/abq6134 $opts
