#! /bin/bash

################################################################################
### DON'T TOUCH (EXCEPT ON RARE OCCASION AS NEEDED)
this_script=`basename "$0"`
today=$(date +'%Y_%m_%d')
ncpus=128       # CPUs per node on Derecho (default 128)
user=$USER      # $USER is an NCAR-exported environment variable
acct=UCUB0166
# $HOME=/glade/u/home/$USER is an NCAR-exported environmnent variable
# $SCRATCH=/glade/derecho/scratch/$USER is an NCAR-exported environment variable
TOP_DIR=${SCRATCH}/NCAR-LES
NOTEBOOK=${HOME}/LAB_NOTEBOOK/NCAR-LES

################################################################################
### EDIT AT WILL FOR EACH RUN
compiler=intel       # options are [cray, intel, gnu]
compile_mode=profile   # options useful on HPC are [debug, profile, fast]
project=classic_smag # creates sub-directory within NCAR-LES folders
job_name=${compiler}_${compile_mode}_no_coriolis  # should be something unique at least for today
RUN_DIR=${TOP_DIR}/${project}/${today}/${job_name}
NOTES=${NOTEBOOK}/${project}/${today}/${job_name}
outfile=logfile.out
ntasks=32 # total tasks wanted (default 32)

# git_remote=origin       # fancier git stuff coming in a future update
# git_branch=classic_smag # fancier git stuff coming in a future update

### AUTOMATIC
# This bash math allows ntasks to not be a multiple of ncpus (for instance ntasks=64 or 192)
(( nnodes= ntasks/ncpus + ((ntasks % ncpus)>0) ))

#######################################################################
echo '---> Creating directory...'
# I'll let you decide if you really want to delete a reused run directory.
# rm -rf $RUN_DIR
mkdir -p ${RUN_DIR}/data

echo '---> Loading modules of specified compiler...'
if [[ "$compiler" == "cray" ]]; then
    module load cce
elif [[ "$compiler" == "gnu" ]]; then
    module load gnu
else
    module load intel
fi

echo '---> Compiling...'
make clean
make COMPILER=$compiler $compile_mode
cp bin/lesmpi.exe $RUN_DIR

echo '---> Creating lab notebook entry...'
mkdir -p $NOTES
cp bin/lesmpi.exe $NOTES/
cp $this_script $NOTES/
echo '------------------------------------------' >> ${NOTES}/version_info
echo '              last git commit             ' >> ${NOTES}/version_info
echo '------------------------------------------' >> ${NOTES}/version_info
git show -q >> ${NOTES}/version_info
echo '------------------------------------------' >> ${NOTES}/version_info
echo '             current git status           ' >> ${NOTES}/version_info
echo '------------------------------------------' >> ${NOTES}/version_info
git status >> ${NOTES}/version_info

#######################################################################
echo '---> Making EXEC_STEP'
cd $RUN_DIR
pwd

# EXECUTE SCRIPT
cat > EXEC_STEP << EXEC
#!/bin/sh
#PBS -N $job_name
#PBS -l walltime=02:00:00
#PBS -l select=$nnodes:ncpus=$ncpus:mpiprocs=$ncpus
#PBS -A $acct
#PBS -q develop
####PBS -l job_priority=economy
#PBS -o exec.out
#PBS -e exec.err

if [[ "$compiler" == "cray" ]]; then
    module load cce
elif [[ "$compiler" == "gnu" ]]; then
    module load gnu
else
    module load intel
fi

cp \${PBS_NODEFILE} $NOTES/

mpiexec -n $ntasks ./lesmpi.exe > $outfile
wait

cp exec.out exec.err $outfile $NOTES/

EXEC
cp EXEC_STEP $NOTES/

echo '---> Running' $RUN_DIR
qsub < EXEC_STEP
