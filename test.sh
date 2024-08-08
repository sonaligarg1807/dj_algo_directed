#$ -cwd -pe nproc 1
#$ -N namd_ham
#$ -l mem=2000M
#$ -e "./$JOB_ID.err"
#$ -o "./$JOB_ID.out"

filename='ham'

echo "Job started on $HOSTNAME at `date` in $scratch"
echo " "

export GMXLIB=/data/fghalami/gromacs-sh-old_Eik/test_plumed/gromacs-sh-old/COUPLED-DYNAMICS/share/top
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/fghalami/gromacs-sh-old_Eik/gromacs-sh-old_Eik/test_plumed/gromacs-sh-old/COUPLED-DYNAMICS/release-tomas-jan2023/lib/
export LD_LIBRARY_PATH=/usr/local/run/plumed-2.5.1/lib:$LD_LIBRARY_PATH


/data/fghalami/gromacs-sh-old_Eik/test_plumed/gromacs-sh-old/COUPLED-DYNAMICS/build-tomas-jan2023/src/kernel/mdrun -ntomp 1 -deffnm $filename > gmx.out
#/data/wxie/Gromacs-SH/build/src/kernel/mdrun-backup -ntomp 1 -deffnm $filename > gmx.out

echo "FINISHED at `date`"

