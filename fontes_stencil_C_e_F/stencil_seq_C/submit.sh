#!/bin/ksh 
# 
# script to run BRAMS on pure MPI at tupa   v82
#

if [ $# -lt 8 ]
then
   echo ""
   echo "Falta parametros..."
   echo ""
   echo "usage: submit.sh total_ranks ranks/no qtde_nos expNum Proc OMP/MPI sourcecode.ext ext"
   echo ""
   echo "where:"
   echo ""
   echo "1. total_ranks:  total number of MPI ranks"
   echo "2. ranks/no:     number of ranks per node"
   echo "3. qtde_nos:     number of nodes"
   echo "4. expNum:       label number of experiment"
   echo "5. Proc:         HASWELL or BROADWELL"
   echo "6. OMP/MPI:      flag to run OMP exec or not: OMP or MPI"
   echo "7. sourcecode:   source.f90 or source.c"
   echo "8. extension:    f90 or c"
   echo ""
   exit
fi

#
# script arguments
#

export energy=1
export ftnopt="-O3"
export TotMpi=${1}
export MpiNod=${2}
export Nodes=${3}
export expnum=${4}
export procfila=${5}
export codename=${7}
export extension=${8}
export comp=${9}
export nxgrid=${10}
export niters=${11}

export px=$(echo "sqrt(${TotMpi})" | bc)
export py=$(echo "sqrt(${TotMpi})" | bc)
export RamsIn=RAMSIN
export Exp=$(date -u +"%Y%m%d%H%M%S")
export Gov=powersave   #  nao aumenta a freq automaticamente, vai operar na freq mais baixa possivel.
#export Gov=performance  # aumenta a frequencia automatticamente

if [ ${6} == "OMP" ]
then
   export ompflag="_OMP"
   export jnameflag="omp"
else
   export ompflag=""
   export jnameflag="mpi"
fi

if [ ${procfila} == "BROADWELL" ]
then
  export fila="small"
else
  export fila="hw16"
fi  
#fila="large"



#
# job name 
#
RunName=${codename}${TotMpi}
OutName=Out_${codename}_${TotMpi}M_${Exp}_${nxgrid}_${TotMpi}-${MpiNod}-${Nodes}_${comp}
DumpOut=Dump${OutName}
#
# directories
# executable full path
#
export DirBase=$(pwd)
cd ${DirBase}
export executable=${DirBase}/${codename}${TotMpi}${MpiNod}${Nodes}${comp}_${Exp}.x

if [ ${extension} == "f90" ]
then
  # Fortran:
  echo "compilando... Fortran"
  rm -f ${executable}
  echo "ftn ${ftnopt} ${codename} -o ${executable} "
  ftn ${ftnopt} ${codename} -o ${executable} 
  echo "compilado!"
  #sleep 5
  if [ ! -s ${executable} ]
  then
     exit
  fi
else
  # C:
  echo "compilando... C"
  rm -f ${executable}
  echo "cc ${ftnopt} ${codename} -o ${executable} -lm"
  cc ${ftnopt} ${codename} -o ${executable} -lm
  echo "compilado!" 
  #sleep 5
  if [ ! -s ${executable} ]
  then
    exit
  fi
fi

mkdir -p ${DirBase}/exp${expnum}

echo "Fila: ${fila} ${procfila}"

#
# starts producing queue script file qsub.sh
#
cat <<EOF0> qsub${Exp}.sh
#!/bin/bash
#PBS -l walltime=01:59:00
##PBS -l nodes=${Nodes}:ppn=${MpiNod}
#PBS -l nodes=${Nodes}:ppn=44
###PBS -l mppdepth=${OmpMpi}
###PBS -l mppnppn=${Nppn}
#PBS -N ${RunName}
#PBS -j oe
#PBS -o ${DirBase}/exp${expnum}/${OutName}
#PBS -q ${fila} 

cd ${DirBase}

export OMP_NUM_THREADS=${OmpMpi}

#export PAT_RT_SUMMARY=1 # off
#export PAT_RT_SUMMARY=0 # enable runtime summarization and data aggreation
#export ATP_ENABLED=1

ulimit -m unlimited
ulimit -s unlimited
ulimit -v unlimited
ulimit -c unlimited



aprun -b -n ${TotMpi} -N ${MpiNod} -d 1 --p-governor ${Gov} ${executable} ${nxgrid} ${energy} ${niters} ${px} ${py}
rm -f ${executable}
EOF0

#
# finishes producing file qsub.sh and moves to executable directory
#
chmod +x qsub${Exp}.sh
#
# qsub with variable # PEs per node
#

qsub qsub${Exp}.sh 
mv qsub${Exp}.sh ${DirBase}/exp${expnum}/qsub.sh_${Exp}
cp -f ${codename} ${DirBase}/exp${expnum}/${codename}_${Exp}
echo "$(date)"

exit 0

