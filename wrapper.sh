#!/bin/bash
#SBATCH -e BIRDbath.err
#SBATCH --mem-per-cpu=2G
#SBATCH -c 1


export DIR=$1 # The directory with the BIRDbath input files, where the BIRDbath output files will be generated

cd $2 # The directory with the BIRDbath model

export SAVEDIR=${DIR}/outputs
export BIRD=$2 # The directory with the BIRDbath model
export PYTHONPATH=$3 # The directory with the GitHub repository: https://github.com/bmajoros/python
export TMPDIR=$DIR/temp
export STAN=$4 # The directory where STAN is installed
export INPUT=$DIR/$5 # The name of the input file
export STARTVAL=$6 # The zero-based index of the first variant
export ENDVAL=$7 # The zero-based index of the last variant

./estimate-conc.py ./estimate-conc 1.25 ${INPUT} ${DIR}/out.txt 1000 ${STARTVAL}-${ENDVAL} | cut -f2 > ${DIR}/tmp.2

./median.pl ${DIR}/tmp.2 > ${DIR}/conc

export CONC=$(cat ${DIR}/conc)

./BIRDbath1.py ./BIRDbath1 1.25 $CONC ${INPUT} $SAVEDIR/${OUTPUT} 1000 ${STARTVAL}-${ENDVAL} > $SAVEDIR/out.txt