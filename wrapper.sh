export INPUT=$1 # The path to the BIRDbath input file
export OUTPUT=$2 # The path to the BIRDbath output file
export STARTVAL=$3 # The 0-based index of the first variant in the BIRDbath input file
export ENDVAL=$4 # The 0-based index of the last variant in the BIRDbath input file
export BIRD=$5 # The path to the directory with the BIRDbath model
export STAN=$6 # The path to the directory where STAN is installed
export CONC=0 # The concentration parameter

export DIR=$(pwd)

cd ${BIRD}

export SAVEDIR=${DIR}
export TMPDIR=${DIR}

./BIRDbath.py ./BIRDbath 1.25 $CONC ${INPUT} $SAVEDIR/${OUTPUT} 1000 ${STARTVAL}-${ENDVAL} > $SAVEDIR/${OUTPUT}
