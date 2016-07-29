#!/bin/bash
date

source $ROOTSYS/bin/thisroot.sh

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./haddPurity.sh <outName>"
	exit
fi

if [ ! -d out/out_$1_purity ]; then
   mkdir out/out_$1_purity
fi

for run in `cat runNumberList_run15pp_transverse`
do
  echo "RUN"
   if [ ! -d out/out_$1/${run} ]; then
      continue;
   fi
   ifile=`find out/out_$1/${run}/ -name "st*${run}*purity.root" | wc -l`
   if [ $ifile -gt 0 ]; then
      echo "run $run"
      echo "found $ifile purity files"
      if [ ! -f out/out_$1_purity/st_physics_${run}.purity.root ]; then
         hadd -f out/out_$1_purity/st_physics_${run}.purity.root `find out/out_$1/${run}/ -name "st*${run}*purity.root"`
#rm -f out/out_$1/${run}/st*${run}*purity.root
      else
         hadd -f out/out_$1_purity/st_physics_${run}.purity_new.root  out/out_$1_purity/st_physics_${run}.purity.root `find out/out_$1/${run}/ -name "st*${run}*purity.root"`
         mv out/out_$1_purity/st_physics_${run}.purity_new.root out/out_$1_purity/st_physics_${run}.purity.root
#rm -f out/out_$1/${run}/st*${run}*purity.root
      fi
   fi

echo "DONE"

done
