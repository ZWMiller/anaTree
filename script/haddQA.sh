#!/bin/bash
date

source $ROOTSYS/bin/thisroot.sh

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./haddQA.sh <outName>"
	exit
fi

if [ ! -d out/out_$1_qa ]; then
   mkdir out/out_$1_qa
fi

for run in `cat runNumberList_run15pp_transverse`
do
  echo "RUN"
   if [ ! -d out/out_$1/${run} ]; then
      continue;
   fi
   ifile=`find out/out_$1/${run}/ -name "st*${run}*qa.root" | wc -l`
   if [ $ifile -gt 0 ]; then
      echo "run $run"
      echo "found $ifile qa files"
      if [ ! -f out/out_$1_qa/st_physics_${run}.qa.root ]; then
         hadd -f out/out_$1_qa/st_physics_${run}.qa.root `find out/out_$1/${run}/ -name "st*${run}*qa.root"`
#rm -f out/out_$1/${run}/st*${run}*qa.root
      else
         hadd -f out/out_$1_qa/st_physics_${run}.qa_new.root  out/out_$1_qa/st_physics_${run}.qa.root `find out/out_$1/${run}/ -name "st*${run}*qa.root"`
         mv out/out_$1_qa/st_physics_${run}.qa_new.root out/out_$1_qa/st_physics_${run}.qa.root
#rm -f out/out_$1/${run}/st*${run}*qa.root
      fi
   fi

echo "DONE"

done
