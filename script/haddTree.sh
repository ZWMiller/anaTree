#!/bin/bash
date

source $ROOTSYS/bin/thisroot.sh

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./haddTree.sh <outName>"
	exit
fi

if [ ! -d out/out_$1_anaTree ]; then
   mkdir out/out_$1_anaTree
fi

for run in `cat runNumberList_run15pp_transverse`
do
   if [ ! -d out/out_$1/${run} ]; then
      continue;
   fi
   ifile=`find out/out_$1/${run}/ -name "st*${run}*anaTree.root" | wc -l`
   if [ $ifile -gt 0 ]; then
      echo "run $run"
      echo "found $ifile anaTree files"
      if [ ! -f out/out_$1_anaTree/st_physics_${run}.anaTree.root ]; then
         hadd -f out/out_$1_anaTree/st_physics_${run}.anaTree.root `find out/out_$1/${run}/ -name "st*${run}*anaTree.root"`
#rm -f out/out_$1/${run}/st*${run}*anaTree.root
      else
         hadd -f out/out_$1_anaTree/st_physics_${run}.anaTree_new.root  out/out_$1_anaTree/st_physics_${run}.anaTree.root `find out/out_$1/${run}/ -name "st*${run}*anaTree.root"`
         mv out/out_$1_anaTree/st_physics_${run}.anaTree_new.root out/out_$1_anaTree/st_physics_${run}.anaTree.root
#rm -f out/out_$1/${run}/st*${run}*anaTree.root
      fi
   fi

done
