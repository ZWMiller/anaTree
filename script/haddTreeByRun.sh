#!/bin/bash
date

source $ROOTSYS/bin/thisroot.sh

if [ $# -lt 2 ]; then
    echo "not enough arguements"
    echo "./haddTreeByRun.sh <outName> <run>"
	exit
fi

#path=/star/u/huangbc/data01/pico/anaTree/20151125_v4/prod
path=/gpfs/mnt/gpfs01/star/pwg/zamiller/run15ppAnaTree/prod/anaTreeMaker_v1_070616
target=$(pwd)

if [ ! -d ${target}/out ]; then
   mkdir ${target}/out
fi

if [ ! -d ${target}/out/out_$1_anaTree ]; then
   mkdir ${target}/out/out_$1_anaTree
fi

run=$2
if [ ! -d ${path}/out/out_$1/${run} ]; then
   echo "no run $run"
   exit
fi
ifile=`find ${path}/out/out_$1/${run} -name "st*${run}*anaTree.root" | wc -l`
if [ $ifile -gt 0 ]; then
   echo "run ${run}"
   echo "found $ifile anaTree files"
   if [ ! -f ${target}/out/out_$1_anaTree/st_physics_${run}.anaTree.root ]; then
      if [ ! -f ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish ]; then
         touch ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish
         hadd -f ${target}/out/out_$1_anaTree/st_physics_${run}.anaTree.root `find ${path}/out/out_$1/${run} -name "st*${run}*anaTree.root"`
         rm -f ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish
#rm -f ${path}/out/out_$1/${run}/st*${run}*anaTree.root
      fi
   else
      if [ ! -f ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish ]; then
         touch ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish
         hadd -f ${target}/out/out_$1_anaTree/st_physics_${run}.anaTree_new.root  ${path}/out/out_$1_anaTree/st_physics_${run}.anaTree.root `find ${path}/out/out_$1/${run}  -name "st*${run}*anaTree.root"`
         mv ${target}/out/out_$1_anaTree/st_physics_${run}.anaTree_new.root ${target}/out/out_$1_anaTree/st_physics_${run}.anaTree.root
         rm -f ${target}/out/out_$1_anaTree/${run}_anaTree_unfinish
#rm -f ${path}/out/out_$1/${run}/st*${run}*anaTree.root
      fi
   fi
fi
