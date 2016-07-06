#!/bin/bash
date

if [ $# -lt 2 ]; then
    echo "not enough arguements"
    echo "./submitAll.sh <outName> <mode: 0-mb, 1-ht, 2-mtd>"
	exit
fi

if [ ! -d out ]; then
    mkdir out
fi
if [ ! -d log ]; then
    mkdir log
fi

dir=$(pwd)
ver=$1
echo $dir
if [ ! -d ${dir}/out/out_$1 ]; then
	mkdir ${dir}/out/out_$1
fi
if [ ! -d ${dir}/log/log_$1 ]; then
	mkdir ${dir}/log/log_$1
fi

if [ ! -d submit ]; then
   mkdir submit
fi

if [ ! -d submit/$1 ]; then
   mkdir submit/$1
fi

cd submit/$1
 
if [ -f runNumber_$1_$2 ]; then
   rm runNumber_$1_$2
fi

get_file_list.pl -keys 'runnumber' -cond 'production=P16id,trgsetupname=production_pp200trans_2015,filetype=daq_reco_mudst,filename~st_physics,storage=nfs' -limit 9 >& runNumber_$1_$2
#get_file_list.pl -keys 'runnumber' -cond 'production=P15ic,collision=auau200,trgsetupname=AuAu_200_production_2014||AuAu_200_production_mid_2014||AuAu_200_production_low_2014,filename~st_physics,filetype=daq_reco_Mudst,storage=local' -limit 0 >& runNumber_$1_$2
cp runNumber_$1_$2 ../.

echo $dir
for runId in `cat runNumber_$1_$2`
do
   if [ ! -d $dir/submit/$1/$runId ]; then
      mkdir $dir/submit/$1/$runId
   fi
   cp $dir/StRoot/script/pp200.xml $dir/submit/$1/$runId/.
   cd $dir/submit/$1/$runId/
   if [ ! -d ${dir}/out/out_$1/$runId ]; then
   	mkdir ${dir}/out/out_$1/$runId
   fi
   if [ ! -d ${dir}/log/log_$1/$runId ]; then
   	mkdir ${dir}/log/log_$1/$runId
   fi
   star-submit-template -template pp200.xml -entities path=${dir},ver=$1,mode=$2,ptype=0,run=$runId
done
