#!/bin/sh


ListDir=CSASamples
prefix=root://cms-xrd-global.cern.ch/
sampledir=/store/mc/Spring14miniaod/
eosdir=/eos/cms${sampledir}

alias eos='/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

for sample in `eos ls $eosdir`;
do
  listfile=$ListDir/$sample
  if [ -e $listfile ]; then rm $listfile ; fi
  touch $listfile
  for dir1 in `eos ls $eosdir/$sample/MINIAODSIM`;
  do
    for dir2 in `eos ls $eosdir/$sample/MINIAODSIM/$dir1`;
    do
      for file in `eos ls $eosdir/$sample/MINIAODSIM/$dir1/$dir2`;
      do
        echo "$prefix$eosdir/$sample/MINIAODSIM/$dir1/$dir2/$file" >> $listfile
      done   
    done
  done
done


