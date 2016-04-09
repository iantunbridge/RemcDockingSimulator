#!/bin/bash
echo "Cluster <fileindex> <symmetric clusters, 0/1> <trajactoryname.pdb> <write first N clusters> <minimum members per cluster> <skip> <molecules per instance> <list of lengths>"
# take the indexfile and generate all instances for clustering.
#./pdbgen $1 $2
#create trajectory file
cat pdbgen_results/* >> $3
#fix trajectory file format
./sanitiseTrajectory.sh $3
#cluster
g_cluster -f $3 -s $3 -wcl $4 -nst $5 -cutoff 0.1 -cl -sz -skip $6
#create input file for reference structures, assume the same as structures in indexfile, otherwise construct your own
cat $1 | grep data | awk '{print $2}' > processTraj.exp
#post process, and calculate drms for each cluster
for i in `ls | egrep "clusters.*\.pdb$"` ; do
	./processTraj $i $7 $8 $9 $10 $11 $12 $13 $14 $15
done
