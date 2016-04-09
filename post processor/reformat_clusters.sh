#!/bin/bash
for i in `ls | egrep "clusters.pdb$"` ; do
	echo "processTraj $i processTraj.exp $@ "
	./processTraj $i processTraj.exp $@ 
done
