#!/bin/bash
sed s/END/ENDMDL/ $1 > $1.processed.pdb
mv $1.processed.pdb $1
