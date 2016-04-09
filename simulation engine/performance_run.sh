#!/bin/bash
echo 100 $1 $2 
./$1 -f performance_configs/UbUIM  -t $3 -g $3 -bx $2 > performance_output/$1_$2_0100
echo 402 $1 $2
./$1 -f performance_configs/2pcc  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_0402
echo 584 $1 $2
./$1 -f performance_configs/2g33  -t 1 -r 1 -e 1 -m 1000 -bx $2 >  performance_output/$1_$2_0584
echo 1136 $1 $2
./$1 -f performance_configs/2g33x4  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_1136
echo 1704 $1 $2
./$1 -f performance_configs/2g33x6  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_1704
echo 2272 $1 $2
./$1 -f performance_configs/2g33x8  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_2272
echo 3408 $1 $2
./$1 -f performance_configs/2g33x12  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_3408
echo 4544 $1 $2
./$1 -f performance_configs/2g33x16  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_4544
echo 5680 $1 $2
./$1 -f performance_configs/2g33x20  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_5680
echo 6818 $1 $2
./$1 -f performance_configs/2g33x24  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_6816
echo 7668 $1 $2
./$1 -f performance_configs/2g33x27  -t 1 -r 1 -e 1 -m 1000 -bx $2 > performance_output/$1_$2_7668

