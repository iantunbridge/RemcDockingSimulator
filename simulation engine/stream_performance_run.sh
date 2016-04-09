#!/bin/bash
echo `date  "+%H:%M:%S"` "Testing kernel: $1 $2 streams $3 Threads"
echo `date  "+%H:%M:%S"` "100"
./$1 -f performance_configs/UbUIM -t $3 -s $2 > performance_output/0100_UIM_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "402"
./$1 -f performance_configs/2pcc -t $3 -s $2 > performance_output/0402_2pcc_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "584"
./$1 -f performance_configs/2g33 -t $3 -s $2  > performance_output/0584_2g33_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "1136"
./$1 -f performance_configs/2g33x4 -t $3 -s $2  > performance_output/1136_2g33x4_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "1704"
./$1 -f performance_configs/2g33x6 -t $3 -s $2 > performance_output/1704_2g33x6_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "2272"
./$1 -f performance_configs/2g33x8 -t $3 -s $2  > performance_output/2272_2g33x8_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "3408"
./$1 -f performance_configs/2g33x12 -t $3 -s $2 0 > performance_output/3408_2g33x12_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "4544"
./$1 -f performance_configs/2g33x16 -t $3 -s $2 > performance_output/4544_2g33x16_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "5680"
./$1 -f performance_configs/2g33x20 -t $3 -s $2 > performance_output/5680_2g33x20_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "6816"
./$1 -f performance_configs/2g33x24 -t $3 -s $2  > performance_output/6816_2g33x24_$3_threads_streams$2
echo `date  "+%H:%M:%S"` "7668"
./$1 -f performance_configs/2g33x27 -t $3 -s $2  > performance_output/7668_2g33x27_$3_threads_streams$2
