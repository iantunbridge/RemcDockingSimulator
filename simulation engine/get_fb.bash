#!/bin/bash
for i in $1*Bound;
do
	
	a=`echo $i | tr "_" " " | egrep -o "[^ ]*microM" | tr -d "microM"`;
	b=`cat $i | tr -d "|" | awk 'BEGIN{} {if(NR>1000){sum+=$12;count++;} } END{print sum/count;}'`
	printf "%4d %0.6f\n" $a $b;
done;

