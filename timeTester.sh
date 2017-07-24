#!/bin/bash

#controllo se input sia un numero
re='^[0-9]+$'
if [[ $1 =~ $re ]] ; then
	nTest=$1 
else
	nTest=50
	echo l\'input non è un numero
fi

#se non è un numero faccio 50 test
echo "il programma verra eseguito $nTest volte"

tot=0

for i in $(seq 1 $nTest)
do
	t=$(./gm --test < input.txt | tail -n1)
	tot=$(echo $tot+$t | bc -l)
done

echo $(echo $tot/$nTest | bc -l)
