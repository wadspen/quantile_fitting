#!/bin/bash


for j in `seq 5713939 5714944` ; do
	scancel $j
	echo $j
done
