#!/bin/bash


for j in `seq 5670582 5670682` ; do
	scancel $j
	echo $j
done
