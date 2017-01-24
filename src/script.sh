#!/bin/bash -l

folder="/inputData/"
executable=biLouvain

datasets=(SouthernWomen)
for i in ${datasets[@];}
do
	fileName=$folder${i}"_bipartite.txt"
	./$executable -i $fileName -d "\t" -fuse 1
done
exit 0

