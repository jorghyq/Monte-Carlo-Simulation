#! /bin/bash
PATH_PRE="$HOME/Dropbox/Project/python/Monte-Carlo-Simulation/results"
echo $PATH_PRE
FILE_NAME="mc-rect-lattice-func-linux6"
FILE_NAME1="call_c++_single.py"
FILE_NAME2="call_c++_n2.py"
for i in 70 71 72 73 74
do
    #echo $PATH_PRE"$i"
    cp -v $FILE_NAME $PATH_PRE"$i"
    cp -v $FILE_NAME1 $PATH_PRE"$i"
    cp -v $FILE_NAME2 $PATH_PRE"$i"
done
