#! /bin/bash
PATH_PRE="$HOME/Dropbox/Project/python/Monte-Carlo-Simulation/results"
echo $PATH_PRE
FILE_NAME="mc-rect-lattice-func-linux4"
for i in 101 102 103 104 105 106 107 108
do
    #echo $PATH_PRE"$i"
    cp -v $FILE_NAME $PATH_PRE$FILENAME"$i"
done
