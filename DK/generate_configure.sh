#!/bin/bash

root_directory=`pwd`
echo $root_directory
protein_lst=$(find $root_directory/shj/superimposed/ -name *.pdbqt)
for i in $protein_lst; do
    if [ ! -d $root_directory/$(basename $i .pdbqt) ]; then
        mkdir $root_directory/$(basename $i .pdbqt)
    fi
done
    cat <<config > config_$(basename $i .pdbqt).txt
receptor = $i
ligand_directory = $root_directory/shj/filter
output_directory = $root_directory/$(basename $i .pdbqt)
opencl_binary_path = /home/kwon/Desktop/vina/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1
center_x = -0.188
center_y = 17.878
center_z = 43.341
size_x = 15
size_y = 15
size_z = 15
thread = 1000
config
done

for i in `ls | grep .txt`; do
    echo "Autodock start job name : $i"
    ./AutoDock-Vina-GPU-2-1  --config ./$i > "${i%.txt}.log"
    echo "Autodock finish job name : $i"
done