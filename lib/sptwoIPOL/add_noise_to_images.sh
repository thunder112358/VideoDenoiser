#!/bin/bash

#/*----------------------------------------------------------------------------
#
#Copyright (c) 2016-2018 A. Buades and J.L. Lisani
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU Affero General Public License for more details.
#
#You should have received a copy of the GNU Affero General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------*/


image=$1
sigma=$2

rm noisyh.i* noisyh.mov

frames=0
for input in `ls -v $image*`; do
    cp $input clean.i${frames}
    (( frames++ ))
done

#add noise
iframemax=$(( frames-1 ))
for i in $(seq 0 $iframemax); do 
  ../src/addGaussianNoise -a $sigma clean.i${i} noisyh.i${i}.tiff
done

# build sequence
echo "cflmovie" > noisyh.mov
echo "$frames" >> noisyh.mov
for i in $(seq 0 $iframemax); do echo noisyh.i${i}.tiff >> noisyh.mov; done

#convert noisy frames to PNG format
for i in $(seq 0 $iframemax); do 
  ../src/convert8bits -c 1 noisyh.i${i}.tiff noisyh.i${i}.png
done








