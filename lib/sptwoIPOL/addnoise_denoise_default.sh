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


input=$1
output=$2
sigma=$3

iframe=-1
if [ "$#" -gt 3 ]; then
    iframe=$4
fi

mkdir work
cd work

../add_noise_to_images.sh ../${input}/ 40

../src/sptwo -s ${sigma} -i ${iframe} noisyh.mov denoised.i
 
for i in denoised.i*.png; do cp $i ../${output}/.; done
for i in noisyh.i*; do cp $i ../${output}/.; done

echo "RMSE" > RMSE.txt
if [ "$#" -gt 3 ]; then
  ../src/computeRMSE -b 15 denoised.i${iframe}.png clean.i${iframe} >> RMSE.txt
else 
  i=0
  for input in `ls -v clean.i*`; do
     ../src/computeRMSE -b 15 denoised.i${i}.png clean.i${i} >> RMSE.txt
     (( i++ ))
  done
fi

cp RMSE.txt ../${output}/.

cd ..
rm -rf work



