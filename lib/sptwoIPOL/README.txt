VIDEO DENOISING WITH OPTICAL FLOW ESTIMATION
--------------------------------------------

Version 2.0 - March 23, 2018
by Antoni Buades <toni.buades@uib.es> 
   and Jose-Luis Lisani <joseluis.lisani@uib.es>


ASSOCIATED PUBLICATION

This code is associated to the following IPOL publication:

A. Buades, J.L. Lisani, M. Miladinovic,
Video Denoising with Optical Flow Estimation, 
Image Processing On Line

An online demo that uses the code is available at:


ORGANIZATION OF THE CODE

The code is organized as follows:

1) 'src' folder: contains the code for denoising a set of frames. 
The folder contains two subfolders:
	1.1) libNoiseVideo: contains the implementation of the main algorithms
of the denoising method
	1.2) library: contains auxiliary functions that manage input/output of 
images 

Four functions are provided:
- sptwo.cpp: implements the denoising method
- computeRMSE.cpp: computes the Root Mean Squared Error between two images
- addGaussianNoise.cpp: adds white Gaussian noise to an image
- convert8bits.cpp: converts float image data to unsigned char (PNG format)

The code described in the associated publication and that has been peer-reviewed
is the one contained in the libNoiseVideo subfolder and in sptwo.cpp.

2) 'data' folder: contains short sequences to test the denoising method.
The folder contains three subfolders:
	2.1) cleanseq: four frames (292 x 194) of a noise-free sequence (PNG 
format)
	2.2) noisyseq: four noisy frames (the result of adding a Gaussian 
noise of standard deviation 40 to the clean frames). The format of the frames 
is TIFF (float data)
	2.3) noisyseq(png): the four noisy frames from folder noiseseq 
converted to 8-bit format (PNG). Remark that the values are been rounded and 
clipped to the [0, 255] range.

3) 'test' folder: empty folder where the results of the denoising are stored

4) scripts: four scripts are provided to help the user apply the denoising 
method to a set of frames. Some examples of usage are presented below.
- add_noise_to_images.sh: adds noise to a set of frames in a given folder
- convert_images_to_mov.sh: creates a text file containing the list of frames 
to be denoised
- denoise_default.sh: denoises a set of frames using the defaults parameters 
described in the associated paper
- addnoise_denoise_default.sh: adds noise to a clean sequence and then 
denoises it using the defaults parameters described in the associated paper



COMPILATION

1) Decompress code and access to folder:
tar xvzf sptwoIPOL.tgz
cd sptwoIPOL

2) Compilation:
cd src
make OMP=1 #compilation using OMP directives, to compile without them use: make
cd ..

The executable files obtained after compilation are 'sptwo', 'computeRMSE', 
'addGaussianNoise' and 'convert8bits'



USAGE

sptwo [-s s] [-i i] [-b b] [-w w] [-t t] [-k k] [-f f] [-h h] [-o o] [-p p] 
      [-q q] [-d d] [-e e] [-l l] [-m m]  input  out 
	-s  s	 noise standard deviation (Default: 5.0)
	-i  i	 frame to denoise (-1: denoise all frames) (Default: -1)
	-b  b	 radius of search region (Default: 12)
	-w  w	 radius of patch (Default: 2)
	-t  t	 radius of temporal neighborhood (Default: 7)
	-k  k	 minimum number of patches (recommended: 55 gray images, 
95 color images) 
	-f  f	 flat parameter (Default: 0.85f)
	-h  h	 occlusion binarization threshold (Default: 0.5f)
	-o  o	 occlusion factor (Default: 5.5f)
	-p  p	 PCA factor 1st step (Default: 1.8f)
	-q  q	 PCA factor 2nd step (Default: 1.45f)
	-d  d	 3D blocks distances 1st step (Default: 0.0f)
	-e  e	 3D blocks distances 2nd step (Default: 2.0f)
	-l  l	 optical flow lambda 1st step (Default: 0.075f)
	-m  m	 optical flow lambda 2nd step (Default: 0.15f)
	input	 input file
	out	 output file

The different usage options are described in the paper accompanying the code.


EXAMPLES OF USE

1) Denoise given sequence with default parameters:

- we are assuming that the frames of the sequence to denoise are 
in sptwoIPOL/data/noisyseq, and that the denoising result will be written in 
folder sptwoIPOL/test/
- we need to provide as input a level of noise, e.g. 40

./denoise_default.sh ./data/noisyseq ./test 40

- The denoised frames (in PNG format) are stored in ./test

Note: it is possible to denoise just a single frame, e.g.
to denoise the frame with index 2 (recall that the first index is 0)
./denoise_default.sh ./data/noisyseq ./test 40 2



2) Add noise to a clean sequence and then denoise, using default parameters:

- we are assuming that the frames of the clean sequence are 
in sptwoIPOL/data/cleanseq, and that the denoising result will be written in 
folder sptwoIPOL/test/
- we need to provide as input a level of noise, in this example is 40

./addnoise_denoise_default.sh ./data/cleanseq ./test 40

- The denoised frames (in PNG format) are stored in ./test
- The noisy frames (in TIFF format) are stored in ./test
- An 8-bit version of the noisy frames (in PNG format) is also stored 
in ./test
- Remark that the 8-bit version is obtained after rounding and clipping the 
float values to [0, 255]
- The RMSE values per frame are in file ./test/RMSE.txt

Note: it is possible to denoise just a single frame, e.g.
to denoise the frame with index 2 (recall that the first index is 0)
./addnoise_denoise_default.sh ./data/cleanseq ./test 40 2



3) Denoise given sequence with parameters different from the default values:

- We are assuming that the frames of the sequence to denoise are 
in sptwoIPOL/data/noisyseq

- make 'work' folder:
mkdir work
cd work

- create text file (noisyh.mov) with the names of the frames:
../convert_images_to_mov.sh ../data/noisyseq/ 

- denoise sequence:
(we need to provide as input the level of noise of the image, e.g. 40)
In this example we use as parameters for denoising Knn=60, radius of patch=3
and a PCA factor for the 1st step=2.0 

../src/sptwo -s 40 -k 60 -w 3 -p 2.0 noisyh.mov denoised.i

The output is a set of denoised frames (in PNG format) with names of the form: 
denoised.i0.png, denoised.i1.png, denoised.i2.png, etc



4) Add noise to a clean sequence and then denoise:

- We are assuming that the frames of the clean sequence are in 
sptwoIPOL/data/cleanseq

- make 'work' folder:
mkdir work
cd work

- add noise to clean sequence and create text file (noisyh.mov) with the 
names of the noisy frames:
(we need to provide as input the level of noise of the image, e.g. 40)

../add_noise_to_images.sh ../data/cleanseq/ 40

- denoise sequence:
(we need to provide as input the level of noise of the image, e.g. 40)
In this example we use as parameters for denoising Knn=60, radius of patch=3
and a PCA factor for the 1st step=2.0 

../src/sptwo -s 40 -k 60 -w 3 -p 2.0 noisyh.mov denoised.i

The output is a set of denoised frames with names of the form: 
denoised.i0.png, denoised.i1.png, denoised.i2.png, etc

- compute denoising error (RMSE), for one of the frames (e.g. frame 2):

../src/computeRMSE -b 15 denoised.i2.png clean.i2 

(the parameter -b 15 excludes from the computation the pixels close to the 
border of the image, at distance less than 15 pixels)



COPYRIGHT AND LICENSE

Copyright (c) 2016-2018 A. Buades and J.L. Lisani

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


THANKS

We would be grateful to receive any comment, especially about errors, bugs,
or strange results.


