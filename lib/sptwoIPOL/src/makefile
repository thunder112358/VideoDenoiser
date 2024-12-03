#/*---------------------------------------------------------------------------
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
#---------------------------------------------------------------------------*/



# C source code
CSRC	= library/libImageFormats.c
# C++ source code
CXXSRC	= library/libImageFormatPM.cpp library/libImage.cpp library/libBasic.cpp \
          libNoiseVideo/libFlow.cpp libNoiseVideo/libNLPCA.cpp libNoiseVideo/libDenoisingVideo.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= sptwo addGaussianNoise computeRMSE convert8bits

default	: $(BIN)

# use DEBUG mode with `make DEBUG=1`
ifdef DEBUG
  # C optimization flags
  COPT = -g
else
  # C optimization flags
  COPT	= -O3 -ftree-vectorize -funroll-loops -fomit-frame-pointer -ffast-math 
endif


# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -ansi -std=c99
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -ansi -pedantic

# link flags
LDFLAGS	= -lpng -lm -lfftw3 -lfftw3f -ltiff  -ljpeg 



# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp 
endif



# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS) 
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS) 
    
# link all the object code
$(BIN): % : %.o $(OBJ) 
	$(CXX) -o $@ $^ $(LDFLAGS) 

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN) *.o *~

