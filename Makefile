# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
FC = nagfor
#  FC = f90             
#  FC = gfortran      
#  FC = ifort
endif
#
#option to choose level of optimization:  make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
  else                           
#   FFLAGS =  -q  -C          # debug=1  (basic gtortran debug option)
#    FFLAGS =  -g     -u       # debug=1  (basic debug option)
    ifeq ($(debug),2)
        FFLAGS =  -C -g        # debug=2  (higher-level debug option)
    endif
#                 can add more options for  debug=3 , etc., as desired
  endif
#
# as usual, list the objects
#
OBJECTS = parfit.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o parfit.x
# To run the code, execute:
# ./level.x < input.5 > fort.6


