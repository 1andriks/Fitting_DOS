#variables
NANOROOT=../
CXX=mpicxx #compilatore parallelo
#mpixlcxx
CXXSR=g++ #compilatore seriale
#xlc++
FLAGS=-O3 -ffast-math
#
#
#-O3 -qmaxmem=-1 
#-std=c++0x
#-ftrapuv -strict-ansi -O0 -Wcheck
#-O1 -xHost -finline #ottimizzazioni
#-O3 -ffast-math
#-ansi -O0 -Wall 
#
EXEPROG=nano
EXEPOST=post
INCL=-I. -I$(NANOROOT)/src  -I$(NANOROOT)/src/lapacke  
FFTWDIR=/usr/local/fftw-3.3.3
INCL_FFTW=-I$(FFTWDIR)/include 
LIBFFTW=#$(FFTWDIR)/lib/libfftw3.a 
LIBLAP=liblapacke.a liblapack.a blas_LINUX.a -L/usr/local/opt/gfortran/gfortran/lib/ -lgfortran
INFILE=input.in
SRC=$(NANOROOT)/src
OBJECTS= line_parser.o p_rand.o lattices.o myfile.o	
UTILS= line_parser_s.o build_input.o	
POST=post.o 	

.SUFFIXES:.cc .o
.cc.o:
	$(CXX) $(INCL) $(FLAGS) -c  $<

code:
	make utils
	make build
	./build $(INFILE)
	make body	
	make nano

post: $(SRC)/post.cc
	$(CXX) $(INCL) $(INCL_FFTW) $(FLAGS) $? -o $(EXEPOST) $(LIBFFTW)	

body: $(OBJECTS)

utils: $(UTILS)

line_parser_s.o: $(SRC)/line_parser.cc 
	$(CXXSR) $(INCL) $(FLAGS) -c  $? -o line_parser_s.o

line_parser.o: $(SRC)/line_parser.cc
	$(CXX) $(INCL) $(FLAGS) -c $?

lattices.o: $(SRC)/lattices.cc
	$(CXX) $(INCL) $(FLAGS) -c  $?

p_rand.o: $(SRC)/p_rand.cc
	$(CXX) $(INCL) $(FLAGS) -c  $?

myfile.o: myfile.cc 
	$(CXX) $(INCL) $(INCL_FFTW) $(FLAGS) -c  $? 

build_input.o: $(SRC)/build_input.cc
	$(CXXSR) $(INCL) $(FLAGS) -c  $?

nano: $(OBJECTS)
	$(CXX) $(INCL) $(INCL_FFTW) $(FLAGS) $(OBJECTS) -o $(EXEPROG) $(LIBFFTW) $(LIBLAP)

build: line_parser_s.o build_input.o
	$(CXXSR) $(INCL) $(FLAGS) line_parser_s.o build_input.o -o ./build

# clean compiled files
clean:
	rm -f *.o 
cl:
	rm nano* myfile.o
