ARMA_PATH=/home/anandps/local/armadillo
CXX=g++


all:
	${CXX} -c main.cpp -I${ARMA_PATH}/include -L${ARMA_PATH}/lib -std=c++14 -larmadillo -g -I../src -L../src
	${CXX} main.o ../src/main.o -o main -I${ARMA_PATH}/include -L${ARMA_PATH}/lib -std=c++14 -larmadillo -g -I../src -L../src 
gems:
	${CXX} main_gems.cpp -o main_gems -I${ARMA_PATH}/include -L${ARMA_PATH}/lib -std=c++14 -larmadillo -g -I../src -L../src -lmain 

swig:
	swig -python -c++ main.i
	${CXX} -c main_wrap.cxx  -I${ARMA_PATH}/include -std=c++14  -I${PYTHON_INC} -fPIC
	ld -shared main_wrap.o -o _main.so -larmadillo  -L${ARMA_PATH}/lib 
