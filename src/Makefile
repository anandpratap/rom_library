ARMA_PATH=/home/anandps/local/armadillo
PYTHON_INC=/home/anandps/local/anaconda2/include/python2.7
CXX=g++


all:
	${CXX} main.cpp -o main -I${ARMA_PATH}/include -L${ARMA_PATH}/lib -std=c++14 -larmadillo -g
swig:
	swig -python -c++ main.i
	${CXX} -c main_wrap.cxx  -I${ARMA_PATH}/include -std=c++14  -I${PYTHON_INC} -fPIC
	ld -shared main_wrap.o -o _main.so -larmadillo  -L${ARMA_PATH}/lib 
