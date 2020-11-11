
#-----------------------------------
# C++ part
#-----------------------------------
# -L/usr/local/lib  -lapronxx -lt1pMPQ -lboxMPQ -loctMPQ  -lpolkaMPQ



all:absBoxD absBoxMPQ absT1pD absT1pMPQ absBoxDSym absBoxMPQSym absT1pDSym absT1pMPQSym

FLAGS = -std=c++11 -I ../ReluplexCav2017/nnet -I ../ReluplexCav2017/common -I ~/ljcLib/apron/apron/taylor1plus -U__STRICT_ANSI__ -UNDEBUG -O3 -g  -Wcast-qual -Wswitch -Wall -Wextra -Wundef -Wcast-align -Wno-unused -Wno-unused-parameter -Wno-unused-function -fPIC -Werror-implicit-function-declaration   -U__STRICT_ANSI__ -UNDEBUG  main.cpp -o abs ../ReluplexCav2017/nnet/nnet.o   -L/usr/local/lib   -lapronxx_debug 	  -lgmp -lgmpxx  -lapron_debug  -lmpfr -lm  -lpolkaMPQ

absBoxD: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D BDOMAIN -D ABST=Abst_apron -D NUM_D -loctD  -lboxD -lt1pD -o $@
absBoxMPQ: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D BDOMAIN -D ABST=Abst_apron -D NUM_MPQ -loctMPQ  -lboxMPQ -lt1pMPQ -o $@
absT1pD: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D TDOMAIN -D ABST=Abst_apron -D NUM_D -loctD  -lboxD -lt1pD -o $@
absT1pMPQ: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D TDOMAIN -D ABST=Abst_apron -D NUM_MPQ -loctMPQ  -lboxMPQ -lt1pMPQ -o $@
absBoxDSym: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D BDOMAIN -D ABST=Abst_symV -D NUM_D -loctD  -lboxD -lt1pD -o $@
absBoxMPQSym: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D BDOMAIN -D ABST=Abst_symV -D NUM_MPQ -loctMPQ  -lboxMPQ -lt1pMPQ -o $@
absT1pDSym: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D TDOMAIN -D ABST=Abst_symV -D NUM_D -loctD  -lboxD -lt1pD -o $@
absT1pMPQSym: main.cpp ../ReluplexCav2017/nnet/nnet.o
	g++  $(FLAGS) -D TDOMAIN -D ABST=Abst_symV -D NUM_MPQ -loctMPQ  -lboxMPQ -lt1pMPQ -o $@

../ReluplexCav2017/nnet/nnet.o: ../ReluplexCav2017/nnet/nnet.cpp
	g++ -I ../ReluplexCav2017/reluplex ../ReluplexCav2017/common   -o ../ReluplexCav2017/nnet/nnet.o -c ../ReluplexCav2017/nnet/nnet.cpp

clean:
	rm absBoxD absBoxDSym absBoxMPQ absBoxMPQSym absT1pD absT1pDSym absT1pMPQ absT1pMPQSym

# gcc -std=c++11 -U__STRICT_ANSI__ -DNDEBUG -O3  -Wcast-qual -Wswitch -Wall -Wextra -Wundef -Wcast-align -Wno-unused -Wno-unused-parameter -Wno-unused-function -fPIC -I /usr/local/lib -I /usr/local/include -I ~/ljcLib/apron/apron/taylor1plus  -I ../nnet -I ../common main.cpp -o abs -lpolkaMPQ_debug  -loctMPQ_debug -L../box -lboxMPQ_debug  -lt1pMPQ_debug /usr/local/lib/libapronxx_debug.a
