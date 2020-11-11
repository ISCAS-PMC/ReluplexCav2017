
#
# Local Variables:
# compile-command: "make -C .. "
# End:
#



#-----------------------------------
# C++ part
#-----------------------------------
# -L/usr/local/lib  -lapronxx -lt1pMPQ -lboxMPQ -loctMPQ  -lpolkaMPQ

FLAGS=  main.cpp -MMD -Wall -Wextra -Wno-deprecated -Wno-unused-but-set-variable -std=c++0x -O3 \
       -I../../common \
       -I../../glpk-4.60/installed/include \
       -I../../reluplex \
       -I../../nnet ./nnet.obj \
       -I../../abs \
       -I. \
       -lglpk -lapron -lapronxx -lgmp -lgmpxx -lmpfr -lm -lpolkaMPQ -loctD -lboxD -lt1pD \
       -L../../glpk-4.60/installed/lib\
       -L /usr/local/lib/

all: adversarialBOXSym.elf adversarialTOPO.elf adversarialPure.elf 
## adversarialTOPO.elf adversarialPure.elf 

adversarialBOXSym.elf: main.cpp nnet.obj
	g++  -D ABST=Abst_symV -D DOMT=BOX -D AION=1  $(FLAGS) -o $@

adversarialTOPO.elf: main.cpp nnet.obj
	g++  -D ABST=Abst_apron -D DOMT=TOPO -D AION=1  $(FLAGS) -o $@

adversarialPure.elf: main.cpp nnet.obj
	g++  -D ABST=Abst_symV -D DOMT=BOX -D AION=0  $(FLAGS) -o $@

nnet.obj: ../../nnet/nnet.cpp
	g++ -I ../../reluplex ../../common   -o nnet.obj -c ../../nnet/nnet.cpp

clean:
	rm adversarialBOXSym.elf adversarialTOPO.elf adversarialPure.elf 
