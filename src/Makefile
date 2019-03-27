# directories
BASEDIR=$(shell pwd)

# Mac
TPLDIR=$(BASEDIR)/libs/seq_mac_dbg_static/_build/

# linux
#TPLDIR=$(BASEDIR)/libs/seq_mac_dbg_static/_build/

## windows
#TPLDIR=$(BASEDIR)/libs/seq_win_dbg_static/_build/

# compilers
EXE_EXT = 
LIB_EXT = .a
OBJ_EXT = .o
CC = clang
CXX = g++

# flags debug
#CFLAGS = -std=c99 -g -O0 -Wall -pedantic -w
#-fdiagnostics-color=auto (entre c99 e -g)
# flags otimizadas
CFLAGS = -std=c99 -O2


# add solver libraries and includes
INCS	= -I${TPLDIR}/hdf5/include -I${TPLDIR}/cgns/include
LDFLAGS	= -L${TPLDIR}/hdf5/lib -L${TPLDIR}/cgns/lib
#INCS    = -I${TPLDIR}/cgns/include
#LDFLAGS = -L${TPLDIR}/cgns/lib
LIBS	= -lcgns -lhdf5_hl -lhdf5 -lm -ldl -llapack -lblas

all: removeFiles datproc

removeFiles:
	@rm -f datproc *.o

datproc:
	${CC} io_cgns_read.c io_cgns_write.c io_bru_write.c \
	io_native_read.c io_native_write.c io_gui.c \
	hdf5_funcs.c solution.c native.c utils.c \
	BSpline.c igesRead.c rbf_functions.c moveMesh.c main.c \
	${CFLAGS} ${INCS} ${LDFLAGS} ${LIBS} -o datproc${EXE_EXT}

clean:
	@echo cleaning ...
	@rm -f datproc *.o
