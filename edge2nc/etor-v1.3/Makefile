#MAKEFILE FOR ETOR (Edge to Radar Software Library)

RSL_INCDIR = .
EDGE=.
EDGE_INCDIR = $(EDGE)/.

CC=gcc -O3 -fomit-frame-pointer -fPIC
CFLAGS= -I$(RSL_INCDIR) -I$(EDGE_INCDIR)

LIB=libetor.a

# Install locations
SRC= \
	edge_to_radar.c \
	RSL_to_EDGE.c \
	get_path.c \
	save_edge_vol.c \
	curr_time.c \
	load_data.c \
	copy_buff.c \
	lzw.c \
	save_data.c \
	run_len_dec.c \
	crc.c \
	bitx.c \
	run_len_enc.c \
	seterr.c \
	culprit.c \
	error.c \
	serror.c \
	set_message_routines.c \
	err_id.c 

OBJ = $(SRC:.c=.o)

LIBFILES = $(OBJ)

libo = $(OBJ)
#--------------------------------------------------------------------
all: lib

lib: $(LIB)

$(LIB): $(OBJ)
	ar ruv $(LIB) $(OBJ)
	ranlib $(LIB)

install: $(LIB)
	cp $(LIB)    $(INSTALL_LIBDIR)
	cp antenna.h $(INSTALL_INCDIR)
	cp vol.h     $(INSTALL_INCDIR)
	cp radar.h   $(INSTALL_INCDIR)
	cp rcu.h     $(INSTALL_INCDIR)

$(OBJ): Makefile *.h

clean:
	rm -f $(LIBFILES) $(LIB) $(FILES)
